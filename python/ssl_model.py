
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_inputs(root):
    manifest = pd.read_csv(root / "export_manifest.tsv", sep="\t")
    for row in manifest.itertuples(index=False):
        assert md5(root / row.filename) == row.md5, row.filename

    shapes = pd.read_csv(root / "shapes.tsv", sep="\t")
    genes = pd.read_csv(root / "union_genes.tsv", sep="\t")
    hierarchy = pd.read_csv(root / "token_hierarchy.tsv", sep="\t")
    members = pd.read_csv(root / "token_gene_membership.tsv", sep="\t")

    assert len(genes) == 1753
    assert genes.ensembl_gene_id.is_unique
    assert len(hierarchy) == 37
    assert hierarchy.token_name.is_unique
    assert np.array_equal(
        genes.gene_index0.to_numpy(), np.arange(len(genes))
    )
    assert np.array_equal(
        hierarchy.token_index0.to_numpy(), np.arange(len(hierarchy))
    )

    concepts = sorted(hierarchy.concept_id.unique().tolist())
    modules = sorted(hierarchy.module_id.unique().tolist())
    assert len(concepts) == 22 and len(modules) == 8

    membership = np.zeros((37, len(genes)), dtype=bool)

    for row in members.itertuples(index=False):
        t, g = int(row.token_index0), int(row.gene_index0)
        assert 0 <= t < 37 and 0 <= g < len(genes)
        assert hierarchy.iloc[t].token_name == row.token_name
        assert genes.iloc[g].ensembl_gene_id == row.ensembl_gene_id
        assert not membership[t, g]
        membership[t, g] = True

    assert membership.any(axis=1).all()
    assert membership.any(axis=0).all()

    tc = np.zeros((37, 22), dtype=np.float32)
    cm = np.zeros((22, 8), dtype=np.float32)

    for i, row in hierarchy.iterrows():
        tc[i, concepts.index(row.concept_id)] = 1.0

    for i, concept in enumerate(concepts):
        parents = hierarchy.loc[
            hierarchy.concept_id == concept, "module_id"
        ].unique()
        assert len(parents) == 1
        cm[i, modules.index(parents[0])] = 1.0

    assert (tc.sum(axis=0) > 0).all()
    assert (cm.sum(axis=0) > 0).all()

    tc /= tc.sum(axis=0, keepdims=True)
    cm /= cm.sum(axis=0, keepdims=True)

    arrays = {}
    ssl_metadata = None

    for row in shapes.itertuples(index=False):
        name = row.dataset
        n, g, t = int(row.n_samples), int(row.n_SSL_genes), int(row.n_tokens)
        assert g == len(genes) and t == 37

        x = np.fromfile(
            root / f"{name}_gene_z.float32.bin", dtype="<f4"
        )
        rank = np.fromfile(
            root / f"{name}_token_rank.float32.bin", dtype="<f4"
        )
        observed = np.fromfile(
            root / f"{name}_gene_observed.uint8.bin", dtype=np.uint8
        )

        assert x.size == n*g and rank.size == n*t
        assert observed.size == n*g
        assert np.isin(observed, [0, 1]).all()

        x = x.reshape(n, g)
        rank = rank.reshape(n, t)
        observed = observed.reshape(n, g).astype(bool)

        assert np.isfinite(x).all() and np.isfinite(rank).all()
        assert np.max(np.abs(x)) <= 6.00001
        assert rank.min() >= 0 and rank.max() <= 1.00001

        if name == "ssl":
            metadata = pd.read_csv(root / "ssl_metadata.tsv", sep="\t")
            ssl_metadata = metadata
            assert observed.all()
            assert (metadata["split"] == "train").sum() == 554
            assert (metadata["split"] == "validation").sum() == 62
        else:
            # Deliberately do not read response values.
            metadata = pd.read_csv(
                root / f"{name}_metadata.tsv", sep="\t",
                usecols=["patient_uid"]
            )

        assert len(metadata) == n and metadata.patient_uid.is_unique
        arrays[name] = {
            "x": x, "rank": rank, "observed": observed,
            "ids": metadata.patient_uid.to_numpy()
        }

    ssl_ids = set(arrays["ssl"]["ids"])
    for name, values in arrays.items():
        if name != "ssl":
            assert not ssl_ids.intersection(values["ids"])

    checks = pd.read_csv(root / "R_Python_value_checks.tsv", sep="\t")
    for row in checks.itertuples(index=False):
        key = "x" if row.kind == "gene_z" else "rank"
        actual = arrays[row.dataset][key][
            int(row.sample_index1)-1, int(row.feature_index1)-1
        ]
        assert np.isclose(actual, row.value, rtol=1e-6, atol=1e-6)

    spec = {
        "gene_ids": genes.ensembl_gene_id.tolist(),
        "token_names": hierarchy.token_name.tolist(),
        "concept_names": concepts,
        "module_names": modules,
        "n_genes": len(genes),
        "n_tokens": 37,
        "n_concepts": 22,
        "n_modules": 8,
        "dim": 24,
        "heads": 4,
        "architecture": "member_masked_gene_to_token_reconstruction_v1"
    }

    return arrays, ssl_metadata, membership, tc, cm, spec


class SSL37(nn.Module):
    def __init__(self, membership, token_to_concept, concept_to_module,
                 dim=24, heads=4):
        super().__init__()
        n_tokens, n_genes = membership.shape
        assert dim % heads == 0

        self.dim = dim
        self.n_genes = n_genes
        self.n_tokens = n_tokens

        self.register_buffer(
            "membership", torch.as_tensor(membership, dtype=torch.bool)
        )
        self.register_buffer(
            "token_to_concept",
            torch.as_tensor(token_to_concept, dtype=torch.float32)
        )
        self.register_buffer(
            "concept_to_module",
            torch.as_tensor(concept_to_module, dtype=torch.float32)
        )

        # For each gene, average information over its member tokens only.
        decoder_membership = membership.astype(np.float32)
        decoder_membership /= decoder_membership.sum(
            axis=0, keepdims=True
        )
        self.register_buffer(
            "decoder_membership",
            torch.as_tensor(decoder_membership)
        )

        self.gene_id = nn.Parameter(torch.randn(n_genes, dim) * 0.03)
        self.gene_value = nn.Linear(1, dim, bias=False)
        self.missing_embedding = nn.Parameter(torch.zeros(dim))
        self.gene_norm = nn.LayerNorm(dim)

        self.token_query = nn.Parameter(torch.randn(n_tokens, dim) * 0.04)
        self.token_id = nn.Parameter(torch.randn(n_tokens, dim) * 0.03)
        self.token_norm = nn.LayerNorm(dim)
        self.token_ff = nn.Sequential(
            nn.Linear(dim, dim), nn.GELU(), nn.Linear(dim, dim)
        )

        self.rank_decoder = nn.Linear(dim, 1)
        self.gene_decoder = nn.Linear(dim, dim)
        self.reconstruction_gene_id = nn.Parameter(
            torch.randn(n_genes, dim) * 0.03
        )
        self.reconstruction_bias = nn.Parameter(torch.zeros(n_genes))

        self.module_attention = nn.MultiheadAttention(
            dim, heads, batch_first=True, dropout=0.1
        )
        self.module_norm = nn.LayerNorm(dim)

    def encode(self, x, observed):
        assert x.ndim == 2 and x.shape[1] == self.n_genes
        assert x.shape == observed.shape

        h = (
            self.gene_id[None]
            + self.gene_value(x.unsqueeze(-1))
            + (~observed).unsqueeze(-1).to(x.dtype)
            * self.missing_embedding
        )
        h = self.gene_norm(h)

        logits = torch.einsum(
            "td,bgd->btg", self.token_query, h
        ) / (self.dim ** 0.5)

        # Nonmembers never contribute to a token.
        logits = logits.masked_fill(
            ~self.membership[None], float("-inf")
        )
        weights = logits.softmax(dim=-1)

        tokens = self.token_norm(
            torch.bmm(weights, h) + self.token_id[None]
        )
        tokens = tokens + self.token_ff(tokens)

        concepts = torch.einsum(
            "tc,btd->bcd", self.token_to_concept, tokens
        )
        modules = torch.einsum(
            "cm,bcd->bmd", self.concept_to_module, concepts
        )
        context, _ = self.module_attention(
            modules, modules, modules, need_weights=False
        )
        module_context = self.module_norm(modules + context)

        # Exported tokens are local, before module interaction.
        return tokens, concepts, module_context

    def forward(self, x, observed):
        tokens, concepts, modules = self.encode(x, observed)

        # Efficient member-restricted gene reconstruction.
        decoded_tokens = self.gene_decoder(tokens)
        gene_context = torch.einsum(
            "btd,tg->bgd", decoded_tokens, self.decoder_membership
        )
        gene_prediction = (
            gene_context * self.reconstruction_gene_id[None]
        ).sum(dim=-1) / (self.dim ** 0.5)
        gene_prediction = gene_prediction + self.reconstruction_bias

        rank_prediction = torch.sigmoid(
            self.rank_decoder(tokens).squeeze(-1)
        )

        return {
            "gene_prediction": gene_prediction,
            "rank_prediction": rank_prediction,
            "tokens": tokens,
            "concepts": concepts,
            "modules": modules,
        }


def make_view(x, observed, generator, rate=0.15):
    artificially_masked = (
        torch.rand(x.shape, generator=generator, device=x.device) < rate
    ) & observed

    view_observed = observed & ~artificially_masked
    corrupted = x.masked_fill(artificially_masked, 0.0)
    return corrupted, view_observed, artificially_masked


def ssl_loss(model, x, rank_target, observed, generator):
    first_x, first_observed, reconstruction_mask = make_view(
        x, observed, generator
    )
    second_x, second_observed, _ = make_view(
        x, observed, generator
    )

    first = model(first_x, first_observed)
    second_tokens, second_concepts, second_modules = model.encode(
        second_x, second_observed
    )

    assert reconstruction_mask.any()

    gene_loss = F.mse_loss(
        first["gene_prediction"][reconstruction_mask],
        x[reconstruction_mask]
    )
    rank_loss = F.mse_loss(
        first["rank_prediction"], rank_target
    )
    consistency = F.mse_loss(first["modules"], second_modules)

    total = gene_loss + 0.5 * rank_loss + 0.05 * consistency

    return total, {
        "total": total,
        "gene_reconstruction": gene_loss,
        "rank_reconstruction": rank_loss,
        "module_consistency": consistency,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    torch.set_num_threads(min(8, torch.get_num_threads()))
    torch.manual_seed(20261116)
    np.random.seed(20261116)

    arrays, metadata, membership, tc, cm, spec = load_inputs(args.data)
    print("INPUT_AND_INDEX_CHECK: PASS", flush=True)

    model = SSL37(membership, tc, cm, dim=24, heads=4)
    args.output.mkdir(parents=True, exist_ok=True)
    init_path = args.output / "SSL37_UNTRAINED_init.pt"

    if init_path.exists():
        raise RuntimeError("Initialization already exists; use a new directory")

    # Save before testing. No optimizer is used in this script.
    torch.save({
        "state_dict": model.state_dict(),
        "spec": spec,
        "seed": 20261116,
        "optimizer_steps": 0,
        "data_manifest_md5": md5(args.data / "export_manifest.tsv"),
    }, init_path)

    train_rows = np.flatnonzero(
        metadata["split"].to_numpy() == "train"
    )[:4]

    x = torch.from_numpy(arrays["ssl"]["x"][train_rows].copy())
    rank_target = torch.from_numpy(
        arrays["ssl"]["rank"][train_rows].copy()
    )
    observed = torch.ones_like(x, dtype=torch.bool)

    model.eval()
    with torch.no_grad():
        result = model(x, observed)
        for key, values in result.items():
            assert torch.isfinite(values).all().item()
            print("OUTPUT:", key, tuple(values.shape), flush=True)

        alone = model.encode(x[:1], observed[:1])[0]
        assert torch.allclose(
            alone, result["tokens"][:1], atol=1e-5, rtol=1e-5
        )

        # A token must be unaffected by perturbing all its nonmember genes.
        token_index = 0
        altered = x[:1].clone()
        outside = ~model.membership[token_index]
        altered[:, outside] = altered[:, outside] + 1.0

        altered_token = model.encode(
            altered, observed[:1]
        )[0][:, token_index]

        assert torch.allclose(
            altered_token, result["tokens"][:1, token_index],
            atol=1e-5, rtol=1e-5
        )

    print("PATIENT_BATCH_INDEPENDENCE: PASS", flush=True)
    print("TOKEN_MEMBER_LOCALITY: PASS", flush=True)

    model.train()
    timings = []
    loss_values = None

    for repeat in range(2):
        model.zero_grad(set_to_none=True)
        generator = torch.Generator().manual_seed(1000 + repeat)
        start = time.perf_counter()

        loss, components = ssl_loss(
            model, x, rank_target, observed, generator
        )
        assert torch.isfinite(loss).item()
        loss.backward()

        for name, parameter in model.named_parameters():
            assert parameter.grad is not None, name
            assert torch.isfinite(parameter.grad).all().item(), name

        timings.append(time.perf_counter() - start)
        loss_values = {
            name: float(value.detach())
            for name, value in components.items()
        }

        print("CPU_SSL_BATCH:", repeat + 1,
              "seconds=", round(timings[-1], 4),
              "losses=", json.dumps(loss_values), flush=True)

    # Explicitly record important gradient groups.
    gradient_groups = {
        "gene_value": model.gene_value.weight,
        "token_query": model.token_query,
        "rank_decoder": model.rank_decoder.weight,
        "gene_decoder": model.gene_decoder.weight,
        "module_attention": model.module_attention.in_proj_weight,
    }
    gradient_norms = {
        name: float(parameter.grad.norm())
        for name, parameter in gradient_groups.items()
    }
    assert all(np.isfinite(v) and v > 0 for v in gradient_norms.values())

    report = {
        "device": "cpu",
        "batch_size": len(train_rows),
        "n_genes": spec["n_genes"],
        "n_tokens": 37,
        "n_concepts": 22,
        "n_modules": 8,
        "parameters": sum(p.numel() for p in model.parameters()),
        "CPU_batch_forward_backward_seconds": timings,
        "loss_components": loss_values,
        "gradient_norms": gradient_norms,
        "optimizer_steps": 0,
        "initialization_md5": md5(init_path),
        "data_manifest_md5": md5(args.data / "export_manifest.tsv"),
        "model_script_md5": md5(Path(__file__)),
        "response_labels_used": False,
    }

    (args.output / "SSL37_model_spec.json").write_text(
        json.dumps(spec, indent=2), encoding="utf-8"
    )
    (args.output / "SSL37_smoke_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )

    print("PARAMETERS:", report["parameters"], flush=True)
    print("GRADIENT_NORMS:", json.dumps(gradient_norms), flush=True)
    print("SSL37_MODEL_CHECK: PASS", flush=True)


if __name__ == "__main__":
    main()

