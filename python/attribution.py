
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_manifest(directory, filename):
    table = pd.read_csv(directory / filename, sep="\t")
    for row in table.itertuples(index=False):
        path = directory / row.filename
        if not path.exists() or md5(path) != row.md5:
            raise RuntimeError(f"Fingerprint mismatch: {path}")


def integrate(model, rank, context, coef, bias, cap, steps):
    nodes, weights = np.polynomial.legendre.leggauss(steps)
    nodes = (nodes + 1) / 2
    weights = weights / 2

    ridge_integral = torch.zeros_like(rank)
    neural_rank_integral = torch.zeros_like(rank)
    context_integral = torch.zeros_like(context)

    for alpha, weight in zip(nodes, weights):
        r = (rank * alpha).detach().requires_grad_(True)
        c = (context * alpha).detach().requires_grad_(True)

        logit = r @ coef + bias + cap * torch.tanh(model(r, c))
        probability = torch.sigmoid(logit)

        grad_r, grad_c = torch.autograd.grad(
            probability.sum(), (r, c)
        )

        # Exact gradient contribution through the Ridge branch.
        ridge_gradient = (
            probability.detach() *
            (1 - probability.detach())
        )[:, None] * coef[None, :]

        ridge_integral += weight * ridge_gradient
        neural_rank_integral += weight * (
            grad_r - ridge_gradient
        )
        context_integral += weight * grad_c

    ridge_ig = (rank * ridge_integral).detach().numpy()
    neural_rank_ig = (
        rank * neural_rank_integral
    ).detach().numpy()
    context_ig = (
        context * context_integral
    ).sum(dim=-1).detach().numpy()

    return ridge_ig, neural_rank_ig, context_ig


def summarize(long_table, level):
    feature_column = {
        "token": "token_name",
        "concept": "concept_id",
        "module": "module_id",
    }[level]

    identity_columns = {
        "token": [
            "token_name", "concept_id", "concept_name",
            "module_id", "module_name",
        ],
        "concept": [
            "concept_id", "concept_name",
            "module_id", "module_name",
        ],
        "module": ["module_id", "module_name"],
    }[level]

    records = []

    def summarize_group(group, scope, cohort):
        row = {
            column: group.iloc[0][column]
            for column in identity_columns
        }

        row.update({
            "scope": scope,
            "cohort": cohort,
            "n_patients": len(group),
            "mean_signed_contribution":
                float(group.total.mean()),
            "median_signed_contribution":
                float(group.total.median()),
            "mean_absolute_contribution":
                float(group.total.abs().mean()),
            "fraction_positive":
                float((group.total > 0).mean()),
            "mean_absolute_Ridge":
                float(group.ridge.abs().mean()),
            "mean_absolute_neural_rank":
                float(group.neural_rank.abs().mean()),
            "mean_absolute_SSL":
                float(group.ssl_context.abs().mean()),
            # Measures gross token contributions before cancellation.
            "mean_sum_absolute_token_contributions":
                float(group.sum_absolute_token_contributions.mean()),
        })

        records.append(row)

    for _, group in long_table.groupby(feature_column, sort=False):
        summarize_group(group, "all616", "ALL")

    for (_, cohort), group in long_table.groupby(
        [feature_column, "cohort"], sort=False
    ):
        summarize_group(group, "within_cohort", cohort)

    return pd.DataFrame(records).sort_values(
        ["scope", "cohort", "mean_absolute_contribution"],
        ascending=[True, True, False],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    torch.set_num_threads(8)
    torch.use_deterministic_algorithms(True)

    verify_manifest(args.model, "all67_checkpoint_manifest.tsv")
    verify_manifest(args.data, "export_manifest.tsv")
    verify_manifest(args.context, "context_files_manifest.tsv")

    protocol = json.loads(
        (args.model / "all67_training_protocol.json").read_text(
            encoding="utf-8"
        )
    )

    fusion_path = args.model / "source_fusion_training.py"
    if md5(fusion_path) != protocol["fusion_module_md5"]:
        raise RuntimeError("Fusion implementation changed")

    if md5(args.data / "export_manifest.tsv") != (
        protocol["data_manifest_md5"]
    ):
        raise RuntimeError("Input version changed")

    provenance = json.loads(
        (args.context / "context_provenance.json").read_text(
            encoding="utf-8"
        )
    )
    if provenance != protocol["context_provenance"]:
        raise RuntimeError("Context provenance mismatch")

    spec = importlib.util.spec_from_file_location(
        "fixed_fusion", fusion_path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    hierarchy = pd.read_csv(
        args.data / "token_hierarchy.tsv", sep="\t"
    )
    metadata = pd.read_csv(
        args.data / "ssl_metadata.tsv", sep="\t"
    )
    reference = pd.read_csv(args.scores, sep="\t")

    if len(hierarchy) != 37 or len(metadata) != 616:
        raise RuntimeError("Unexpected dimensions")

    if (
        metadata.patient_uid.duplicated().any()
        or reference.patient_uid.duplicated().any()
        or set(metadata.patient_uid) != set(reference.patient_uid)
    ):
        raise RuntimeError("Patient manifest mismatch")

    reference = reference.set_index("patient_uid").loc[
        metadata.patient_uid
    ]

    token_names = hierarchy.token_name.tolist()

    rank = np.fromfile(
        args.data / "ssl_token_rank.float32.bin",
        dtype="<f4",
    ).reshape(616, 37).astype(np.float64)

    with np.load(
        args.context / "ssl_pretrained_context.npz",
        allow_pickle=False,
    ) as saved:
        if saved["token_names"].tolist() != token_names:
            raise RuntimeError("Token order mismatch")
        if not np.array_equal(
            saved["patient_uid"].astype(str),
            metadata.patient_uid.to_numpy(str),
        ):
            raise RuntimeError("Context patient order mismatch")
        context = saved["token_context"].astype(np.float64)

    components = np.zeros((3, 616, 37), dtype=np.float64)
    seed_predictions = []
    seed_baselines = []
    diagnostic_rows = []
    seed_importance_rows = []

    with threadpool_limits(limits=8):
        for seed in [24, 42, 64]:
            checkpoint_path = (
                args.model / f"P009_BA37_Early_all67_seed{seed}.pt"
            )
            checkpoint = torch.load(
                checkpoint_path,
                map_location="cpu",
                weights_only=False,
            )

            if checkpoint["token_names"] != token_names:
                raise RuntimeError("Checkpoint token order mismatch")

            model = module.Fusion(37, "Early").double()
            model.load_state_dict(
                checkpoint["state_dict"], strict=True
            )
            model.eval()

            for parameter in model.parameters():
                parameter.requires_grad_(False)

            state = checkpoint["scaler_and_Ridge"]
            cap = float(checkpoint["selected_cap"])

            rank_z = (
                (rank - state["rank_mean"]) /
                state["rank_scale"]
            )

            context_z = (
                (
                    context.reshape(616, -1) -
                    state["context_mean"]
                ) / state["context_scale"]
            ).reshape(616, 37, 24)

            coef = torch.tensor(
                np.asarray(state["ridge_coef"]).reshape(-1),
                dtype=torch.float64,
            )
            bias = float(
                np.asarray(state["ridge_intercept"]).reshape(-1)[0]
            )

            with torch.no_grad():
                baseline = float(torch.sigmoid(
                    torch.tensor(bias, dtype=torch.float64) +
                    cap * torch.tanh(model(
                        torch.zeros((1, 37), dtype=torch.float64),
                        torch.zeros((1, 37, 24), dtype=torch.float64),
                    ))[0]
                ))

            seed_components = np.zeros((3, 616, 37))
            predictions = np.zeros(616)

            for start in range(0, 616, 32):
                end = min(start + 32, 616)

                r = torch.tensor(
                    rank_z[start:end], dtype=torch.float64
                )
                c = torch.tensor(
                    context_z[start:end], dtype=torch.float64
                )

                with torch.no_grad():
                    p = torch.sigmoid(
                        r @ coef + bias +
                        cap * torch.tanh(model(r, c))
                    ).numpy()

                for steps in [64, 128, 256]:
                    parts = integrate(
                        model, r, c, coef, bias, cap, steps
                    )

                    total = sum(parts)
                    error = np.abs(
                        total.sum(axis=1) - (p - baseline)
                    )

                    if error.max() <= 1e-5:
                        break

                if error.max() > 1e-4:
                    raise RuntimeError(
                        "IG completeness failed: "
                        f"seed={seed}, batch={start}, "
                        f"max_error={error.max()}"
                    )

                for component_index, part in enumerate(parts):
                    seed_components[
                        component_index, start:end, :
                    ] = part

                predictions[start:end] = p

                diagnostic_rows.append({
                    "seed": seed,
                    "batch_start": start,
                    "integration_nodes": steps,
                    "maximum_completeness_error":
                        float(error.max()),
                })

            # Check inference against previously saved seed scores.
            reference_seed = reference[
                f"score_seed{seed}"
            ].to_numpy()

            if np.max(np.abs(predictions - reference_seed)) > 1e-6:
                raise RuntimeError("Fixed-score reproduction failed")

            components += seed_components / 3
            seed_predictions.append(predictions)
            seed_baselines.append(baseline)

            seed_total = seed_components.sum(axis=0)
            for token_index, token in enumerate(token_names):
                seed_importance_rows.append({
                    "seed": seed,
                    "token_name": token,
                    "mean_absolute_contribution": float(
                        np.abs(seed_total[:, token_index]).mean()
                    ),
                })

            print(f"IG_SEED_COMPLETED: {seed}", flush=True)

    prediction = np.mean(seed_predictions, axis=0)
    baseline = float(np.mean(seed_baselines))
    total_contribution = components.sum(axis=0)

    completeness_error = np.abs(
        total_contribution.sum(axis=1) -
        (prediction - baseline)
    )

    # Construct patient-token long table.
    token_table = metadata[
        ["patient_uid", "cohort"]
    ].loc[
        metadata.index.repeat(37)
    ].reset_index(drop=True)

    for column in [
        "token_name", "concept_id", "concept_name",
        "module_id", "module_name",
    ]:
        token_table[column] = np.tile(
            hierarchy[column].to_numpy(), 616
        )

    token_table["ridge"] = components[0].reshape(-1)
    token_table["neural_rank"] = components[1].reshape(-1)
    token_table["ssl_context"] = components[2].reshape(-1)
    token_table["total"] = total_contribution.reshape(-1)
    token_table["sum_absolute_token_contributions"] = (
        token_table.total.abs()
    )

    values = [
        "ridge", "neural_rank", "ssl_context", "total",
        "sum_absolute_token_contributions",
    ]

    concept_keys = [
        "patient_uid", "cohort",
        "concept_id", "concept_name",
        "module_id", "module_name",
    ]
    module_keys = [
        "patient_uid", "cohort", "module_id", "module_name"
    ]

    concept_table = token_table.groupby(
        concept_keys, as_index=False, sort=False
    )[values].sum()

    module_table = token_table.groupby(
        module_keys, as_index=False, sort=False
    )[values].sum()

    # Verify hierarchical aggregation preserves signed contributions.
    target_sums = total_contribution.sum(axis=1)

    for level_table in [concept_table, module_table]:
        sums = level_table.groupby("patient_uid").total.sum()
        if not np.allclose(
            sums.reindex(metadata.patient_uid).to_numpy(),
            target_sums,
            atol=1e-10,
            rtol=1e-10,
        ):
            raise RuntimeError("Hierarchical aggregation failed")

    for level, table in [
        ("token", token_table),
        ("concept", concept_table),
        ("module", module_table),
    ]:
        table.to_csv(
            args.output / f"{level}_patient_contributions.tsv",
            sep="\t", index=False,
        )
        summarize(table, level).to_csv(
            args.output / f"{level}_importance_summary.tsv",
            sep="\t", index=False,
        )

    patient_check = metadata[
        ["patient_uid", "cohort"]
    ].copy()
    patient_check["P009_BA37_model_score"] = prediction
    patient_check["reference_score"] = baseline
    patient_check["sum_contributions"] = target_sums
    patient_check["completeness_error"] = completeness_error

    patient_check.to_csv(
        args.output / "patient_attribution_completeness.tsv",
        sep="\t", index=False,
    )

    pd.DataFrame(diagnostic_rows).to_csv(
        args.output / "integration_diagnostics.tsv",
        sep="\t", index=False,
    )
    pd.DataFrame(seed_importance_rows).to_csv(
        args.output / "token_importance_by_seed.tsv",
        sep="\t", index=False,
    )

    report = {
        "model": "fixed P009_BA37-only Early all67 ensemble",
        "patients": 616,
        "token_count": 37,
        "concept_count": int(hierarchy.concept_id.nunique()),
        "module_count": int(hierarchy.module_id.nunique()),
        "target": "mean predicted probability over three seeds",
        "baseline": "zero standardized rank and context",
        "baseline_score": baseline,
        "maximum_completeness_error":
            float(completeness_error.max()),
        "primary_importance": "mean absolute patient-level net contribution",
        "group_contributions": "sum signed member token contributions",
        "interpretation": (
            "Model attribution, not causal biology. "
            "Rank and context are jointly varied along a synthetic path."
        ),
        "model_protocol_md5": md5(
            args.model / "all67_training_protocol.json"
        ),
        "source_scores_md5": md5(args.scores),
        "script_md5": md5(__file__),
    }

    (args.output / "attribution_report.json").write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(report, indent=2), flush=True)
    print("HIERARCHICAL_ATTRIBUTION_COMPLETED", flush=True)


if __name__ == "__main__":
    main()

