
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch


def effective_rank(values):
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    normalized = values / np.maximum(norms, 1e-12)
    centered = normalized - normalized.mean(axis=0, keepdims=True)

    eigenvalues = np.linalg.svd(
        centered, compute_uv=False
    ) ** 2

    if eigenvalues.sum() <= 1e-20:
        return 0.0

    probabilities = eigenvalues / eigenvalues.sum()
    probabilities = probabilities[probabilities > 0]

    return float(np.exp(
        -(probabilities * np.log(probabilities)).sum()
    ))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--initialization", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    torch.set_num_threads(8)
    torch.use_deterministic_algorithms(True)

    spec = importlib.util.spec_from_file_location(
        "p009_ssl_model", args.model
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    summary = json.loads(
        (args.training / "training_summary.json").read_text(
            encoding="utf-8"
        )
    )

    expected_checkpoint_md5 = summary["best_checkpoint_md5"]

    if module.md5(args.checkpoint) != expected_checkpoint_md5:
        raise RuntimeError("Unexpected SSL checkpoint version")

    if (
        summary["best_checkpoint_md5"] != expected_checkpoint_md5
        or not summary["completed"]
        or summary["response_labels_used"]
    ):
        raise RuntimeError("Training summary validation failed")

    # This loader validates hashes and indices.
    # For ICI cohorts it reads patient_uid only, not response values.
    arrays, metadata, membership, tc, cm, model_spec = (
        module.load_inputs(args.data)
    )

    if (
        model_spec["n_tokens"] != 37
        or model_spec["n_genes"] != 1753
        or model_spec["n_concepts"] != 22
        or model_spec["n_modules"] != 8
    ):
        raise RuntimeError("Unexpected P009_BA37 dimensions")

    best = torch.load(
        args.checkpoint,
        map_location="cpu",
        weights_only=False,
    )

    initial = torch.load(
        args.initialization,
        map_location="cpu",
        weights_only=False,
    )

    if best["spec"] != model_spec or initial["spec"] != model_spec:
        raise RuntimeError("Checkpoint spec mismatch")

    if initial["optimizer_steps"] != 0:
        raise RuntimeError("Initialization has optimizer updates")

    manifest_md5 = module.md5(args.data / "export_manifest.tsv")

    if best["config"]["data_manifest_md5"] != manifest_md5:
        raise RuntimeError("Training data version mismatch")

    if initial["data_manifest_md5"] != manifest_md5:
        raise RuntimeError("Initialization data version mismatch")

    if best["config"]["model_script_md5"] != module.md5(args.model):
        raise RuntimeError("Model script version mismatch")

    if best["config"]["initialization_md5"] != module.md5(
        args.initialization
    ):
        raise RuntimeError("Initialization version mismatch")

    args.output.mkdir(parents=True, exist_ok=True)

    if (args.output / "context_provenance.json").exists():
        raise RuntimeError("Use a new context output directory")

    provenance = {
        "experiment": "P009_BA37 SSL",
        "pretrained_checkpoint": str(args.checkpoint),
        "pretrained_checkpoint_md5": module.md5(args.checkpoint),
        "pretrained_epoch": int(best["epoch"]),
        "initialization": str(args.initialization),
        "initialization_md5": module.md5(args.initialization),
        "model_script": str(args.model),
        "model_script_md5": module.md5(args.model),
        "data_directory": str(args.data),
        "data_manifest_md5": manifest_md5,
        "context": "member-local token vectors before module attention",
        "n_tokens": 37,
        "n_genes": 1753,
        "context_dim": 24,
        "ICI_response_values_read": False,
        "token_selection_applied": False,
        "panel_fixed_before_this_experiment": "P009_BA37",
    }

    # Write provenance before extraction, recording the planned source.
    (args.output / "context_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    diagnostic_rows = []
    output_rows = []

    for encoder_name, checkpoint in [
        ("pretrained", best),
        ("random", initial),
    ]:
        model = module.SSL37(
            membership, tc, cm, dim=24, heads=4
        )
        model.load_state_dict(
            checkpoint["state_dict"], strict=True
        )
        model.eval()

        for dataset, values in arrays.items():
            pieces = []

            with torch.no_grad():
                for start in range(0, len(values["x"]), 16):
                    end = min(start + 16, len(values["x"]))

                    x = torch.from_numpy(
                        values["x"][start:end].copy()
                    )
                    observed = torch.from_numpy(
                        values["observed"][start:end].copy()
                    )

                    tokens, _, _ = model.encode(x, observed)
                    pieces.append(tokens.numpy())

            context = np.concatenate(pieces).astype(np.float32)

            if context.shape != (len(values["ids"]), 37, 24):
                raise RuntimeError("Incorrect context shape")

            if not np.isfinite(context).all():
                raise RuntimeError("Nonfinite context values")

            # Verify clean inference is independent of patient batching.
            with torch.no_grad():
                one_x = torch.from_numpy(
                    values["x"][:1].copy()
                )
                one_mask = torch.from_numpy(
                    values["observed"][:1].copy()
                )
                one = model.encode(one_x, one_mask)[0].numpy()

            if not np.allclose(
                one, context[:1], atol=1e-5, rtol=1e-5
            ):
                raise RuntimeError("Batch independence failed")

            path = args.output / (
                f"{dataset}_{encoder_name}_context.npz"
            )

            np.savez_compressed(
                path,
                token_context=context,
                patient_uid=np.asarray(
                    values["ids"], dtype=str
                ),
                token_names=np.asarray(
                    model_spec["token_names"], dtype=str
                ),
            )

            # Read the saved artifact back and verify exact values/order.
            with np.load(path, allow_pickle=False) as saved:
                if not np.array_equal(
                    saved["token_context"], context
                ):
                    raise RuntimeError("Saved context value mismatch")

                if not np.array_equal(
                    saved["patient_uid"],
                    np.asarray(values["ids"], dtype=str),
                ):
                    raise RuntimeError("Saved patient order mismatch")

                if saved["token_names"].tolist() != (
                    model_spec["token_names"]
                ):
                    raise RuntimeError("Saved token order mismatch")

            patient_sd = context.std(axis=0)
            token_ranks = [
                effective_rank(context[:, token, :])
                for token in range(37)
            ]

            diagnostic_rows.append({
                "encoder": encoder_name,
                "dataset": dataset,
                "n_samples": len(context),
                "n_tokens": 37,
                "context_dim": 24,
                "mean_token_patient_SD":
                    float(patient_sd.mean()),
                "minimum_token_patient_SD":
                    float(patient_sd.mean(axis=1).min()),
                "median_token_effective_rank":
                    float(np.median(token_ranks)),
                "batch_independence_pass": True,
                "saved_values_and_order_pass": True,
            })

            output_rows.append({
                "filename": path.name,
                "size_bytes": path.stat().st_size,
                "md5": module.md5(path),
            })

            print(
                f"EXTRACTED {encoder_name} {dataset} "
                f"shape={tuple(context.shape)}",
                flush=True,
            )

    diagnostics = pd.DataFrame(diagnostic_rows)

    diagnostics.to_csv(
        args.output / "context_diagnostics.tsv",
        sep="\t", index=False,
    )

    pd.DataFrame(output_rows).to_csv(
        args.output / "context_files_manifest.tsv",
        sep="\t", index=False,
    )

    (args.output / "BA37_model_spec.json").write_text(
        json.dumps(model_spec, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("\nCONTEXT_DIAGNOSTICS:", flush=True)
    print(diagnostics.to_string(index=False))
    print("BA37_CONTEXT_EXTRACTION_COMPLETED", flush=True)


if __name__ == "__main__":
    main()

