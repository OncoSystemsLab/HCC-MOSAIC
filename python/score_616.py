
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
    manifest = pd.read_csv(directory / filename, sep="\t")
    for row in manifest.itertuples(index=False):
        path = directory / row.filename
        if not path.exists() or md5(path) != row.md5:
            raise RuntimeError(f"Fingerprint mismatch: {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
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
        raise RuntimeError("Fusion implementation version mismatch")

    if md5(args.data / "export_manifest.tsv") != (
        protocol["data_manifest_md5"]
    ):
        raise RuntimeError("Input version mismatch")

    provenance = json.loads(
        (args.context / "context_provenance.json").read_text(
            encoding="utf-8"
        )
    )

    if provenance != protocol["context_provenance"]:
        raise RuntimeError("SSL context provenance mismatch")

    if provenance["data_manifest_md5"] != (
        protocol["data_manifest_md5"]
    ):
        raise RuntimeError("Context and expression versions differ")

    spec = importlib.util.spec_from_file_location(
        "fixed_fusion", fusion_path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    hierarchy = pd.read_csv(
        args.data / "token_hierarchy.tsv", sep="\t"
    )
    token_names = hierarchy.token_name.tolist()

    if len(token_names) != 37:
        raise RuntimeError("Expected 37 tokens")

    shapes = pd.read_csv(args.data / "shapes.tsv", sep="\t")

    checkpoints = {}
    for seed in [24, 42, 64]:
        checkpoint = torch.load(
            args.model / f"P009_BA37_Early_all67_seed{seed}.pt",
            map_location="cpu",
            weights_only=False,
        )

        if (
            checkpoint["architecture"] != "Early"
            or checkpoint["token_names"] != token_names
            or checkpoint["seed"] != seed
            or not np.isclose(checkpoint["selected_cap"], protocol["cap"])
        ):
            raise RuntimeError("Fixed checkpoint configuration mismatch")

        checkpoints[seed] = checkpoint

    outputs = {}
    source_records = []

    with threadpool_limits(limits=8):
        for dataset in [
            "ssl", "GSE181946", "GSE202069", "PRJEB34724"
        ]:
            row = shapes.loc[shapes.dataset == dataset].iloc[0]
            n = int(row.n_samples)

            metadata_path = args.data / f"{dataset}_metadata.tsv"
            metadata = pd.read_csv(metadata_path, sep="\t")

            if len(metadata) != n or not metadata.patient_uid.is_unique:
                raise RuntimeError("Invalid metadata")

            rank_path = args.data / f"{dataset}_token_rank.float32.bin"
            rank = np.fromfile(
                rank_path, dtype="<f4"
            ).reshape(n, 37).astype(np.float64)

            context_path = (
                args.context / f"{dataset}_pretrained_context.npz"
            )

            with np.load(context_path, allow_pickle=False) as saved:
                if saved["token_names"].tolist() != token_names:
                    raise RuntimeError("Token order mismatch")
                if not np.array_equal(
                    saved["patient_uid"].astype(str),
                    metadata.patient_uid.to_numpy(str),
                ):
                    raise RuntimeError("Patient order mismatch")
                context = saved["token_context"].astype(np.float64)

            if (
                context.shape != (n, 37, 24)
                or not np.isfinite(context).all()
                or not np.isfinite(rank).all()
            ):
                raise RuntimeError("Invalid inference values")

            scored = metadata.copy()
            seed_predictions = []

            for seed, checkpoint in checkpoints.items():
                state = checkpoint["scaler_and_Ridge"]

                rank_z = (
                    (rank - state["rank_mean"]) /
                    state["rank_scale"]
                ).astype(np.float32)

                context_z = (
                    (
                        context.reshape(n, -1) -
                        state["context_mean"]
                    ) / state["context_scale"]
                ).reshape(n, 37, 24).astype(np.float32)

                base_logit = (
                    rank_z @ state["ridge_coef"].T +
                    state["ridge_intercept"]
                ).reshape(-1).astype(np.float32)

                model = module.Fusion(37, "Early")
                model.load_state_dict(
                    checkpoint["state_dict"], strict=True
                )
                model.eval()

                with torch.no_grad():
                    raw = model(
                        torch.from_numpy(rank_z),
                        torch.from_numpy(context_z),
                    )
                    correction = float(checkpoint["selected_cap"]) * torch.tanh(raw)
                    probability = torch.sigmoid(
                        torch.from_numpy(base_logit) + correction
                    ).numpy()

                if not np.isfinite(probability).all():
                    raise RuntimeError("Nonfinite predictions")

                scored[f"score_seed{seed}"] = probability
                seed_predictions.append(probability)

            scored["P009_BA37_model_score"] = np.mean(
                seed_predictions, axis=0
            )
            outputs[dataset] = scored

            source_records.append({
                "dataset": dataset,
                "rank_file_md5": md5(rank_path),
                "context_file_md5": md5(context_path),
                "metadata_file_md5": md5(metadata_path),
            })

    # Inference must reproduce the previously saved all-67 scores.
    ici = pd.concat(
        [outputs[c] for c in [
            "GSE181946", "GSE202069", "PRJEB34724"
        ]],
        ignore_index=True,
    )

    reference = pd.read_csv(
        args.model / "P009_BA37_all67_fixed_model_scores.tsv",
        sep="\t",
    )

    if (
        len(reference) != 67
        or reference.patient_uid.duplicated().any()
    ):
        raise RuntimeError("Invalid reference patient manifest")

    checked = ici.merge(
        reference[
            ["patient_uid", "P009_BA37_fixed_model_score"]
        ],
        on="patient_uid",
        how="left",
        validate="one_to_one",
    )

    if checked.P009_BA37_fixed_model_score.isna().any():
        raise RuntimeError("Reference patient mismatch")

    error = float(np.max(np.abs(
        checked.P009_BA37_model_score -
        checked.P009_BA37_fixed_model_score
    )))

    if error > 1e-6:
        raise RuntimeError(
            f"Fixed-score reproduction failed: {error}"
        )

    scored616 = outputs["ssl"]

    if len(scored616) != 616:
        raise RuntimeError("Expected 616 HCC patients")

    expected = {"TCGA": 361, "ICGC": 96, "CPTAC": 159}
    if scored616.cohort.value_counts().to_dict() != expected:
        raise RuntimeError("Unexpected HCC cohort counts")

    scored616.to_csv(
        args.output / "P009_BA37_fixed_model_scores_616.tsv",
        sep="\t", index=False,
    )

    pd.DataFrame(source_records).to_csv(
        args.output / "inference_source_fingerprints.tsv",
        sep="\t", index=False,
    )

    report = {
        "model": "P009_BA37-only Early all67 fixed ensemble",
        "n_HCC_scored": 616,
        "reference67_max_prediction_error": error,
        "supervised_model_refitted": False,
        "scalers_refitted": False,
        "SSL_encoder_updated": False,
        "HCC_participated_in_SSL_pretraining": True,
        "purpose": "mutation-score association exploration",
        "script_md5": md5(__file__),
        "fixed_model_protocol_md5": md5(
            args.model / "all67_training_protocol.json"
        ),
    }

    (args.output / "inference_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )

    print("REFERENCE67_MAX_ERROR:", error, flush=True)
    print("HCC_COHORT_COUNTS:", expected, flush=True)
    print("FIXED_P009_BA37_616_SCORING_COMPLETED", flush=True)


if __name__ == "__main__":
    main()

