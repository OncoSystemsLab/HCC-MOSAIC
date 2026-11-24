
import argparse
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--loco", type=Path, required=True)
    parser.add_argument("--fusion", type=Path, required=True)
    parser.add_argument("--loader", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    torch.set_num_threads(8)
    torch.use_deterministic_algorithms(True)

    fusion = load_module("fixed_fusion", args.fusion)
    loader = load_module("fixed_loader", args.loader)

    data = loader.load_data(args.data, args.context)

    # Check that the fixed cap matches all previous outer folds.
    previous_metrics = pd.read_csv(
        args.loco / "outer_cohort_metrics.tsv", sep="\t"
    )

    selected = previous_metrics.loc[
        (previous_metrics.panel_id == "P009_BA37") &
        (previous_metrics["mode"] == "Early")
    ]

    if len(selected) != 3:
        raise RuntimeError("Expected three completed LOCO folds")
    cap_counts = selected.selected_cap.value_counts().to_dict()
    selected_cap = float(sorted(cap_counts, key=lambda x: (-cap_counts[x], x))[0])
    if selected_cap not in [0.25, 0.75] or cap_counts[selected_cap] < 2:
        raise RuntimeError("No valid majority cap")

    if (
        len(data["hierarchy"]) != 37
        or len(data["meta"]) != 67
    ):
        raise RuntimeError("Unexpected panel or patient count")

    indices = list(range(37))
    all_rows = np.ones(len(data["y"]), dtype=bool)

    protocol = {
        "model_id": "P009_BA37only_Early_all67_refit_v1",
        "training_n": 67,
        "panel": "P009_BA37",
        "n_tokens": 37,
        "SSL_pretraining_n": 616,
        "SSL_encoder_updated": False,
        "mode": "Early",
        "cap": selected_cap,
        "cap_rule": "majority cap among prior LOCO folds; fixed before refit",
        "epochs": 60,
        "seeds": [24, 42, 64],
        "prediction_rule": "mean of three seed probabilities",
        "evaluation": "training-set resubstitution, not LOCO",
        "fusion_module_md5": loader.md5(args.fusion),
        "loader_module_md5": loader.md5(args.loader),
        "data_manifest_md5": loader.md5(
            args.data / "export_manifest.tsv"
        ),
        "context_provenance": data["provenance"],
        "context_files_md5": {
            cohort: loader.md5(
                args.context /
                f"{cohort}_pretrained_context.npz"
            )
            for cohort in [
                "GSE181946", "GSE202069", "PRJEB34724"
            ]
        },
    }

    protocol_path = args.output / "all67_training_protocol.json"

    if protocol_path.exists():
        raise RuntimeError("Use a new output directory")

    protocol_path.write_text(
        json.dumps(protocol, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    # Snapshot executable dependencies.
    shutil.copy2(
        args.fusion,
        args.output / "source_fusion_training.py",
    )
    shutil.copy2(
        args.loader,
        args.output / "source_data_loader.py",
    )

    with threadpool_limits(limits=8):
        prepared = fusion.prepare(
            data, indices, all_rows, all_rows
        )

        predictions = []
        diagnostic_rows = []
        seed_rows = []

        for seed in [24, 42, 64]:
            p, state, context_change, correction = (
                fusion.fit_predict(
                    prepared,
                    mode="Early",
                    cap=selected_cap,
                    seed=seed,
                )
            )

            predictions.append(p)

            checkpoint_path = args.output / (
                f"P009_BA37_Early_all67_seed{seed}.pt"
            )

            torch.save({
                "state_dict": state,
                "architecture": "Early",
                "panel_id": "P009_BA37",
                "token_indices": indices,
                "token_names":
                    data["hierarchy"].token_name.tolist(),
                "selected_cap": selected_cap,
                "seed": seed,
                "training_patient_uids":
                    data["meta"].patient_uid.tolist(),
                "scaler_and_Ridge":
                    prepared["scaler_state"],
                "protocol": protocol,
            }, checkpoint_path)

            for i, patient in data["meta"].iterrows():
                seed_rows.append({
                    "patient_uid": patient.patient_uid,
                    "cohort": patient.cohort,
                    "response": int(patient.response),
                    "seed": seed,
                    "score": float(p[i]),
                })

            diagnostic_rows.append({
                "seed": seed,
                "mean_abs_SSL_ablation_logit_change":
                    context_change,
                "mean_abs_logit_correction": correction,
            })

            print(f"ALL67_REFIT seed={seed}: COMPLETED", flush=True)

    final_prediction = np.mean(predictions, axis=0)

    scored = data["meta"][
        ["patient_uid", "cohort", "response"]
    ].copy()

    scored["P009_BA37_fixed_model_score"] = final_prediction
    scored["prediction_type"] = "all67_training_set_prediction"

    # Compare with previous LOCO scores for the same patients.
    old = pd.read_csv(
        args.loco / "patient_predictions.tsv", sep="\t"
    )

    old = old.loc[
        (old.panel_id == "P009_BA37") &
        (old["mode"] == "Early") &
        (old.seed.astype(str) == "ensemble"),
        ["patient_uid", "prediction"],
    ].rename(columns={"prediction": "P009_BA37_LOCO_score"})

    if len(old) != 67 or old.patient_uid.duplicated().any():
        raise RuntimeError("Invalid previous LOCO predictions")

    scored = scored.merge(
        old, on="patient_uid", how="left", validate="one_to_one"
    )

    if scored.P009_BA37_LOCO_score.isna().any():
        raise RuntimeError("Missing LOCO scores")

    scored["fixed_minus_LOCO"] = (
        scored.P009_BA37_fixed_model_score -
        scored.P009_BA37_LOCO_score
    )

    scored.to_csv(
        args.output / "P009_BA37_all67_fixed_model_scores.tsv",
        sep="\t", index=False,
    )

    pd.DataFrame(seed_rows).to_csv(
        args.output / "P009_BA37_all67_seed_scores.tsv",
        sep="\t", index=False,
    )

    pd.DataFrame(diagnostic_rows).to_csv(
        args.output / "all67_SSL_input_diagnostics.tsv",
        sep="\t", index=False,
    )

    source_files = [
        args.output / f"P009_BA37_Early_all67_seed{s}.pt"
        for s in [24, 42, 64]
    ]

    pd.DataFrame([
        {
            "filename": p.name,
            "size_bytes": p.stat().st_size,
            "md5": loader.md5(p),
        }
        for p in source_files
    ]).to_csv(
        args.output / "all67_checkpoint_manifest.tsv",
        sep="\t", index=False,
    )

    print("P009_BA37_ALL67_REFIT_COMPLETED", flush=True)


if __name__ == "__main__":
    main()

