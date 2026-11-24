
import argparse
import copy
import importlib.util
import json
import random
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    brier_score_loss,
)
from threadpoolctl import threadpool_limits


COHORTS = ["GSE181946", "GSE202069", "PRJEB34724"]
MODES = ["Early"]


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def weights(y, cohorts):
    w = np.zeros(len(y), dtype=np.float32)

    for cohort in np.unique(cohorts):
        for label in [0, 1]:
            hit = (cohorts == cohort) & (y == label)

            if hit.sum() == 0:
                raise ValueError("A training cohort lacks a class")

            w[hit] = 1.0 / hit.sum()

    return w * len(w) / w.sum()


class Fusion(nn.Module):
    def __init__(self, n_tokens, mode):
        super().__init__()
        self.mode = mode

        self.context_projection = nn.Sequential(
            nn.Linear(24, 8),
            nn.GELU(),
        )

        if mode == "Early":
            self.body = nn.Sequential(
                nn.Linear(n_tokens * 9, 16),
                nn.GELU(),
                nn.Dropout(0.15),
            )

        elif mode == "Late":
            self.rank_branch = nn.Sequential(
                nn.Linear(n_tokens, 8),
                nn.GELU(),
            )
            self.context_branch = nn.Sequential(
                nn.Linear(n_tokens * 8, 8),
                nn.GELU(),
            )
            self.body = nn.Sequential(
                nn.Linear(16, 16),
                nn.GELU(),
                nn.Dropout(0.15),
            )

        elif mode == "Gated":
            self.rank_projection = nn.Linear(1, 8)
            self.gate = nn.Linear(16, 1)
            self.body = nn.Sequential(
                nn.Linear(n_tokens * 8, 16),
                nn.GELU(),
                nn.Dropout(0.15),
            )

        else:
            raise ValueError(mode)

        self.output = nn.Linear(16, 1)

        # The initial model exactly equals the Ridge backbone.
        nn.init.zeros_(self.output.weight)
        nn.init.zeros_(self.output.bias)

    def forward(self, rank, context):
        projected = self.context_projection(context)

        if self.mode == "Early":
            joined = torch.cat(
                [rank.unsqueeze(-1), projected], dim=-1
            )
            hidden = self.body(joined.flatten(1))

        elif self.mode == "Late":
            rank_hidden = self.rank_branch(rank)
            context_hidden = self.context_branch(
                projected.flatten(1)
            )
            hidden = self.body(torch.cat(
                [rank_hidden, context_hidden], dim=1
            ))

        else:
            rank_hidden = self.rank_projection(
                rank.unsqueeze(-1)
            )
            gate = torch.sigmoid(self.gate(torch.cat(
                [rank_hidden, projected], dim=-1
            )))
            fused = (
                (1 - gate) * rank_hidden +
                gate * projected
            )
            hidden = self.body(fused.flatten(1))

        return self.output(hidden).squeeze(-1)


def metric(y, p):
    return {
        "AUC": float(roc_auc_score(y, p)),
        "average_precision": float(
            average_precision_score(y, p)
        ),
        "Brier_score": float(brier_score_loss(y, p)),
    }


def prepare(data, indices, train, test):
    r = data["rank"][:, indices].astype(np.float64)
    c = data["ssl"][:, indices].astype(np.float64)

    rank_scaler = StandardScaler().fit(r[train])
    context_scaler = StandardScaler().fit(
        c[train].reshape(train.sum(), -1)
    )

    rt = rank_scaler.transform(r[train]).astype(np.float32)
    rv = rank_scaler.transform(r[test]).astype(np.float32)

    ct = context_scaler.transform(
        c[train].reshape(train.sum(), -1)
    ).reshape(train.sum(), len(indices), 24).astype(np.float32)

    cv = context_scaler.transform(
        c[test].reshape(test.sum(), -1)
    ).reshape(test.sum(), len(indices), 24).astype(np.float32)

    w = weights(data["y"][train], data["cohort"][train])

    ridge = LogisticRegression(
        C=0.1,
        solver="lbfgs",
        max_iter=10000,
        tol=1e-8,
        random_state=42,
    )

    ridge.fit(rt, data["y"][train], sample_weight=w)

    if ridge.n_iter_.max() >= 10000:
        raise RuntimeError("Ridge did not converge")

    result = {
        "rank_train": torch.from_numpy(rt),
        "rank_test": torch.from_numpy(rv),
        "context_train": torch.from_numpy(ct),
        "context_test": torch.from_numpy(cv),
        "base_train": torch.tensor(
            ridge.decision_function(rt), dtype=torch.float32
        ),
        "base_test": torch.tensor(
            ridge.decision_function(rv), dtype=torch.float32
        ),
        "y": torch.tensor(
            data["y"][train], dtype=torch.float32
        ),
        "weights": torch.from_numpy(w),
        "scaler_state": {
            "rank_mean": rank_scaler.mean_,
            "rank_scale": rank_scaler.scale_,
            "context_mean": context_scaler.mean_,
            "context_scale": context_scaler.scale_,
            "ridge_coef": ridge.coef_,
            "ridge_intercept": ridge.intercept_,
        },
    }

    return result


def fit_predict(prepared, mode, cap, seed):
    seed_all(seed)

    model = Fusion(
        prepared["rank_train"].shape[1],
        mode,
    )

    model.eval()
    with torch.no_grad():
        initial = model(
            prepared["rank_test"],
            prepared["context_test"],
        )
        if not torch.equal(initial, torch.zeros_like(initial)):
            raise RuntimeError("Initial correction is not zero")

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=0.001,
        weight_decay=0.05,
    )

    for epoch in range(60):
        model.train()
        optimizer.zero_grad(set_to_none=True)

        raw = model(
            prepared["rank_train"],
            prepared["context_train"],
        )

        correction = cap * torch.tanh(raw)
        logits = prepared["base_train"] + correction

        loss = (
            F.binary_cross_entropy_with_logits(
                logits, prepared["y"], reduction="none"
            ) * prepared["weights"]
        ).mean()

        loss = loss + 0.01 * correction.square().mean()

        if not torch.isfinite(loss):
            raise RuntimeError("Nonfinite training loss")

        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

    model.eval()

    with torch.no_grad():
        raw = model(
            prepared["rank_test"],
            prepared["context_test"],
        )
        correction = cap * torch.tanh(raw)
        prediction = torch.sigmoid(
            prepared["base_test"] + correction
        ).numpy()

        # Verify actual dependence on the SSL input.
        zero_context_raw = model(
            prepared["rank_test"],
            torch.zeros_like(prepared["context_test"]),
        )

        context_change = (
            cap * torch.tanh(raw) -
            cap * torch.tanh(zero_context_raw)
        ).abs().mean().item()

    if not np.isfinite(prediction).all():
        raise RuntimeError("Nonfinite predictions")

    return (
        prediction,
        copy.deepcopy(model.state_dict()),
        float(context_change),
        float(correction.abs().mean()),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--context", required=True)
    parser.add_argument("--panels", required=True)
    parser.add_argument("--loader", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    started = time.time()
    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)

    torch.set_num_threads(8)
    torch.use_deterministic_algorithms(True)

    spec = importlib.util.spec_from_file_location(
        "previous_supervised", args.loader
    )
    previous = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(previous)

    # Includes exported-input hashes, SSL checkpoint identity,
    # patient order and token order checks.
    data = previous.load_data(
        Path(args.data), Path(args.context)
    )

    data["cohort"] = data["meta"].cohort.to_numpy()

    hierarchy = data["hierarchy"]
    lookup = {
        name: i
        for i, name in enumerate(hierarchy.token_name)
    }

    roster_path = (
        Path(args.panels) / "candidate_panel_membership.tsv"
    )
    roster = pd.read_csv(roster_path, sep="\t")

    panels = {}
    for name, expected_size in [("P009_BA37", 37)]:
        tokens = roster.loc[
            roster.panel_id == name, "token_name"
        ].tolist()

        if (
            len(tokens) != expected_size
            or len(set(tokens)) != expected_size
        ):
            raise ValueError("Panel roster mismatch")

        panels[name] = sorted(lookup[t] for t in tokens)

    protocol = {
        "panels": ["P009_BA37"],
        "architectures": MODES,
        "inner_seeds": [24, 42],
        "outer_seeds": [24, 42, 64],
        "correction_caps": [0.25, 0.75],
        "epochs": 60,
        "learning_rate": 0.001,
        "weight_decay": 0.05,
        "correction_penalty": 0.01,
        "Ridge_C": 0.1,
        "encoder_frozen": True,
        "SSL_source": "616 HCC pretrained token contexts",
        "inner_selection": "0.7*mean_AUC + 0.3*minimum_AUC",
        "tie_break": "smaller correction cap",
        "final_prediction": "mean probabilities over three seeds",
        "development_label_informed": True,
        "independent_validation": False,
        "script_md5": previous.md5(Path(__file__)),
        "loader_md5": previous.md5(Path(args.loader)),
        "panel_roster_md5": previous.md5(roster_path),
        "context_files_md5": {
            cohort: previous.md5(
                Path(args.context) /
                f"{cohort}_pretrained_context.npz"
            )
            for cohort in COHORTS
        },
        "context_provenance": json.loads(
            (Path(args.context) /
             "context_provenance.json").read_text(
                 encoding="utf-8"
             )
        ),
    }

    protocol_path = root / "training_protocol.json"
    if protocol_path.exists():
        raise RuntimeError("Use a new output directory")

    protocol_path.write_text(
        json.dumps(protocol, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    inner_rows = []
    outer_rows = []
    prediction_rows = []
    seed_metric_rows = []

    cohorts = data["cohort"]
    y = data["y"]

    with threadpool_limits(limits=8):
        for panel, indices in panels.items():
            for held in COHORTS:
                outer_train = cohorts != held
                outer_test = cohorts == held

                remaining = [
                    c for c in COHORTS if c != held
                ]

                outer_prepared = prepare(
                    data, indices, outer_train, outer_test
                )

                inner_prepared = []
                for validation_cohort in remaining:
                    inner_train = (
                        outer_train &
                        (cohorts != validation_cohort)
                    )
                    inner_test = cohorts == validation_cohort

                    inner_prepared.append((
                        validation_cohort,
                        inner_test,
                        prepare(
                            data, indices,
                            inner_train, inner_test
                        ),
                    ))

                # Baseline on the exact same outer training split.
                base_prediction = torch.sigmoid(
                    outer_prepared["base_test"]
                ).numpy()

                baseline_metrics = metric(
                    y[outer_test], base_prediction
                )

                outer_rows.append({
                    "panel_id": panel,
                    "mode": "Ridge",
                    "cohort": held,
                    "selected_cap": 0.0,
                    **baseline_metrics,
                })

                test_meta = data["meta"].loc[
                    outer_test
                ].reset_index(drop=True)

                for k, patient in test_meta.iterrows():
                    prediction_rows.append({
                        "panel_id": panel,
                        "mode": "Ridge",
                        "cohort": held,
                        "patient_uid": patient.patient_uid,
                        "response": int(patient.response),
                        "seed": "ensemble",
                        "prediction": float(base_prediction[k]),
                    })

                for mode in MODES:
                    choices = []

                    for cap in [0.25, 0.75]:
                        fold_aucs = []

                        for val_cohort, mask, prepared in inner_prepared:
                            probabilities = []

                            for seed in [24, 42]:
                                p, _, _, _ = fit_predict(
                                    prepared, mode, cap, seed
                                )
                                probabilities.append(p)

                            auc = roc_auc_score(
                                y[mask],
                                np.mean(probabilities, axis=0),
                            )
                            fold_aucs.append(float(auc))

                            inner_rows.append({
                                "panel_id": panel,
                                "outer_heldout": held,
                                "mode": mode,
                                "cap": cap,
                                "inner_validation": val_cohort,
                                "AUC": float(auc),
                            })

                        score = (
                            0.7 * np.mean(fold_aucs) +
                            0.3 * np.min(fold_aucs)
                        )
                        choices.append((float(score), cap))

                    selected_score, selected_cap = sorted(
                        choices,
                        key=lambda x: (-x[0], x[1]),
                    )[0]

                    probabilities = []
                    changes = []
                    corrections = []

                    for seed in [24, 42, 64]:
                        p, state, change, correction = fit_predict(
                            outer_prepared,
                            mode,
                            selected_cap,
                            seed,
                        )

                        probabilities.append(p)
                        changes.append(change)
                        corrections.append(correction)

                        seed_metric_rows.append({
                            "panel_id": panel,
                            "mode": mode,
                            "cohort": held,
                            "seed": seed,
                            **metric(y[outer_test], p),
                        })

                        for k, patient in test_meta.iterrows():
                            prediction_rows.append({
                                "panel_id": panel,
                                "mode": mode,
                                "cohort": held,
                                "patient_uid": patient.patient_uid,
                                "response": int(patient.response),
                                "seed": str(seed),
                                "prediction": float(p[k]),
                            })

                        torch.save({
                            "state_dict": state,
                            "architecture": mode,
                            "panel_id": panel,
                            "token_indices": indices,
                            "token_names": [
                                hierarchy.iloc[i].token_name
                                for i in indices
                            ],
                            "selected_cap": selected_cap,
                            "seed": seed,
                            "heldout_cohort": held,
                            "training_patient_uids":
                                data["meta"].loc[
                                    outer_train, "patient_uid"
                                ].tolist(),
                            "scaler_and_Ridge":
                                outer_prepared["scaler_state"],
                            "protocol": protocol,
                        }, root / (
                            f"{panel}_{mode}_{held}_seed{seed}.pt"
                        ))

                    prediction = np.mean(
                        probabilities, axis=0
                    )

                    fold_metrics = metric(
                        y[outer_test], prediction
                    )

                    outer_rows.append({
                        "panel_id": panel,
                        "mode": mode,
                        "cohort": held,
                        "selected_cap": selected_cap,
                        "selected_inner_score": selected_score,
                        "mean_abs_SSL_ablation_logit_change":
                            float(np.mean(changes)),
                        "mean_abs_logit_correction":
                            float(np.mean(corrections)),
                        **fold_metrics,
                    })

                    for k, patient in test_meta.iterrows():
                        prediction_rows.append({
                            "panel_id": panel,
                            "mode": mode,
                            "cohort": held,
                            "patient_uid": patient.patient_uid,
                            "response": int(patient.response),
                            "seed": "ensemble",
                            "prediction": float(prediction[k]),
                        })

                    print(
                        f"{panel} {mode} held={held} "
                        f"cap={selected_cap} "
                        f"AUC={fold_metrics['AUC']:.6f}",
                        flush=True,
                    )

                    pd.DataFrame(outer_rows).to_csv(
                        root / "outer_cohort_metrics.tsv",
                        sep="\t", index=False,
                    )
                    pd.DataFrame(inner_rows).to_csv(
                        root / "inner_selection.tsv",
                        sep="\t", index=False,
                    )
                    pd.DataFrame(prediction_rows).to_csv(
                        root / "patient_predictions.tsv",
                        sep="\t", index=False,
                    )
                    pd.DataFrame(seed_metric_rows).to_csv(
                        root / "individual_seed_metrics.tsv",
                        sep="\t", index=False,
                    )

    metrics = pd.DataFrame(outer_rows)
    predictions = pd.DataFrame(prediction_rows)

    ensemble = predictions.loc[
        predictions.seed == "ensemble"
    ]

    if ensemble.duplicated(
        ["panel_id", "mode", "patient_uid"]
    ).any():
        raise RuntimeError("Duplicate ensemble predictions")

    if not (
        ensemble.groupby(["panel_id", "mode"]).size() == 67
    ).all():
        raise RuntimeError("Incomplete ensemble predictions")

    summaries = []

    for (panel, mode), rows in metrics.groupby(
        ["panel_id", "mode"]
    ):
        if len(rows) != 3:
            raise RuntimeError("Incomplete outer metrics")

        macro = float(rows.AUC.mean())
        minimum = float(rows.AUC.min())

        record = {
            "panel_id": panel,
            "mode": mode,
            "macro_AUC": macro,
            "minimum_cohort_AUC": minimum,
            "macro_average_precision":
                float(rows.average_precision.mean()),
            "target_met_development_only":
                bool(macro > 0.8 and minimum > 0.7),
        }

        for cohort in COHORTS:
            record[f"AUC_{cohort}"] = float(
                rows.loc[
                    rows.cohort == cohort, "AUC"
                ].iloc[0]
            )

        summaries.append(record)

    summary = pd.DataFrame(summaries).sort_values(
        [
            "target_met_development_only",
            "minimum_cohort_AUC",
            "macro_AUC",
        ],
        ascending=[False, False, False],
    )

    summary.to_csv(
        root / "fusion_method_summary.tsv",
        sep="\t", index=False,
    )

    print("\nFUSION_METHOD_SUMMARY:", flush=True)
    print(summary.to_string(index=False))
    print(
        f"\nTOTAL_SECONDS: {time.time()-started:.1f}",
        flush=True,
    )
    print("FIXED_PANEL_SSL_FUSION_COMPLETED", flush=True)


if __name__ == "__main__":
    main()

