
import argparse
import csv
import importlib.util
import json
import math
import random
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F


def import_model(path):
    spec = importlib.util.spec_from_file_location("ssl230_model", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--initialization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if (args.output / "training_log.csv").exists():
        raise RuntimeError("Training output exists; use a new directory")

    torch.set_num_threads(min(8, torch.get_num_threads()))
    torch.manual_seed(20261116)
    np.random.seed(20261116)
    random.seed(20261116)
    torch.use_deterministic_algorithms(True)

    module = import_model(args.model)
    arrays, meta, membership, tc, cm, model_spec = module.load_inputs(
        args.data
    )

    # Own project checkpoint; generated in the preceding validated step.
    initial = torch.load(
        args.initialization, map_location="cpu", weights_only=False
    )
    assert initial["optimizer_steps"] == 0
    assert initial["spec"] == model_spec
    assert initial["data_manifest_md5"] == module.md5(
        args.data / "export_manifest.tsv"
    )

    model = module.SSL37(membership, tc, cm, dim=24, heads=4)
    model.load_state_dict(initial["state_dict"], strict=True)

    x = torch.from_numpy(arrays["ssl"]["x"].copy())
    ranks = torch.from_numpy(arrays["ssl"]["rank"].copy())
    observed = torch.from_numpy(arrays["ssl"]["observed"].copy())

    train_rows = np.flatnonzero(meta["split"].to_numpy() == "train")
    val_rows = np.flatnonzero(meta["split"].to_numpy() == "validation")
    assert len(train_rows) == 554 and len(val_rows) == 62
    assert set(train_rows).isdisjoint(val_rows)
    assert observed.all()

    cohorts = ["TCGA", "ICGC", "CPTAC"]
    train_pools = {
        cohort: np.flatnonzero(
            (meta["split"].to_numpy() == "train")
            & (meta.cohort.to_numpy() == cohort)
        )
        for cohort in cohorts
    }
    val_pools = {
        cohort: np.flatnonzero(
            (meta["split"].to_numpy() == "validation")
            & (meta.cohort.to_numpy() == cohort)
        )
        for cohort in cohorts
    }
    assert all(len(v) > 0 for v in train_pools.values())
    assert all(len(v) > 0 for v in val_pools.values())

    config = {
        "seed": 20261116,
        "device": "cpu",
        "batch_size": 16,
        "max_epochs": 60,
        "patience": 12,
        "minimum_improvement": 1e-5,
        "learning_rate": 0.001,
        "weight_decay": 0.02,
        "gradient_clip_norm": 1.0,
        "mask_rate": 0.15,
        "rank_weight": 0.5,
        "consistency_weight": 0.05,
        "training_sampling": "rotating_equal_cohort_batches_with_replacement",
        "selection_metric": "macro_mean_validation_composite_loss_by_cohort",
        "validation_views": "fixed_two_views_per_validation_patient",
        "response_labels_used": False,
        "data_manifest_md5": module.md5(args.data / "export_manifest.tsv"),
        "model_script_md5": module.md5(args.model),
        "initialization_md5": module.md5(args.initialization),
        "torch_version": torch.__version__,
        "numpy_version": np.__version__,
    }

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "training_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    meta.to_csv(
        args.output / "SSL_fixed_split.tsv", sep="\t", index=False
    )

    # Fixed validation views generated once, shared across all epochs.
    fixed_validation = {}
    validation_arrays = {}

    for ci, cohort in enumerate(cohorts):
        rows = val_pools[cohort]
        generator = torch.Generator().manual_seed(900000 + ci)

        a, am, reconstruction_mask = module.make_view(
            x[rows], observed[rows], generator, rate=0.15
        )
        b, bm, _ = module.make_view(
            x[rows], observed[rows], generator, rate=0.15
        )
        fixed_validation[cohort] = (
            rows, a, am, reconstruction_mask, b, bm
        )
        validation_arrays[f"{cohort}_rows"] = rows
        validation_arrays[f"{cohort}_mask1"] = am.numpy()
        validation_arrays[f"{cohort}_mask2"] = bm.numpy()

    np.savez_compressed(
        args.output / "fixed_validation_masks.npz",
        **validation_arrays
    )

    # Baselines are fitted on SSL training patients only.
    mean_gene = x[train_rows].mean(dim=0)
    mean_rank = ranks[train_rows].mean(dim=0)

    def evaluate():
        model.eval()
        cohort_metrics = []

        with torch.no_grad():
            for cohort in cohorts:
                rows, a, am, masked, b, bm = fixed_validation[cohort]

                # Aggregate whole-cohort errors, not an average of batch MSEs.
                gene_se = 0.0
                gene_count = 0
                rank_se = 0.0
                rank_count = 0
                consistency_se = 0.0
                consistency_count = 0
                baseline_gene_se = 0.0
                baseline_rank_se = 0.0

                for start in range(0, len(rows), 16):
                    end = min(start + 16, len(rows))
                    rr = rows[start:end]

                    first = model(a[start:end], am[start:end])
                    _, _, second_modules = model.encode(
                        b[start:end], bm[start:end]
                    )
                    mask = masked[start:end]

                    errors = (first["gene_prediction"] - x[rr]).square()
                    gene_se += float(errors[mask].sum())
                    gene_count += int(mask.sum())

                    rank_errors = (
                        first["rank_prediction"] - ranks[rr]
                    ).square()
                    rank_se += float(rank_errors.sum())
                    rank_count += rank_errors.numel()

                    consistency_errors = (
                        first["modules"] - second_modules
                    ).square()
                    consistency_se += float(consistency_errors.sum())
                    consistency_count += consistency_errors.numel()

                    baseline_errors = (x[rr] - mean_gene).square()
                    baseline_gene_se += float(
                        baseline_errors[mask].sum()
                    )
                    baseline_rank_se += float(
                        (ranks[rr] - mean_rank).square().sum()
                    )

                gene_loss = gene_se / gene_count
                rank_loss = rank_se / rank_count
                consistency = consistency_se / consistency_count

                cohort_metrics.append({
                    "cohort": cohort,
                    "total": gene_loss + 0.5*rank_loss + 0.05*consistency,
                    "gene": gene_loss,
                    "rank": rank_loss,
                    "consistency": consistency,
                    "mean_predictor_gene": baseline_gene_se/gene_count,
                    "mean_predictor_rank": baseline_rank_se/rank_count,
                })

        macro = {
            key: float(np.mean([row[key] for row in cohort_metrics]))
            for key in cohort_metrics[0] if key != "cohort"
        }
        return macro, cohort_metrics

    def representation_diagnostics():
        model.eval()
        token_parts = []
        module_parts = []

        with torch.no_grad():
            for start in range(0, len(val_rows), 16):
                rows = val_rows[start:start+16]
                tokens, _, modules = model.encode(x[rows], observed[rows])
                token_parts.append(tokens.numpy())
                module_parts.append(modules.flatten(1).numpy())

        token_values = np.concatenate(token_parts)
        module_values = np.concatenate(module_parts)

        normalized = module_values / np.maximum(
            np.linalg.norm(module_values, axis=1, keepdims=True), 1e-12
        )
        centered = normalized - normalized.mean(axis=0, keepdims=True)
        eigenvalues = np.linalg.svd(
            centered, compute_uv=False
        ) ** 2

        if eigenvalues.sum() <= 1e-20:
            rank = 0.0
        else:
            p = eigenvalues / eigenvalues.sum()
            p = p[p > 0]
            rank = float(np.exp(-(p*np.log(p)).sum()))

        return {
            "validation_token_std_mean":
                float(token_values.std(axis=0).mean()),
            "validation_module_effective_rank": rank,
            "validation_module_normalized_std_mean":
                float(normalized.std(axis=0).mean()),
        }

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config["learning_rate"],
        weight_decay=config["weight_decay"]
    )
    sampling_rng = np.random.default_rng(config["seed"])
    augmentation_rng = torch.Generator().manual_seed(
        config["seed"] + 1
    )

    initial_metrics, initial_by_cohort = evaluate()
    best_loss = initial_metrics["total"]
    best_epoch = 0
    stale = 0

    def checkpoint(epoch, metrics):
        return {
            "state_dict": model.state_dict(),
            "spec": model_spec,
            "epoch": epoch,
            "validation_metrics": metrics,
            "optimizer_state": optimizer.state_dict(),
            "config": config,
            "torch_rng_state": torch.get_rng_state(),
            "augmentation_rng_state": augmentation_rng.get_state(),
            "sampling_rng_state": sampling_rng.bit_generator.state,
            "response_labels_used": False,
        }

    torch.save(
        checkpoint(0, initial_metrics),
        args.output / "SSL37_best.pt"
    )

    print("INITIAL_VALIDATION:", json.dumps(initial_metrics), flush=True)

    fields = [
        "epoch", "train_total", "train_gene", "train_rank",
        "train_consistency", "val_total", "val_gene", "val_rank",
        "val_consistency", "mean_predictor_gene",
        "mean_predictor_rank", "epoch_seconds", "best_epoch",
        "stale_epochs"
    ]
    steps_per_epoch = math.ceil(len(train_rows)/config["batch_size"])

    with open(
        args.output / "training_log.csv", "w",
        newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow({
            "epoch": 0,
            **{f"val_{k}": initial_metrics[k] for k in
               ["total", "gene", "rank", "consistency"]},
            "mean_predictor_gene": initial_metrics["mean_predictor_gene"],
            "mean_predictor_rank": initial_metrics["mean_predictor_rank"],
            "epoch_seconds": 0,
            "best_epoch": 0,
            "stale_epochs": 0,
        })
        handle.flush()

        for epoch in range(1, config["max_epochs"]+1):
            start_time = time.perf_counter()
            model.train()
            training_values = []

            for step in range(steps_per_epoch):
                cohort = cohorts[(step + epoch - 1) % len(cohorts)]
                rows = sampling_rng.choice(
                    train_pools[cohort],
                    size=config["batch_size"], replace=True
                )

                optimizer.zero_grad(set_to_none=True)
                loss, values = module.ssl_loss(
                    model, x[rows], ranks[rows],
                    observed[rows], augmentation_rng
                )

                if not torch.isfinite(loss).item():
                    raise RuntimeError(f"Nonfinite loss at epoch {epoch}")

                loss.backward()
                torch.nn.utils.clip_grad_norm_(
                    model.parameters(), 1.0, error_if_nonfinite=True
                )
                optimizer.step()

                training_values.append({
                    name: float(value.detach())
                    for name, value in values.items()
                })

            metrics, by_cohort = evaluate()
            assert all(np.isfinite(v) for v in metrics.values())

            if metrics["total"] < best_loss - config["minimum_improvement"]:
                best_loss = metrics["total"]
                best_epoch = epoch
                stale = 0
                torch.save(
                    checkpoint(epoch, metrics),
                    args.output / "SSL37_best.pt"
                )
            else:
                stale += 1

            torch.save(
                checkpoint(epoch, metrics),
                args.output / "SSL37_last.pt"
            )

            means = {
                key: float(np.mean([row[key] for row in training_values]))
                for key in training_values[0]
            }
            elapsed = time.perf_counter() - start_time

            writer.writerow({
                "epoch": epoch,
                "train_total": means["total"],
                "train_gene": means["gene_reconstruction"],
                "train_rank": means["rank_reconstruction"],
                "train_consistency": means["module_consistency"],
                **{f"val_{k}": metrics[k] for k in
                   ["total", "gene", "rank", "consistency"]},
                "mean_predictor_gene": metrics["mean_predictor_gene"],
                "mean_predictor_rank": metrics["mean_predictor_rank"],
                "epoch_seconds": elapsed,
                "best_epoch": best_epoch,
                "stale_epochs": stale,
            })
            handle.flush()

            print(
                f"EPOCH {epoch:03d} "
                f"train={means['total']:.5f} "
                f"val={metrics['total']:.5f} "
                f"gene={metrics['gene']:.5f} "
                f"rank={metrics['rank']:.5f} "
                f"best={best_epoch} seconds={elapsed:.1f}",
                flush=True
            )

            if stale >= config["patience"]:
                print("EARLY_STOP", flush=True)
                break

    # Always reload the best checkpoint for final diagnostics.
    best = torch.load(
        args.output / "SSL37_best.pt",
        map_location="cpu", weights_only=False
    )
    model.load_state_dict(best["state_dict"], strict=True)

    best_metrics, best_by_cohort = evaluate()
    diagnostics = representation_diagnostics()

    pd.DataFrame(best_by_cohort).to_csv(
        args.output / "best_validation_by_cohort.tsv",
        sep="\t", index=False
    )

    summary = {
        "best_epoch": best_epoch,
        "epochs_run": epoch,
        "initial_validation": initial_metrics,
        "best_validation": best_metrics,
        "representation_diagnostics": diagnostics,
        "gene_reconstruction_better_than_mean_predictor":
            best_metrics["gene"] < best_metrics["mean_predictor_gene"],
        "rank_reconstruction_better_than_mean_predictor":
            best_metrics["rank"] < best_metrics["mean_predictor_rank"],
        "improved_over_initialization": best_epoch > 0,
        "best_checkpoint_md5":
            module.md5(args.output / "SSL37_best.pt"),
        "response_labels_used": False,
        "ICB_performance_evaluated": False,
        "completed": True,
    }

    (args.output / "training_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print("TRAINING_SUMMARY:", json.dumps(summary), flush=True)
    print("SSL37_PRETRAINING_COMPLETED", flush=True)


if __name__ == "__main__":
    main()

