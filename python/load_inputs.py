
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_manifest(root, filename):
    manifest = pd.read_csv(root / filename, sep="\t")

    for row in manifest.itertuples(index=False):
        path = root / row.filename
        if not path.exists() or md5(path) != row.md5:
            raise RuntimeError(f"File fingerprint mismatch: {path}")


def load_data(root, context_dir):

    verify_manifest(root, "export_manifest.tsv")
    verify_manifest(context_dir, "context_files_manifest.tsv")

    provenance = json.loads(
        (context_dir / "context_provenance.json").read_text(
            encoding="utf-8"
        )
    )

    if provenance["data_manifest_md5"] != md5(
        root / "export_manifest.tsv"
    ):
        raise RuntimeError("Context and expression versions differ")

    for path_key, hash_key in [
        ("pretrained_checkpoint", "pretrained_checkpoint_md5"),
        ("initialization", "initialization_md5"),
        ("model_script", "model_script_md5"),
    ]:
        if md5(Path(provenance[path_key])) != provenance[hash_key]:
            raise RuntimeError(f"Source changed: {path_key}")

    hierarchy = pd.read_csv(root / "token_hierarchy.tsv", sep="\t")
    shapes = pd.read_csv(root / "shapes.tsv", sep="\t")

    if len(hierarchy) != 37 or not np.array_equal(
        hierarchy.token_index0.to_numpy(), np.arange(37)
    ):
        raise RuntimeError("Invalid 37-token hierarchy")

    names = hierarchy.token_name.tolist()

    if len(set(names)) != 37:
        raise RuntimeError("Duplicate token names")


    rank_parts = []
    context_parts = []
    metadata_parts = []

    for cohort in ["GSE181946", "GSE202069", "PRJEB34724"]:
        row = shapes.loc[shapes.dataset == cohort].iloc[0]

        n = int(row.n_samples)

        if (
            int(row.n_tokens) != 37
            or int(row.n_SSL_genes) != 1753
        ):
            raise RuntimeError("Input dimensions differ")

        metadata = pd.read_csv(
            root / f"{cohort}_metadata.tsv", sep="\t"
        )

        if (
            len(metadata) != n
            or not metadata.cohort.eq(cohort).all()
            or not metadata.response.isin([0, 1]).all()
            or metadata.response.nunique() != 2
        ):
            raise RuntimeError("Invalid clinical metadata")

        rank = np.fromfile(
            root / f"{cohort}_token_rank.float32.bin",
            dtype="<f4",
        ).reshape(n, 37)

        with np.load(
            context_dir / f"{cohort}_pretrained_context.npz",
            allow_pickle=False,
        ) as saved:
            if saved["token_names"].tolist() != names:
                raise RuntimeError("Context token order mismatch")

            if not np.array_equal(
                saved["patient_uid"].astype(str),
                metadata.patient_uid.to_numpy(str),
            ):
                raise RuntimeError("Context patient order mismatch")

            context = saved["token_context"].copy()

        if context.shape != (n, 37, 24):
            raise RuntimeError("Incorrect context dimensions")

        if not np.isfinite(context).all():
            raise RuntimeError("Nonfinite context")

        if not np.isfinite(rank).all() or (rank < 0).any() or (rank > 1).any():
            raise RuntimeError("Invalid token rank values")
        rank_parts.append(rank)
        context_parts.append(context)
        metadata_parts.append(metadata)

    metadata = pd.concat(metadata_parts, ignore_index=True)

    if (
        len(metadata) != 67
        or not metadata.patient_uid.is_unique
        or int(metadata.response.sum()) != 20
    ):
        raise RuntimeError("Unexpected patient manifest")

    print("LOCAL_INPUT_AND_PATIENT_CHECK: PASS", flush=True)
    print("BA37_CONTEXT_INPUT_CHECK: PASS", flush=True)

    return {
        "rank": np.concatenate(rank_parts),
        "ssl": np.concatenate(context_parts),
        "meta": metadata,
        "y": metadata.response.to_numpy(dtype=int),
        "cohort": metadata.cohort.to_numpy(str),
        "hierarchy": hierarchy,
        "provenance": provenance,
    }

