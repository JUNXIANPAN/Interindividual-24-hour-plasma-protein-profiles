#!/usr/bin/env python3
"""Build a hashed ScientificResultManifest for the protein-rhythm workflow."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--protocol-version", required=True)
    parser.add_argument("--input-matrix", type=Path, required=True)
    parser.add_argument("--screen-results", type=Path, required=True)
    parser.add_argument("--screen-summary", type=Path, required=True)
    parser.add_argument("--validation-results", type=Path, required=True)
    parser.add_argument("--validation-summary", type=Path, required=True)
    parser.add_argument("--fitted-plot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path, role: str) -> dict:
    return {
        "path": str(path),
        "role": role,
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
    }


def git_version() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except Exception:
        return "unknown"


def main() -> None:
    args = parse_args()
    screen_summary = json.loads(args.screen_summary.read_text(encoding="utf-8"))
    validation_summary = json.loads(args.validation_summary.read_text(encoding="utf-8"))
    validation = pd.read_csv(args.validation_results, sep="\t")
    robust = validation.loc[validation["passes_internal_robustness"], "gene"].tolist()
    clean = validation.loc[validation["high_priority_clean"], "gene"].tolist()
    manifest = {
        "schema_version": "1.0.0",
        "run_id": args.run_id,
        "job_id": args.job_id,
        "workflow_name": "protein_rhythm_internal_validation",
        "protocol_version": args.protocol_version,
        "status": "COMPLETED",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "code_version": git_version(),
        "inputs": [artifact(args.input_matrix, "raw_protein_matrix")],
        "outputs": [
            artifact(args.screen_results, "mixed_cosinor_screen"),
            artifact(args.screen_summary, "screen_summary"),
            artifact(args.validation_results, "candidate_validation"),
            artifact(args.validation_summary, "validation_summary"),
            artifact(args.fitted_plot, "candidate_plot"),
        ],
        "methods": {
            "sample_qc": "detected protein groups >= 457",
            "protein_qc": "missing fraction <= 0.40 after sample QC",
            "transformation": "log2",
            "imputation": "none",
            "null_model": "fixed cos24/sin24/plate + subject random intercept",
            "full_model": "null + independent subject random cosine and sine variances",
            "validation": [
                "leave-one-subject-out",
                "time-origin sensitivity",
                "linear elapsed-time sensitivity",
                "parametric bootstrap",
                "variance-boundary diagnostic",
                "Empty-run carryover flag",
            ],
        },
        "summary": {
            **screen_summary,
            **validation_summary,
            "internally_robust_candidates": robust,
            "high_priority_clean_candidates": clean,
            "evidence_level": "internally robust preliminary evidence",
            "independent_replication": False,
            "run_order_adjustment": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
