#!/usr/bin/env python3
"""Internal robustness checks for simplified mixed-cosinor candidates."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
from first_pass_mixed_cosinor import sample_manifest  # noqa: E402
from simplified_mixed_cosinor import bootstrap, fit_formula, fit_one, protein_data  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--screen-results", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-reps", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20250903)
    return parser.parse_args()


def compare(data: pd.DataFrame, fixed_formula: str | None = None) -> dict:
    kwargs = {} if fixed_formula is None else {"fixed_formula": fixed_formula}
    try:
        null, _ = fit_formula(data, full=False, maxiter=500, **kwargs)
        full, _ = fit_formula(data, full=True, maxiter=900, **kwargs)
        lrt = max(0.0, 2 * (full.llf - null.llf))
        return {
            "converged": bool(null.converged and full.converged),
            "lrt": lrt,
            "p": float(chi2.sf(lrt, 2)),
            "boundary": bool(np.any(np.asarray(full.vcomp) < 1e-6)),
        }
    except Exception:
        return {"converged": False, "lrt": np.nan, "p": np.nan, "boundary": True}


def largest_individual_amplitude(data: pd.DataFrame) -> tuple[str, float]:
    estimates = []
    for subject, group in data.groupby("subject_id"):
        observed = group.dropna(subset=["y"])
        if len(observed) < 5:
            continue
        design = np.column_stack(
            [np.ones(len(observed)), observed["cos24"], observed["sin24"]]
        )
        coef = np.linalg.lstsq(design, observed["y"].to_numpy(), rcond=None)[0]
        estimates.append((subject, float(np.hypot(coef[1], coef[2]))))
    return max(estimates, key=lambda item: item[1])


def validate_candidate(
    gene: str,
    row_index: int,
    matrix: pd.DataFrame,
    manifest: pd.DataFrame,
    bootstrap_reps: int,
    seed: int,
) -> tuple[dict, list[dict], list[dict]]:
    data = protein_data(matrix.loc[row_index, manifest["sample_column"]], manifest)
    original_fit, cache = fit_one(matrix.loc[row_index, manifest["sample_column"]], manifest)
    largest_subject, largest_amplitude = largest_individual_amplitude(data)

    loo_rows = []
    for subject in sorted(data["subject_id"].unique(), key=int):
        check = compare(data.loc[data["subject_id"].ne(subject)].copy())
        loo_rows.append({"gene": gene, "excluded_subject": subject, **check})

    origin_rows = []
    for shift in (0, 3, 6, 9, 12):
        shifted = data.copy()
        angle = 2 * np.pi * (shifted["elapsed_hour"].to_numpy() - shift) / 24
        shifted["cos24"] = np.cos(angle)
        shifted["sin24"] = np.sin(angle)
        check = compare(shifted)
        origin_rows.append({"gene": gene, "origin_shift_hours": shift, **check})

    trend = data.copy()
    trend["elapsed_z"] = (trend["elapsed_hour"] - trend["elapsed_hour"].mean()) / trend[
        "elapsed_hour"
    ].std()
    trend_check = compare(
        trend, fixed_formula="y ~ cos24 + sin24 + elapsed_z + C(plate)"
    )
    boot = bootstrap(cache, float(original_fit["lrt_stat"]), bootstrap_reps, np.random.default_rng(seed))

    loo = pd.DataFrame(loo_rows)
    origins = pd.DataFrame(origin_rows)
    summary = {
        "gene": gene,
        "row_index": row_index,
        "original_lrt": original_fit["lrt_stat"],
        "original_chi2_p": original_fit["lrt_chi2_df2_p"],
        "largest_individual_amplitude_subject": largest_subject,
        "largest_individual_amplitude_log2": largest_amplitude,
        "loo_success": int(loo["converged"].sum()),
        "loo_total": len(loo),
        "loo_p_lt_0_05_fraction": float((loo.loc[loo["converged"], "p"] < 0.05).mean()),
        "loo_min_lrt": float(loo.loc[loo["converged"], "lrt"].min()),
        "loo_max_p": float(loo.loc[loo["converged"], "p"].max()),
        "loo_without_largest_subject_p": float(
            loo.loc[loo["excluded_subject"].eq(largest_subject), "p"].iloc[0]
        ),
        "origin_success": int(origins["converged"].sum()),
        "origin_p_lt_0_05_fraction": float(
            (origins.loc[origins["converged"], "p"] < 0.05).mean()
        ),
        "origin_max_p": float(origins.loc[origins["converged"], "p"].max()),
        "linear_trend_converged": trend_check["converged"],
        "linear_trend_lrt": trend_check["lrt"],
        "linear_trend_p": trend_check["p"],
        **boot,
    }
    return summary, loo_rows, origin_rows


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix, sep="\t")
    manifest = sample_manifest(matrix.columns[4:].tolist(), matrix)
    screen = pd.read_csv(args.screen_results, sep="\t")
    candidates = screen.loc[screen["heterogeneity_chi2_q"] < 0.05].sort_values(
        "heterogeneity_chi2_q"
    )
    summaries, loo_rows, origin_rows = [], [], []
    for number, (_, candidate) in enumerate(candidates.iterrows(), start=1):
        gene = str(candidate["Genes"])
        summary, loo, origins = validate_candidate(
            gene,
            int(candidate["row_index"]),
            matrix,
            manifest,
            args.bootstrap_reps,
            args.seed + number,
        )
        summary["screen_q"] = float(candidate["heterogeneity_chi2_q"])
        summary["carryover_flag"] = bool(candidate["carryover_flag"])
        summary["boundary_full"] = bool(candidate["boundary_full"])
        summaries.append(summary)
        loo_rows.extend(loo)
        origin_rows.extend(origins)
        print(f"validated {number}/{len(candidates)}: {gene}", flush=True)

    summary_frame = pd.DataFrame(summaries)
    summary_frame["passes_internal_robustness"] = (
        summary_frame["loo_p_lt_0_05_fraction"].eq(1.0)
        & summary_frame["origin_p_lt_0_05_fraction"].eq(1.0)
        & summary_frame["linear_trend_p"].lt(0.05)
        & summary_frame["bootstrap_p"].lt(0.05)
        & summary_frame["bootstrap_success"].ge(0.90 * args.bootstrap_reps)
        & ~summary_frame["boundary_full"]
    )
    summary_frame["high_priority_clean"] = (
        summary_frame["passes_internal_robustness"] & ~summary_frame["carryover_flag"]
    )
    summary_frame.to_csv(args.output_dir / "candidate_validation_summary.tsv", sep="\t", index=False)
    pd.DataFrame(loo_rows).to_csv(
        args.output_dir / "leave_one_subject_out.tsv", sep="\t", index=False
    )
    pd.DataFrame(origin_rows).to_csv(
        args.output_dir / "time_origin_sensitivity.tsv", sep="\t", index=False
    )
    run_summary = {
        "candidates": len(summary_frame),
        "bootstrap_reps_per_candidate": args.bootstrap_reps,
        "passes_internal_robustness": int(summary_frame["passes_internal_robustness"].sum()),
        "high_priority_clean": int(summary_frame["high_priority_clean"].sum()),
        "criteria": {
            "all_24_leave_one_subject_out_nominal_p_lt_0_05": True,
            "all_time_origins_nominal_p_lt_0_05": True,
            "linear_trend_adjusted_nominal_p_lt_0_05": True,
            "bootstrap_p_lt_0_05": True,
            "bootstrap_success_fraction": ">=0.90",
            "not_on_variance_boundary": True,
        },
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(run_summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(run_summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
