#!/usr/bin/env python3
"""Reproduce the paper QC and run a first-pass mixed-effects cosinor analysis.

The primary analysis intentionally does not impute missing protein intensities.
It removes the eight low-quality samples using the paper's protein-count rule,
and treats timestamp-suffixed files as technical re-injections rather than new
biological time points.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests


SAMPLE_RE = re.compile(
    r"Plate(?P<plate>\d+)_(?P<subject>\d+)-Dag(?P<day>\d+)-t(?P<hour>\d+)"
    r"(?:_(?P<reinjection>\d+))?\.mzML$",
    re.IGNORECASE,
)
PLATE_RE = re.compile(r"Plate(?P<plate>\d+)_(?P<label>.+)\.mzML$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-top", type=int, default=10)
    parser.add_argument("--bootstrap-reps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20250901)
    return parser.parse_args()


def sample_manifest(columns: list[str], matrix: pd.DataFrame) -> pd.DataFrame:
    records = []
    for column in columns:
        basename = column.replace("\\", "/").split("/")[-1]
        sample_match = SAMPLE_RE.search(basename)
        plate_match = PLATE_RE.search(basename)
        is_empty = bool(plate_match and plate_match.group("label").lower().startswith("empty"))
        values = pd.to_numeric(matrix[column], errors="coerce")
        record = {
            "sample_column": column,
            "basename": basename,
            "plate": int(plate_match.group("plate")) if plate_match else np.nan,
            "subject_id": sample_match.group("subject") if sample_match else "",
            "day": int(sample_match.group("day")) if sample_match else np.nan,
            "clock_hour": int(sample_match.group("hour")) if sample_match else np.nan,
            "reinjection_timestamp": sample_match.group("reinjection") if sample_match else "",
            "is_empty": is_empty,
            "protein_count": int(values.notna().sum()),
            "missing_fraction": float(values.isna().mean()),
            "median_log2_intensity": float(np.log2(values.dropna()).median())
            if values.notna().any()
            else np.nan,
        }
        record["is_technical_reinjection"] = bool(record["reinjection_timestamp"])
        record["passes_paper_sample_qc"] = bool(
            not is_empty and record["protein_count"] >= 457
        )
        record["include_primary"] = bool(
            record["passes_paper_sample_qc"] and not record["is_technical_reinjection"]
        )
        reasons = []
        if is_empty:
            reasons.append("empty_run")
        if not is_empty and record["protein_count"] < 457:
            reasons.append("paper_low_protein_count")
        if record["is_technical_reinjection"]:
            reasons.append("timestamped_reinjection")
        record["exclusion_reason"] = ";".join(reasons)
        records.append(record)
    manifest = pd.DataFrame(records)
    manifest["elapsed_hour"] = np.where(
        manifest["day"].eq(1), manifest["clock_hour"] - 9, manifest["clock_hour"] + 15
    )
    return manifest


def paper_protein_qc(
    raw: pd.DataFrame, manifest: pd.DataFrame, annotation_columns: list[str]
) -> pd.DataFrame:
    paper_columns = manifest.loc[manifest["passes_paper_sample_qc"], "sample_column"].tolist()
    values = raw[paper_columns].apply(pd.to_numeric, errors="coerce")
    result = raw[annotation_columns].copy()
    result["observed_paper_qc"] = values.notna().sum(axis=1)
    result["missing_fraction_paper_qc"] = values.isna().mean(axis=1)
    result["passes_paper_protein_qc"] = result["missing_fraction_paper_qc"] <= 0.40
    primary_columns = manifest.loc[manifest["include_primary"], "sample_column"].tolist()
    primary = raw[primary_columns].apply(pd.to_numeric, errors="coerce")
    result["observed_primary"] = primary.notna().sum(axis=1)
    result["missing_fraction_primary"] = primary.isna().mean(axis=1)
    empty_columns = manifest.loc[manifest["is_empty"], "sample_column"].tolist()
    empty = raw[empty_columns].apply(pd.to_numeric, errors="coerce")
    result["empty_detection_rate"] = empty.notna().mean(axis=1)
    result["empty_median_intensity"] = empty.median(axis=1, skipna=True)
    result["primary_median_intensity"] = primary.median(axis=1, skipna=True)
    result["empty_to_primary_median_ratio"] = (
        result["empty_median_intensity"] / result["primary_median_intensity"]
    )
    result["carryover_flag"] = (
        result["empty_detection_rate"].ge(0.5)
        & result["empty_to_primary_median_ratio"].ge(0.1)
    )
    return result


def design_for_protein(
    values: pd.Series, manifest: pd.DataFrame
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    meta = manifest.loc[manifest["include_primary"]].copy().set_index("sample_column")
    meta["y"] = np.log2(pd.to_numeric(values.loc[meta.index], errors="coerce"))
    meta = meta.replace([np.inf, -np.inf], np.nan).dropna(subset=["y"])
    angle = 2 * np.pi * meta["elapsed_hour"].to_numpy(dtype=float) / 24.0
    meta["cos24"] = np.cos(angle)
    meta["sin24"] = np.sin(angle)
    plate = pd.get_dummies(meta["plate"].astype(int), prefix="plate", drop_first=True, dtype=float)
    exog = np.column_stack(
        [np.ones(len(meta)), meta["cos24"], meta["sin24"], plate.to_numpy()]
    )
    names = np.array(["intercept", "cos24", "sin24", *plate.columns.tolist()])
    exog_re_null = np.ones((len(meta), 1))
    exog_re_full = np.column_stack(
        [np.ones(len(meta)), meta["cos24"].to_numpy(), meta["sin24"].to_numpy()]
    )
    return meta, exog, names, exog_re_null, exog_re_full


def fit_mixedlm(endog, exog, groups, exog_re, maxiter: int):
    """Try multiple optimizers and retain the best converged ML solution."""
    fits = []
    errors = []
    for method in ("lbfgs", "bfgs", "cg"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = MixedLM(endog, exog, groups=groups, exog_re=exog_re).fit(
                    reml=False, method=method, maxiter=maxiter, disp=False
                )
            if np.isfinite(fit.llf):
                fits.append((fit, method))
            if fit.converged:
                break
        except Exception as exc:
            errors.append(f"{method}:{type(exc).__name__}")
    if not fits:
        raise RuntimeError("all optimizers failed: " + ",".join(errors))
    converged = [item for item in fits if item[0].converged]
    pool = converged if converged else fits
    return max(pool, key=lambda item: item[0].llf)


def fit_one(values: pd.Series, manifest: pd.DataFrame) -> tuple[dict, dict | None]:
    meta, exog, names, re_null, re_full = design_for_protein(values, manifest)
    base = {
        "n_obs": len(meta),
        "n_subjects": meta["subject_id"].nunique(),
        "converged_null": False,
        "converged_full": False,
        "singular_full": False,
        "optimizer_null": "",
        "optimizer_full": "",
        "fit_error": "",
    }
    if len(meta) < 100 or meta["subject_id"].nunique() < 20:
        base["fit_error"] = "insufficient_observations"
        return base, None
    try:
        null, null_method = fit_mixedlm(meta["y"], exog, meta["subject_id"], re_null, 600)
        full, full_method = fit_mixedlm(meta["y"], exog, meta["subject_id"], re_full, 1200)
        base["converged_null"] = bool(null.converged)
        base["converged_full"] = bool(full.converged)
        base["optimizer_null"] = null_method
        base["optimizer_full"] = full_method
        eig = np.linalg.eigvalsh(np.asarray(full.cov_re))
        base["singular_full"] = bool(np.min(eig) < 1e-6)
        lrt = max(0.0, 2.0 * (full.llf - null.llf))
        beta = np.asarray(full.fe_params)
        cov = np.asarray(full.cov_params())[: len(beta), : len(beta)]
        rhythm_beta = beta[1:3]
        rhythm_cov = cov[1:3, 1:3]
        wald = float(rhythm_beta @ np.linalg.pinv(rhythm_cov) @ rhythm_beta)
        amplitude = float(np.hypot(beta[1], beta[2]))
        peak_elapsed = float((math.atan2(beta[2], beta[1]) * 24 / (2 * np.pi)) % 24)
        peak_clock = float((peak_elapsed + 9) % 24)
        cov_re = np.asarray(full.cov_re)
        base.update(
            {
                "ll_null": float(null.llf),
                "ll_full": float(full.llf),
                "lrt_stat": lrt,
                "lrt_chi2_df5_p": float(chi2.sf(lrt, 5)),
                "population_wald_stat": wald,
                "population_rhythm_p": float(chi2.sf(wald, 2)),
                "population_amplitude_log2": amplitude,
                "population_peak_clock": peak_clock,
                "random_intercept_variance": float(cov_re[0, 0]),
                "random_cos_variance": float(cov_re[1, 1]),
                "random_sin_variance": float(cov_re[2, 2]),
                "residual_variance": float(full.scale),
            }
        )
        cache = {"meta": meta, "exog": exog, "re_null": re_null, "re_full": re_full, "null": null}
        return base, cache
    except Exception as exc:  # individual proteins must not abort the screen
        base["fit_error"] = f"{type(exc).__name__}: {exc}"
        return base, None


def bootstrap_lrt(cache: dict, observed_lrt: float, reps: int, rng: np.random.Generator) -> dict:
    meta, exog = cache["meta"], cache["exog"]
    re_null, re_full, null = cache["re_null"], cache["re_full"], cache["null"]
    subjects = meta["subject_id"].astype(str).to_numpy()
    unique_subjects = pd.unique(subjects)
    fixed_mean = exog @ np.asarray(null.fe_params)
    tau = math.sqrt(max(float(np.asarray(null.cov_re)[0, 0]), 0.0))
    sigma = math.sqrt(max(float(null.scale), 0.0))
    simulated = []
    failed = 0
    for _ in range(reps):
        offsets = dict(zip(unique_subjects, rng.normal(0, tau, len(unique_subjects))))
        y = fixed_mean + np.array([offsets[x] for x in subjects]) + rng.normal(0, sigma, len(meta))
        try:
            nfit, _ = fit_mixedlm(y, exog, subjects, re_null, 400)
            ffit, _ = fit_mixedlm(y, exog, subjects, re_full, 700)
            if np.isfinite(nfit.llf) and np.isfinite(ffit.llf):
                simulated.append(max(0.0, 2 * (ffit.llf - nfit.llf)))
            else:
                failed += 1
        except Exception:
            failed += 1
    exceed = sum(x >= observed_lrt for x in simulated)
    return {
        "bootstrap_p": (exceed + 1) / (len(simulated) + 1),
        "bootstrap_success": len(simulated),
        "bootstrap_failed": failed,
    }


def plot_qc(manifest: pd.DataFrame, output: Path) -> None:
    real = manifest.loc[~manifest["is_empty"]]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].hist(real["protein_count"], bins=24, color="#4472C4", edgecolor="white")
    axes[0].axvline(457, color="#C00000", linestyle="--", label="paper cutoff = 457")
    axes[0].set(xlabel="Detected protein groups", ylabel="Samples", title="Sample QC")
    axes[0].legend(frameon=False)
    for plate, group in real.groupby("plate"):
        axes[1].scatter(group["protein_count"], group["median_log2_intensity"], label=f"Plate {plate}")
    axes[1].axvline(457, color="#C00000", linestyle="--")
    axes[1].set(xlabel="Detected protein groups", ylabel="Median log2 intensity", title="Plate structure")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(args.input, sep="\t")
    annotations = raw.columns[:4].tolist()
    manifest = sample_manifest(raw.columns[4:].tolist(), raw)
    protein_qc = paper_protein_qc(raw, manifest, annotations)
    manifest.to_csv(args.output_dir / "sample_manifest.tsv", sep="\t", index=False)
    protein_qc.to_csv(args.output_dir / "protein_qc.tsv", sep="\t", index=False)
    plot_qc(manifest, args.output_dir / "sample_qc.png")

    results = []
    caches: dict[int, dict] = {}
    selected = protein_qc.index[protein_qc["passes_paper_protein_qc"]]
    for number, index in enumerate(selected, start=1):
        fit, cache = fit_one(raw.loc[index, manifest["sample_column"]], manifest)
        fit.update({name: raw.loc[index, name] for name in annotations})
        for name in (
            "missing_fraction_primary",
            "empty_detection_rate",
            "empty_to_primary_median_ratio",
            "carryover_flag",
        ):
            fit[name] = protein_qc.loc[index, name]
        fit["row_index"] = int(index)
        results.append(fit)
        if cache is not None:
            caches[int(index)] = cache
        if number % 50 == 0:
            print(f"fitted {number}/{len(selected)} proteins", flush=True)
    result = pd.DataFrame(results)
    valid_model = (
        result["converged_null"]
        & result["converged_full"]
        & ~result["singular_full"]
        & result["fit_error"].fillna("").eq("")
    )
    valid_population = valid_model & result["population_rhythm_p"].notna()
    result["population_rhythm_q"] = np.nan
    if valid_population.any():
        result.loc[valid_population, "population_rhythm_q"] = multipletests(
            result.loc[valid_population, "population_rhythm_p"], method="fdr_bh"
        )[1]
    valid_lrt = valid_model & result["lrt_chi2_df5_p"].notna()
    result["heterogeneity_chi2_q"] = np.nan
    if valid_lrt.any():
        result.loc[valid_lrt, "heterogeneity_chi2_q"] = multipletests(
            result.loc[valid_lrt, "lrt_chi2_df5_p"], method="fdr_bh"
        )[1]
    result["bootstrap_p"] = np.nan
    result["bootstrap_success"] = np.nan
    result["bootstrap_failed"] = np.nan
    candidates = result.loc[
        valid_model & result["lrt_stat"].notna()
    ].nsmallest(args.bootstrap_top, "lrt_chi2_df5_p")
    rng = np.random.default_rng(args.seed)
    for index in candidates.index:
        row_index = int(result.loc[index, "row_index"])
        boot = bootstrap_lrt(
            caches[row_index], float(result.loc[index, "lrt_stat"]), args.bootstrap_reps, rng
        )
        for name, value in boot.items():
            result.loc[index, name] = value
        print(f"bootstrapped {result.loc[index, 'Genes']}: {boot}", flush=True)
    result.to_csv(args.output_dir / "mixed_cosinor_results.tsv", sep="\t", index=False)

    summary = {
        "input": str(args.input),
        "raw_protein_groups": len(raw),
        "raw_real_samples": int((~manifest["is_empty"]).sum()),
        "empty_runs": int(manifest["is_empty"].sum()),
        "paper_low_quality_samples": int(
            ((~manifest["is_empty"]) & ~manifest["passes_paper_sample_qc"]).sum()
        ),
        "paper_qc_samples": int(manifest["passes_paper_sample_qc"].sum()),
        "timestamped_reinjections_excluded_primary": int(
            manifest["is_technical_reinjection"].sum()
        ),
        "primary_sample_columns": int(manifest["include_primary"].sum()),
        "paper_qc_proteins": int(protein_qc["passes_paper_protein_qc"].sum()),
        "models_attempted": len(result),
        "full_models_converged": int(result["converged_full"].sum()),
        "full_models_nonsingular": int(
            (result["converged_full"] & ~result["singular_full"]).sum()
        ),
        "population_rhythmic_fdr_0_05": int((result["population_rhythm_q"] < 0.05).sum()),
        "heterogeneity_chi2_fdr_0_05_diagnostic": int(
            (result["heterogeneity_chi2_q"] < 0.05).sum()
        ),
        "bootstrap_top": args.bootstrap_top,
        "bootstrap_reps": args.bootstrap_reps,
        "seed": args.seed,
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
