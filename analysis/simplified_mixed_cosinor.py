#!/usr/bin/env python3
"""Mixed cosinor with independent subject baseline/cosine/sine variances."""

from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from first_pass_mixed_cosinor import paper_protein_qc, sample_manifest  # noqa: E402


FIXED_FORMULA = "y ~ cos24 + sin24 + C(plate)"
VC_FORMULA = {"random_cos": "0 + cos24", "random_sin": "0 + sin24"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-top", type=int, default=10)
    parser.add_argument("--bootstrap-reps", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20250902)
    return parser.parse_args()


def protein_data(values: pd.Series, manifest: pd.DataFrame) -> pd.DataFrame:
    meta = manifest.loc[manifest["include_primary"]].copy().set_index("sample_column")
    meta["y"] = np.log2(pd.to_numeric(values.loc[meta.index], errors="coerce"))
    meta = meta.replace([np.inf, -np.inf], np.nan).dropna(subset=["y"])
    angle = 2 * np.pi * meta["elapsed_hour"].to_numpy(dtype=float) / 24
    meta["cos24"] = np.cos(angle)
    meta["sin24"] = np.sin(angle)
    meta["subject_id"] = meta["subject_id"].astype(str)
    meta["plate"] = meta["plate"].astype(int).astype(str)
    return meta.reset_index(drop=True)


def fit_formula(
    data: pd.DataFrame,
    full: bool,
    maxiter: int,
    fixed_formula: str = FIXED_FORMULA,
):
    fits = []
    failures = []
    for method in ("lbfgs", "bfgs", "cg", "powell"):
        try:
            model = MixedLM.from_formula(
                fixed_formula,
                groups="subject_id",
                re_formula="1",
                vc_formula=VC_FORMULA if full else None,
                data=data,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = model.fit(reml=False, method=method, maxiter=maxiter, disp=False)
            if np.isfinite(fit.llf):
                fits.append((fit, method))
            if fit.converged:
                break
        except Exception as exc:
            failures.append(f"{method}:{type(exc).__name__}")
    if not fits:
        raise RuntimeError("all optimizers failed: " + ",".join(failures))
    converged = [item for item in fits if item[0].converged]
    pool = converged if converged else fits
    return max(pool, key=lambda item: item[0].llf)


def fit_one(values: pd.Series, manifest: pd.DataFrame) -> tuple[dict, dict | None]:
    data = protein_data(values, manifest)
    result = {
        "n_obs": len(data),
        "n_subjects": data["subject_id"].nunique(),
        "converged_null": False,
        "converged_full": False,
        "boundary_full": False,
        "optimizer_null": "",
        "optimizer_full": "",
        "fit_error": "",
    }
    if len(data) < 100 or data["subject_id"].nunique() < 20:
        result["fit_error"] = "insufficient_observations"
        return result, None
    try:
        null, null_method = fit_formula(data, full=False, maxiter=600)
        full, full_method = fit_formula(data, full=True, maxiter=1200)
        result["converged_null"] = bool(null.converged)
        result["converged_full"] = bool(full.converged)
        result["optimizer_null"] = null_method
        result["optimizer_full"] = full_method
        vcomp = np.asarray(full.vcomp, dtype=float)
        result["boundary_full"] = bool(np.any(vcomp < 1e-6))
        lrt = max(0.0, 2 * (full.llf - null.llf))
        beta = full.fe_params
        cov = full.cov_params().loc[beta.index, beta.index]
        rhythm_beta = beta.loc[["cos24", "sin24"]].to_numpy()
        rhythm_cov = cov.loc[["cos24", "sin24"], ["cos24", "sin24"]].to_numpy()
        wald = float(rhythm_beta @ np.linalg.pinv(rhythm_cov) @ rhythm_beta)
        amplitude = float(np.hypot(beta["cos24"], beta["sin24"]))
        peak_elapsed = float(
            (math.atan2(beta["sin24"], beta["cos24"]) * 24 / (2 * np.pi)) % 24
        )
        result.update(
            {
                "ll_null": float(null.llf),
                "ll_full": float(full.llf),
                "lrt_stat": lrt,
                "lrt_chi2_df2_p": float(chi2.sf(lrt, 2)),
                "population_wald_stat": wald,
                "population_rhythm_p": float(chi2.sf(wald, 2)),
                "population_amplitude_log2": amplitude,
                "population_peak_clock": float((peak_elapsed + 9) % 24),
                "random_intercept_variance": float(full.cov_re.iloc[0, 0]),
                "random_cos_variance": float(vcomp[0]),
                "random_sin_variance": float(vcomp[1]),
                "residual_variance": float(full.scale),
            }
        )
        return result, {"data": data, "null": null}
    except Exception as exc:
        result["fit_error"] = f"{type(exc).__name__}: {exc}"
        return result, None


def bootstrap(cache: dict, observed: float, reps: int, rng: np.random.Generator) -> dict:
    data = cache["data"].copy()
    null = cache["null"]
    fixed_mean = np.asarray(null.model.exog) @ np.asarray(null.fe_params)
    tau = math.sqrt(max(float(null.cov_re.iloc[0, 0]), 0))
    sigma = math.sqrt(max(float(null.scale), 0))
    subjects = data["subject_id"].to_numpy()
    unique = pd.unique(subjects)
    simulated = []
    failed = 0
    for _ in range(reps):
        offsets = dict(zip(unique, rng.normal(0, tau, len(unique))))
        data["y"] = fixed_mean + np.array([offsets[x] for x in subjects]) + rng.normal(
            0, sigma, len(data)
        )
        try:
            null_fit, _ = fit_formula(data, full=False, maxiter=400)
            full_fit, _ = fit_formula(data, full=True, maxiter=700)
            if null_fit.converged and full_fit.converged:
                simulated.append(max(0.0, 2 * (full_fit.llf - null_fit.llf)))
            else:
                failed += 1
        except Exception:
            failed += 1
    exceed = sum(value >= observed for value in simulated)
    return {
        "bootstrap_p": (exceed + 1) / (len(simulated) + 1),
        "bootstrap_success": len(simulated),
        "bootstrap_failed": failed,
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(args.input, sep="\t")
    annotations = raw.columns[:4].tolist()
    manifest = sample_manifest(raw.columns[4:].tolist(), raw)
    qc = paper_protein_qc(raw, manifest, annotations)
    selected = qc.index[qc["passes_paper_protein_qc"]]
    rows = []
    caches = {}
    for number, index in enumerate(selected, start=1):
        result, cache = fit_one(raw.loc[index, manifest["sample_column"]], manifest)
        result.update({name: raw.loc[index, name] for name in annotations})
        result["row_index"] = int(index)
        for name in (
            "missing_fraction_primary",
            "empty_detection_rate",
            "empty_to_primary_median_ratio",
            "carryover_flag",
        ):
            result[name] = qc.loc[index, name]
        rows.append(result)
        if cache is not None:
            caches[int(index)] = cache
        if number % 50 == 0:
            print(f"fitted {number}/{len(selected)} proteins", flush=True)
    results = pd.DataFrame(rows)
    valid = (
        results["converged_null"]
        & results["converged_full"]
        & results["fit_error"].fillna("").eq("")
    )
    for p_column, q_column in (
        ("population_rhythm_p", "population_rhythm_q"),
        ("lrt_chi2_df2_p", "heterogeneity_chi2_q"),
    ):
        results[q_column] = np.nan
        use = valid & results[p_column].notna()
        results.loc[use, q_column] = multipletests(results.loc[use, p_column], method="fdr_bh")[1]
    for column in ("bootstrap_p", "bootstrap_success", "bootstrap_failed"):
        results[column] = np.nan
    candidates = results.loc[valid & results["lrt_stat"].notna()].nsmallest(
        args.bootstrap_top, "lrt_chi2_df2_p"
    )
    rng = np.random.default_rng(args.seed)
    for index in candidates.index:
        row_index = int(results.loc[index, "row_index"])
        boot = bootstrap(caches[row_index], float(results.loc[index, "lrt_stat"]), args.bootstrap_reps, rng)
        for name, value in boot.items():
            results.loc[index, name] = value
        print(f"bootstrapped {results.loc[index, 'Genes']}: {boot}", flush=True)
    results.to_csv(args.output_dir / "simplified_mixed_cosinor_results.tsv", sep="\t", index=False)
    manifest.to_csv(args.output_dir / "sample_manifest.tsv", sep="\t", index=False)
    qc.to_csv(args.output_dir / "protein_qc.tsv", sep="\t", index=False)
    summary = {
        "model": "independent subject random intercept, cosine, and sine variances",
        "random_effect_parameters": 3,
        "primary_sample_columns": int(manifest["include_primary"].sum()),
        "paper_qc_proteins": int(qc["passes_paper_protein_qc"].sum()),
        "models_attempted": len(results),
        "both_models_converged": int(valid.sum()),
        "full_models_on_boundary": int((valid & results["boundary_full"]).sum()),
        "population_rhythmic_fdr_0_05": int((results["population_rhythm_q"] < 0.05).sum()),
        "heterogeneity_chi2_fdr_0_05_diagnostic": int(
            (results["heterogeneity_chi2_q"] < 0.05).sum()
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
