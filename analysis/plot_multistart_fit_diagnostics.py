"""Plot multi-start likelihood and variance diagnostics for candidate proteins."""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.regression.mixed_linear_model import MixedLM, MixedLMParams

sys.path.insert(0, str(Path(__file__).resolve().parent))
from first_pass_mixed_cosinor import sample_manifest
from simplified_mixed_cosinor import VC_FORMULA, protein_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit each candidate's null/full mixed-cosinor models from multiple "
            "randomized starting values and plot the results."
        )
    )
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--validation-results", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--genes",
        nargs="+",
        help=(
            "Protein names to analyze. By default, use all rows marked "
            "passes_internal_robustness in the validation results."
        ),
    )
    parser.add_argument("--starts", type=int, default=20)
    parser.add_argument("--maxiter", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=20261015)
    return parser.parse_args()


def select_candidates(validation: pd.DataFrame, genes: list[str] | None) -> pd.DataFrame:
    required = {"gene", "row_index"}
    missing = required.difference(validation.columns)
    if missing:
        raise ValueError(f"validation results are missing columns: {sorted(missing)}")

    if genes is None:
        if "passes_internal_robustness" not in validation.columns:
            raise ValueError(
                "validation results need passes_internal_robustness when --genes is omitted"
            )
        selected = validation.loc[validation["passes_internal_robustness"].astype(bool)]
    else:
        selected = validation.loc[validation["gene"].astype(str).isin(genes)]
        missing_genes = sorted(set(genes).difference(selected["gene"].astype(str)))
        if missing_genes:
            raise ValueError(f"genes not found in validation results: {missing_genes}")
        gene_order = {gene: order for order, gene in enumerate(genes)}
        selected = (
            selected.assign(_gene_order=selected["gene"].astype(str).map(gene_order))
            .sort_values("_gene_order")
            .drop(columns="_gene_order")
        )

    if selected.empty:
        raise ValueError("no candidate proteins selected")
    if selected["gene"].astype(str).duplicated().any():
        raise ValueError("validation results contain duplicate gene entries")
    return selected.reset_index(drop=True)


def starting_params(model: MixedLM, rng: np.random.Generator) -> MixedLMParams:
    intercept_variance = float(np.exp(rng.uniform(np.log(1e-4), np.log(0.1))))
    cov_re = np.array([[intercept_variance]])
    vcomp = (
        np.exp(rng.uniform(np.log(1e-6), np.log(0.03), size=model.k_vc))
        if model.k_vc
        else None
    )
    return MixedLMParams.from_components(
        fe_params=np.zeros(model.k_fe),
        cov_re=cov_re,
        vcomp=vcomp,
    )


def fit_multistart(
    data: pd.DataFrame,
    gene: str,
    starts: int,
    maxiter: int,
    rng: np.random.Generator,
) -> list[dict]:
    records = []
    specifications = (
        ("H0", False, 600),
        ("H1", True, maxiter),
    )
    for model_name, full, iteration_limit in specifications:
        for start in range(starts):
            model = MixedLM.from_formula(
                "y ~ cos24 + sin24 + C(plate)",
                groups="subject_id",
                re_formula="1",
                vc_formula=VC_FORMULA if full else None,
                data=data,
            )
            record = {
                "gene": gene,
                "model": model_name,
                "start": start + 1,
                "method": "lbfgs",
                "converged": False,
                "log_likelihood": np.nan,
                "random_cos_variance": np.nan,
                "random_sin_variance": np.nan,
                "fit_error": "",
            }
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit = model.fit(
                        reml=False,
                        method="lbfgs",
                        maxiter=iteration_limit,
                        start_params=starting_params(model, rng),
                        disp=False,
                    )
                record["converged"] = bool(fit.converged)
                record["log_likelihood"] = float(fit.llf)
                if full:
                    record["random_cos_variance"] = float(fit.vcomp[0])
                    record["random_sin_variance"] = float(fit.vcomp[1])
            except (ValueError, np.linalg.LinAlgError, RuntimeError, FloatingPointError) as exc:
                record["fit_error"] = f"{type(exc).__name__}: {exc}"
            records.append(record)
    return records


def draw_diagnostics(results: pd.DataFrame, gene: str, output: Path) -> None:
    plt.rcParams.update({"font.family": "Noto Sans CJK JP", "axes.unicode_minus": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax_h0, ax_h1, ax_variance = axes[0, 0], axes[0, 1], axes[1, 0]
    ax_h1.sharey(ax_h0)
    axes[1, 1].axis("off")
    fig.suptitle(f"{gene}：多重初期値による混合コサイナーモデル診断", fontsize=15)

    panel_notes = {
        "H0": "被験者ごとのリズム差を含まないモデル",
        "H1": "被験者ごとのcos・sin変動を加えたモデル",
    }
    for model_name, ax in (("H0", ax_h0), ("H1", ax_h1)):
        subset = results.loc[results["model"].eq(model_name)]
        points = subset.loc[subset["log_likelihood"].notna()].copy()
        if not points.empty:
            converged_points = points.loc[points["converged"]]
            reference_points = converged_points if not converged_points.empty else points
            points["delta_log_likelihood"] = (
                points["log_likelihood"] - reference_points["log_likelihood"].max()
            )
            ax.scatter(
                points["start"],
                points["delta_log_likelihood"],
                color="tab:blue",
                marker="o",
            )
            ax.set_xticks(sorted(subset["start"].unique()))
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.5)
        model_title = (
            "H0：個人別リズム差なし"
            if model_name == "H0"
            else "H1：個人別リズム差あり"
        )
        ax.set_title(f"{model_title}\n{panel_notes[model_name]}", fontsize=10)
        ax.set_xlabel("初期値番号")
        ax.set_ylabel("最良値からの Δ 対数尤度")

    full = results.loc[
        results["model"].eq("H1") & results["converged"]
    ]
    for column, color in (
        ("random_cos_variance", "tab:blue"),
        ("random_sin_variance", "tab:orange"),
    ):
        points = full.loc[full[column].notna()]
        if not points.empty:
            ax_variance.scatter(
                points["start"],
                points[column],
                color=color,
            )
    ax_variance.set_xticks(sorted(full["start"].unique()))
    ax_variance.set_title("H1：被験者別リズム変動の分散", fontsize=10)
    ax_variance.set_xlabel("初期値番号")
    ax_variance.set_ylabel("分散（log2強度の二乗）")

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.starts < 1:
        raise ValueError("--starts must be at least 1")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matrix = pd.read_csv(args.matrix, sep="\t")
    validation = pd.read_csv(args.validation_results, sep="\t")
    candidates = select_candidates(validation, args.genes)
    manifest = sample_manifest(matrix.columns[4:].tolist(), matrix)
    rng = np.random.default_rng(args.seed)
    all_results = []

    for _, candidate in candidates.iterrows():
        gene = str(candidate["gene"])
        row_index = int(candidate["row_index"])
        if row_index not in matrix.index:
            raise ValueError(f"row_index {row_index} for {gene} is not present in the matrix")
        values = matrix.loc[row_index, manifest["sample_column"]]
        data = protein_data(values, manifest)
        if len(data) < 100 or data["subject_id"].nunique() < 20:
            raise ValueError(
                f"{gene} has insufficient data for fitting "
                f"({len(data)} observations, {data['subject_id'].nunique()} subjects)"
            )
        gene_results = pd.DataFrame(
                fit_multistart(data, gene, args.starts, args.maxiter, rng)
        )
        all_results.append(gene_results)
        draw_diagnostics(
            gene_results,
            gene,
            args.output_dir / f"{gene}_multistart_fit_diagnostics.png",
        )
        converged = int(gene_results["converged"].sum())
        print(f"{gene}: {converged}/{len(gene_results)} fits converged", flush=True)

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(args.output_dir / "multistart_fit_results.tsv", sep="\t", index=False)
    if combined["fit_error"].ne("").any():
        print(
            f"Warning: {int(combined['fit_error'].ne('').sum())} fits failed; "
            "see fit_error in multistart_fit_results.tsv",
            flush=True,
        )


if __name__ == "__main__":
    main()
