#!/usr/bin/env python3
"""Plot one observed 24-hour trajectory per participant for candidate proteins."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager


CANDIDATES = ["DCD", "APOC2", "ADIPOQ", "APOA4", "LMAN2"]
PROTEIN_COLORS = {
    "DCD": "#E64B35",
    "APOC2": "#4DBBD5",
    "ADIPOQ": "#00A087",
    "APOA4": "#F39B7F",
    "LMAN2": "#3C5488",
}
TIME_GRID = np.arange(0, 25, 3)
TIME_LABELS = ["09", "12", "15", "18", "21", "00", "03", "06", "09"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def configure_fonts() -> None:
    font_path = "/usr/share/fonts/google-noto-cjk/NotoSansCJK-Regular.ttc"
    if Path(font_path).exists():
        font_manager.fontManager.addfont(font_path)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=font_path).get_name()
    plt.rcParams["axes.unicode_minus"] = False


def long_candidate_data(matrix: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    meta = manifest.loc[manifest["include_primary"].astype(bool)].copy()
    protein_rows = matrix.loc[matrix["Genes"].isin(CANDIDATES)].set_index("Genes")
    missing = set(CANDIDATES) - set(protein_rows.index)
    if missing:
        raise ValueError(f"Candidate proteins absent from matrix: {sorted(missing)}")
    records = []
    for gene in CANDIDATES:
        values = pd.to_numeric(protein_rows.loc[gene, meta["sample_column"]], errors="coerce")
        for (_, sample), value in zip(meta.iterrows(), values.to_numpy()):
            records.append(
                {
                    "gene": gene,
                    "subject_id": str(int(sample["subject_id"])),
                    "plate": int(sample["plate"]),
                    "elapsed_hour": int(sample["elapsed_hour"]),
                    "log2_intensity": np.log2(value) if pd.notna(value) and value > 0 else np.nan,
                }
            )
    data = pd.DataFrame(records)
    data["centered_log2"] = data.groupby(["gene", "subject_id"])["log2_intensity"].transform(
        lambda x: x - x.mean()
    )
    return data


def draw_panel(ax, data: pd.DataFrame, value_column: str, gene: str) -> None:
    vivid = PROTEIN_COLORS[gene]
    trajectories = []
    for subject, group in data.groupby("subject_id"):
        series = group.set_index("elapsed_hour")[value_column].reindex(TIME_GRID)
        ax.plot(
            TIME_GRID,
            series,
            marker="o",
            markersize=3.0,
            linewidth=1.15,
            alpha=0.23,
            color=vivid,
        )
        trajectories.append(series.to_numpy(dtype=float))
    mean_curve = np.nanmean(np.vstack(trajectories), axis=0)
    ax.plot(TIME_GRID, mean_curve, color=vivid, linewidth=3.2, marker="o", markersize=4.0)
    ax.axvspan(14, 23, color="#243447", alpha=0.08, linewidth=0)
    ax.set_xticks(TIME_GRID, TIME_LABELS)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)


def make_combined(data: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(5, 1, figsize=(11, 15), sharex=True)
    for ax, gene in zip(axes, CANDIDATES):
        draw_panel(ax, data.loc[data["gene"].eq(gene)], "centered_log2", gene)
        ax.set_title(gene, loc="left", fontsize=14, fontweight="bold")
        ax.set_ylabel("個人平均からの差\n(log2)")
    axes[-1].set_xlabel("採血時刻（1日目09:00 → 2日目09:00）")
    fig.suptitle("5候補タンパク質の個人別24時間変動", fontsize=18, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.003,
        "淡色の細線：各被験者　鮮色の太線：被験者平均　灰色領域：夜間の睡眠機会　欠測値は補完していない。",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 0.98))
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_cosinor_fits(data: pd.DataFrame, output: Path) -> None:
    smooth = np.linspace(0, 24, 241)
    smooth_design = np.column_stack(
        [np.ones(len(smooth)), np.cos(2 * np.pi * smooth / 24), np.sin(2 * np.pi * smooth / 24)]
    )
    fig, axes = plt.subplots(5, 1, figsize=(11, 15), sharex=True)
    for ax, gene in zip(axes, CANDIDATES):
        vivid = PROTEIN_COLORS[gene]
        subset = data.loc[data["gene"].eq(gene)]
        fitted_curves = []
        for subject, group in subset.groupby("subject_id"):
            observed = group.dropna(subset=["log2_intensity"])
            if len(observed) < 5:
                continue
            angle = 2 * np.pi * observed["elapsed_hour"].to_numpy(dtype=float) / 24
            design = np.column_stack([np.ones(len(observed)), np.cos(angle), np.sin(angle)])
            coef = np.linalg.lstsq(design, observed["log2_intensity"].to_numpy(), rcond=None)[0]
            curve = smooth_design @ coef - coef[0]
            fitted_curves.append(curve)
            ax.plot(smooth, curve, color=vivid, linewidth=1.2, alpha=0.22)
            centered_points = observed["log2_intensity"] - observed["log2_intensity"].mean()
            ax.scatter(
                observed["elapsed_hour"], centered_points, color=vivid, s=10, alpha=0.18
            )
        ax.plot(smooth, np.mean(fitted_curves, axis=0), color=vivid, linewidth=3.2)
        ax.axvspan(14, 23, color="#243447", alpha=0.08, linewidth=0)
        ax.set_title(gene, loc="left", fontsize=14, fontweight="bold")
        ax.set_ylabel("個人平均からの差\n(log2)")
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xticks(TIME_GRID, TIME_LABELS)
    axes[-1].set_xlabel("採血時刻（1日目09:00 → 2日目09:00）")
    fig.suptitle("5候補タンパク質の個人別24時間フィット曲線", fontsize=18, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.003,
        "淡色の曲線：各被験者のcosinorフィット　淡色の点：実測値　鮮色の太線：個人フィットの平均。",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 0.98))
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_individual(data: pd.DataFrame, output_dir: Path) -> None:
    for gene in CANDIDATES:
        subset = data.loc[data["gene"].eq(gene)]
        fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
        draw_panel(axes[0], subset, "log2_intensity", gene)
        axes[0].set_title("実測強度：個人平均と曲線差を同時に表示")
        axes[0].set_ylabel("タンパク質強度（log2）")
        draw_panel(axes[1], subset, "centered_log2", gene)
        axes[1].set_title("個人内中心化：曲線形状の差を強調")
        axes[1].set_ylabel("個人平均からの差（log2）")
        for ax in axes:
            ax.set_xlabel("採血時刻（1日目09:00 → 2日目09:00）")
        fig.suptitle(f"{gene}：淡色の各線は1名の被験者を表す", fontsize=17, fontweight="bold")
        fig.text(
            0.5,
            0.01,
            "鮮色の太線：被験者平均　灰色領域：夜間の睡眠機会　欠測値は補完していない。",
            ha="center",
            fontsize=10,
        )
        fig.tight_layout(rect=(0, 0.04, 1, 0.92))
        fig.savefig(output_dir / f"{gene}_individual_trajectories.png", dpi=220, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_fonts()
    matrix = pd.read_csv(args.matrix, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    data = long_candidate_data(matrix, manifest)
    data.to_csv(args.output_dir / "candidate_trajectory_data.tsv", sep="\t", index=False)
    make_combined(data, args.output_dir / "five_candidates_centered_trajectories.png")
    make_cosinor_fits(data, args.output_dir / "five_candidates_individual_cosinor_fits.png")
    make_individual(data, args.output_dir)
    print(f"wrote plots for {', '.join(CANDIDATES)} to {args.output_dir}")


if __name__ == "__main__":
    main()
