"""Export observed individual and mean trajectories for the BH-screened candidates."""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/rhythm-candidate-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/ten_candidate_curves"
GRID = np.arange(0, 25, 3)
LABELS = ["09", "12", "15", "18", "21", "00", "03", "06", "09*"]
BLUE = "#176A99"
GRAY = "#8295AA"


def load_data():
    screen = pd.read_csv(
        ROOT / "outputs/simplified_mixed_cosinor/simplified_mixed_cosinor_results.tsv",
        sep="\t",
    )
    candidates = screen.loc[screen.heterogeneity_chi2_q.lt(0.05)].sort_values(
        "lrt_chi2_df2_p"
    )
    assert len(candidates) == 10
    raw = pd.read_csv(ROOT / "data/report.pg_matrix.tsv", sep="\t")
    manifest = pd.read_csv(
        ROOT / "outputs/simplified_mixed_cosinor/sample_manifest.tsv", sep="\t"
    )
    meta = manifest.loc[manifest.include_primary.eq(True)].copy()
    assert len(meta) == 206 and meta.subject_id.nunique() == 24
    records = []
    for _, candidate in candidates.iterrows():
        row = raw.iloc[int(candidate.row_index)]
        assert row.Genes == candidate.Genes
        values = pd.to_numeric(row[meta.sample_column], errors="coerce").to_numpy(float)
        part = meta[["subject_id", "plate", "elapsed_hour", "sample_column"]].copy()
        part["gene"] = candidate.Genes
        part["log2_intensity"] = np.log2(np.where(values > 0, values, np.nan))
        records.append(part)
    data = pd.concat(records, ignore_index=True)
    data["subject_id"] = data.subject_id.astype(int).astype(str)
    data["centered_log2"] = data.log2_intensity - data.groupby(
        ["gene", "subject_id"]
    ).log2_intensity.transform("mean")
    assert not data.duplicated(["gene", "subject_id", "elapsed_hour"]).any()
    assert np.allclose(
        data.groupby(["gene", "subject_id"]).centered_log2.mean().dropna(), 0, atol=1e-12
    )
    return candidates, data


def draw(ax, subset, column, candidate, compact=False):
    curves = subset.pivot(index="elapsed_hour", columns="subject_id", values=column)
    curves = curves.reindex(GRID)
    for subject in curves.columns:
        ax.plot(GRID, curves[subject], color=GRAY, alpha=0.48, lw=0.85,
                marker="o", ms=2, zorder=1)
    mean = curves.mean(axis=1)
    ax.plot(GRID, mean, color=BLUE, lw=2.8, marker="o", ms=3.5, zorder=3)
    if column == "centered_log2":
        ax.axhline(0, color="#425466", lw=0.7, alpha=0.45, zorder=0)
        limit = float(np.nanmax(np.abs(curves.to_numpy()))) * 1.12
        ax.set_ylim(-max(limit, 0.05), max(limit, 0.05))
    ax.set_title(candidate.Genes, loc="left", fontsize=14 if compact else 18,
                 fontweight="bold", color="#16324A", pad=18)
    ax.text(0, 1.025,
            f"异质性初筛 q = {candidate.heterogeneity_chi2_q:.3g}",
            transform=ax.transAxes, fontsize=8 if compact else 10, color="#53677B")
    ax.set_xticks(GRID, LABELS, fontsize=8 if compact else 10)
    ax.tick_params(axis="y", labelsize=8 if compact else 10, length=0)
    ax.set_xlim(-0.5, 24.5)
    ax.set_xlabel("采血时刻（末次为次日09时）", fontsize=8 if compact else 10)
    ax.set_ylabel("相对个人均值（log2）" if column == "centered_log2"
                  else "观测强度（log2）", fontsize=9 if compact else 11)
    ax.grid(axis="y", color="#E2E8EE", lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#C8D2DD")


def decorate(fig, title, subtitle):
    fig.suptitle(title, x=0.055, ha="left", y=0.97, fontsize=23,
                 fontweight="bold", color="#16324A")
    fig.text(0.055, 0.915, subtitle, fontsize=12, color="#53677B")
    fig.legend(handles=[
        Line2D([0], [0], color=GRAY, lw=1.2, marker="o", ms=3, label="细灰线：每位受试者的观测轨迹"),
        Line2D([0], [0], color=BLUE, lw=3, marker="o", ms=4, label="粗蓝线：该时刻有观测者的算术平均"),
    ], loc="lower left", bbox_to_anchor=(0.05, 0.035), frameon=False,
        ncol=2, fontsize=11)
    fig.text(0.055, 0.022,
             "缺失值不填补，个人曲线在缺失处断开；各面板纵轴范围不同。观测均值未经plate调整，不是混合模型固定效应曲线。",
             fontsize=9, color="#53677B")


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=240, facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", facecolor="white")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Noto Sans CJK JP", "axes.unicode_minus": False,
                         "pdf.fonttype": 42})
    candidates, data = load_data()
    data.to_csv(OUT / "observed_trajectories.tsv", sep="\t", index=False)
    summary = data.groupby(["gene", "elapsed_hour"]).agg(
        n_observed=("log2_intensity", "count"),
        mean_log2=("log2_intensity", "mean"),
        mean_centered_log2=("centered_log2", "mean"),
    ).reset_index()
    summary.to_csv(OUT / "timepoint_means.tsv", sep="\t", index=False)
    candidates[["Genes", "lrt_stat", "lrt_chi2_df2_p", "heterogeneity_chi2_q",
                "population_rhythm_q"]].to_csv(OUT / "candidate_screening.tsv", sep="\t", index=False)
    for column, label in [("centered_log2", "centered"), ("log2_intensity", "raw_log2")]:
        fig, axes = plt.subplots(2, 5, figsize=(20, 11.25))
        for ax, (_, candidate) in zip(axes.flat, candidates.iterrows()):
            draw(ax, data.loc[data.gene.eq(candidate.Genes)], column, candidate, compact=True)
        decorate(fig, "10个初筛候选：群体平均与个人24小时轨迹",
                 "个人内中心化：减去每位受试者自身的观测均值，展示动态差异" if label == "centered"
                 else "原始log2强度：同时展示个人平均水平与时间变化")
        fig.subplots_adjust(left=0.055, right=0.985, bottom=0.15, top=0.82,
                            wspace=0.39, hspace=0.63)
        save(fig, f"ten_candidates_{label}")
    for page in range(2):
        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        for ax, (_, candidate) in zip(axes.flat, candidates.iloc[page*5:page*5+5].iterrows()):
            draw(ax, data.loc[data.gene.eq(candidate.Genes)], "centered_log2", candidate)
        axes.flat[-1].axis("off")
        axes.flat[-1].text(0.04, 0.85,
            "如何读图\n\n每条细线对应一位受试者\n粗线为各时刻观测均值\n\n个人内中心化去除平均水平差异\n曲线用于描述，不替代LRT检验\n\n候选按初筛LRT p值排序",
            va="top", fontsize=13, color="#53677B", linespacing=1.6)
        decorate(fig, f"10个初筛候选的个人动态（{page+1}/2）",
                 "24名受试者 · 206个质控后样本列 · 24小时观测 · log2强度个人内中心化")
        fig.subplots_adjust(left=0.065, right=0.98, bottom=0.18, top=0.79,
                            wspace=0.3, hspace=0.91)
        save(fig, f"ppt_candidates_page_{page+1}")
    for _, candidate in candidates.iterrows():
        fig, axes = plt.subplots(1, 2, figsize=(16, 9))
        for ax, column in zip(axes, ["log2_intensity", "centered_log2"]):
            draw(ax, data.loc[data.gene.eq(candidate.Genes)], column, candidate)
        decorate(fig, f"{candidate.Genes}：群体平均与个人24小时轨迹",
                 "左：原始log2强度，保留个人平均水平差异    |    右：个人内中心化，突出时间变化")
        fig.subplots_adjust(left=0.075, right=0.98, bottom=0.18, top=0.77, wspace=0.23)
        save(fig, f"{candidate.Genes}_mean_and_individual")
    (OUT / "README.md").write_text(
        "# 10个候选的PPT曲线图\n\n"
        "建议正文使用 `ppt_candidates_page_1.png` 和 `ppt_candidates_page_2.png`；"
        "总览及每个蛋白的双面板图同时提供PNG与PDF。\n\n"
        "候选条件：简化模型Null与Full均收敛、无fit_error；将有效模型的LRT近似卡方(df=2)"
        "p值在402个蛋白内作BH校正，heterogeneity_chi2_q < 0.05，共10个。"
        "随机方差边界未在初筛时排除，因此包括IGHG3。200次bootstrap用于后续内部验证。\n\n"
        "细灰线是每位受试者的观测点连线，不是cosinor拟合。粗蓝线是在每个时点对"
        "可用受试者的log2观测值求算术平均；中心化版本先减去各人全部可用时点的均值。"
        "个人缺测点处断线，无填补；各时点实际人数见timepoint_means.tsv。"
        "末次09时为次日，未与首次09时合并。\n\n"
        "图中曲线不进行plate校正，不代表混合模型固定效应或个人随机效应预测；"
        "不显示置信区间，不据此单独推断节律显著或振幅/相位机制。各面板纵轴范围不同。"
        "完整缺失观测也保留在observed_trajectories.tsv中。\n\n"
        "复现：`/home/pan/miniconda3/envs/rhythm-agent/bin/python analysis/plot_ten_candidate_curves.py`\n",
        encoding="utf-8",
    )
    print(f"Saved 14 figures (PNG + PDF) and source tables to {OUT}")


if __name__ == "__main__":
    main()
