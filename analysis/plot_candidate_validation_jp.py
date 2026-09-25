"""Create a Japanese slide summarizing the existing candidate validation results."""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/rhythm-candidate-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/simplified_candidate_validation"
ORDER = ["APOC4", "LMAN2", "APOC2", "APOA4", "DCD", "IGHV2-5",
         "IGLV1-51", "ADIPOQ", "IGLC3", "IGHG3"]
PASS = "#E1F1EB"
FAIL = "#F7E3E1"
RISK = "#FFF0D6"
INK = "#17344B"


def main():
    data = pd.read_csv(OUT / "candidate_validation_summary.tsv", sep="\t").set_index("gene")
    assert len(data) == 10
    assert data.passes_internal_robustness.sum() == 6
    assert data.high_priority_clean.sum() == 2
    plt.rcParams.update({"font.family": "Noto Sans CJK JP", "axes.unicode_minus": False,
                         "pdf.fonttype": 42})
    fig = plt.figure(figsize=(16, 9), facecolor="white")
    fig.text(0.045, 0.945, "内部検証と技術的リスク評価による候補の絞り込み",
             fontsize=23, weight="bold", color=INK)
    fig.text(0.045, 0.895, "10候補の曲線を確認した後、個人依存性・モデル感度・ブランク信号を評価",
             fontsize=12, color="#526779")
    cards = [
        (0.045, "初期スクリーニング", "10候補", "近似LRTのBH補正 q < 0.05", "#EDF3F8"),
        (0.365, "内部頑健性基準を満たす", "6候補", "4候補は全基準を満たさず", "#E8F1F7"),
        (0.685, "優先検証候補", "APOC4 / LMAN2", "残る4候補は技術的リスクを付記", PASS),
    ]
    for x, title, result, note, color in cards:
        ax = fig.add_axes((x, 0.70, 0.27, 0.16), facecolor=color)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.06, 0.76, title, fontsize=12, color=INK)
        ax.text(0.06, 0.39, result, fontsize=22, weight="bold", color=INK)
        ax.text(0.06, 0.12, note, fontsize=10, color="#526779")
    for x in (0.333, 0.653):
        fig.text(x, 0.768, "→", fontsize=25, color="#8297AA", ha="center")
    headers = ["候補", "被験者を1名ずつ除外\n有意 / 全24回", "時間原点を変更\n有意 / 全5条件",
               "線形トレンド\n調整後", "Bootstrap\np値（200回）", "ランダム分散\n境界解", "ブランク\n検出", "総合評価"]
    rows, colors = [], []
    for gene in ORDER:
        row = data.loc[gene]
        loo = round(row.loo_p_lt_0_05_fraction * row.loo_success)
        origins = round(row.origin_p_lt_0_05_fraction * row.origin_success)
        clean = bool(row.high_priority_clean)
        robust = bool(row.passes_internal_robustness)
        status = "優先候補" if clean else "技術的リスクあり" if robust else "基準未達"
        rows.append([gene, f"{loo} / 24", f"{origins} / 5",
                     "p < 0.05" if row.linear_trend_p < 0.05 else "p ≥ 0.05",
                     f"{row.bootstrap_p:.5f}", "あり" if row.boundary_full else "なし",
                     "検出" if row.carryover_flag else "未検出", status])
        colors.append(["#F4F7FA", PASS if loo == 24 else FAIL,
                       PASS if origins == 5 else FAIL,
                       PASS if row.linear_trend_p < 0.05 else FAIL,
                       PASS if row.bootstrap_p < 0.05 else FAIL,
                       FAIL if row.boundary_full else PASS,
                       RISK if row.carryover_flag else PASS,
                       PASS if clean else RISK if robust else FAIL])
    ax = fig.add_axes((0.045, 0.225, 0.91, 0.435))
    ax.axis("off")
    table = ax.table(cellText=rows, cellColours=colors, colLabels=headers,
                     colColours=[INK] * len(headers), cellLoc="center", bbox=[0, 0, 1, 1],
                     colWidths=[0.11, 0.16, 0.16, 0.115, 0.13, 0.10, 0.08, 0.145])
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_linewidth(2)
        cell.get_text().set_color("white" if r == 0 else INK)
        if r == 0:
            cell.set_height(cell.get_height() * 1.4)
            cell.get_text().set_fontsize(10)
        if c == 0:
            cell.get_text().set_weight("bold")
    fig.text(0.045, 0.165, "主な絞り込み要因：被験者除外で4候補が不安定。技術的リスク評価で優先候補を2つに限定。",
             fontsize=13, weight="bold", color=INK)
    fig.text(0.045, 0.113,
             "感度分析は名目上の p < 0.05。全候補でBootstrap 200/200回が成功し、p = 1/201（最小分解能）。",
             fontsize=10, color="#526779")
    fig.text(0.045, 0.082,
             "時間原点：0・3・6・9・12時間。境界解：ランダムcos/sin分散 < 1e-6。ブランク検出だけで偽陽性とは断定しない。",
             fontsize=10, color="#526779")
    fig.text(0.045, 0.05,
             "同一データ内の予備的検証。Bootstrapの全タンパク質に対する多重検定補正・独立集団での再現検証は未実施。",
             fontsize=10, color="#526779")
    for extension in ("png", "pdf"):
        path = OUT / f"candidate_selection_jp.{extension}"
        fig.savefig(path, dpi=240, facecolor="white")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
