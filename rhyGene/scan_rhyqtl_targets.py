#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Scan Chen et al. rhyQTL SupplementaryData2 txt files for six target genes.

Targets:
    APOA4, LMAN2, APOC4, IGHV2-5, DCD, APOC2

Outputs:
    1) rhyQTL_target_hits.csv
       All matching rows across all txt files.

    2) rhyQTL_target_summary.csv
       One summary row per Gene × source file/tissue:
       Gene, Tissue/File, # rows, Lead SNP, min pval,
       Chosen.model, G.test_pval, phase/amp values and differences.

Usage:
    Put this script in the extracted SupplementaryData2 folder, then run:

        python scan_rhyqtl_targets.py

    Or specify another folder:

        python scan_rhyqtl_targets.py /path/to/SupplementaryData2

Requirements:
    pip install pandas
"""

import sys
from pathlib import Path
import pandas as pd

TARGET_GENES = ["APOA4", "LMAN2", "APOC4", "IGHV2-5", "DCD", "APOC2"]

# Columns we care about if they exist
KEEP_COLUMNS = [
    "ID",
    "Chromosome",
    "Position",
    "REF",
    "ALT",
    "rhyGene.ID",
    "rhyGene.name",
    "Sample.size.0",
    "Sample.size.1",
    "Sample.size.2",
    "pval_0",
    "phase_0",
    "amp_0",
    "pval_1",
    "phase_1",
    "amp_1",
    "Chosen.model",
    "G.test_pval",
    "pval",
]


def read_table(path: Path) -> pd.DataFrame:
    """Read a tab-delimited rhyQTL result file robustly."""
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        # Fallback: let pandas infer separator
        return pd.read_csv(path, sep=None, engine="python", low_memory=False)


def to_numeric_if_present(df: pd.DataFrame, cols):
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def derive_tissue_name(path: Path) -> str:
    """
    Use the filename stem as tissue/source name.
    This is intentionally conservative because SupplementaryData2 filenames
    may encode tissue names differently.
    """
    return path.stem


def main():
    root = Path(sys.argv[1]).expanduser().resolve() if len(sys.argv) > 1 else Path(".").resolve()

    if not root.exists():
        raise SystemExit(f"Folder does not exist: {root}")

    txt_files = sorted(root.rglob("*.txt"))

    if not txt_files:
        raise SystemExit(
            f"No .txt files found under:\n{root}\n\n"
            "Put this script inside the extracted SupplementaryData2 folder "
            "or pass the folder path as the first argument."
        )

    print(f"Scanning folder: {root}")
    print(f"Found {len(txt_files)} txt files.")
    print("Target genes:", ", ".join(TARGET_GENES))
    print()

    hits = []
    skipped = []

    for i, file in enumerate(txt_files, start=1):
        try:
            df = read_table(file)
        except Exception as e:
            skipped.append((str(file), f"read error: {e}"))
            continue

        if "rhyGene.name" not in df.columns:
            skipped.append((str(file), "missing column: rhyGene.name"))
            continue

        # Exact symbol match only
        gene_col = df["rhyGene.name"].astype(str).str.strip()
        hit = df[gene_col.isin(TARGET_GENES)].copy()

        if hit.empty:
            continue

        hit.insert(0, "source_file", file.name)
        hit.insert(1, "source_path", str(file.relative_to(root)))
        hit.insert(2, "tissue_or_source", derive_tissue_name(file))

        # Keep original columns, but make sure target columns are near the front
        front = ["source_file", "source_path", "tissue_or_source"]
        cols = front + [c for c in KEEP_COLUMNS if c in hit.columns]
        other_cols = [c for c in hit.columns if c not in cols]
        hit = hit[cols + other_cols]

        hits.append(hit)

        genes_here = ", ".join(sorted(hit["rhyGene.name"].astype(str).unique()))
        print(f"[{i:>3}/{len(txt_files)}] HIT  {file.name}: {genes_here} ({len(hit)} rows)")

    if not hits:
        print("\nNo matching rows found for the six target genes.")
        if skipped:
            print(f"{len(skipped)} files were skipped.")
        return

    all_hits = pd.concat(hits, ignore_index=True)

    # Numeric conversion for ranking / differences
    numeric_cols = [
        "pval",
        "G.test_pval",
        "pval_0",
        "pval_1",
        "phase_0",
        "phase_1",
        "amp_0",
        "amp_1",
    ]
    all_hits = to_numeric_if_present(all_hits, numeric_cols)

    # Calculate interpretable differences when both columns exist
    if "phase_0" in all_hits.columns and "phase_1" in all_hits.columns:
        all_hits["phase_diff_1_minus_0"] = all_hits["phase_1"] - all_hits["phase_0"]

    if "amp_0" in all_hits.columns and "amp_1" in all_hits.columns:
        all_hits["amp_diff_1_minus_0"] = all_hits["amp_1"] - all_hits["amp_0"]

    # Save every matching row
    detailed_out = root / "rhyQTL_target_hits.csv"
    all_hits.to_csv(detailed_out, index=False)

    # Build one row per Gene × file/tissue.
    summary_rows = []

    group_cols = ["rhyGene.name", "source_file", "source_path", "tissue_or_source"]

    for keys, g in all_hits.groupby(group_cols, dropna=False, sort=True):
        gene, source_file, source_path, tissue = keys
        g = g.copy()

        # Lead row = smallest final pval if available; otherwise smallest G.test_pval
        ranking_col = None
        if "pval" in g.columns and g["pval"].notna().any():
            ranking_col = "pval"
        elif "G.test_pval" in g.columns and g["G.test_pval"].notna().any():
            ranking_col = "G.test_pval"

        if ranking_col is not None:
            lead_idx = g[ranking_col].idxmin()
            lead = g.loc[lead_idx]
        else:
            lead = g.iloc[0]

        row = {
            "Gene": gene,
            "Tissue_or_source": tissue,
            "Source_file": source_file,
            "Source_path": source_path,
            "Num_rows": len(g),
            "Lead_SNP": lead.get("ID", pd.NA),
            "Chromosome": lead.get("Chromosome", pd.NA),
            "Position": lead.get("Position", pd.NA),
            "REF": lead.get("REF", pd.NA),
            "ALT": lead.get("ALT", pd.NA),
            "min_pval": g["pval"].min() if "pval" in g.columns else pd.NA,
            "Lead_Chosen.model": lead.get("Chosen.model", pd.NA),
            "Lead_G.test_pval": lead.get("G.test_pval", pd.NA),
            "Lead_pval_0": lead.get("pval_0", pd.NA),
            "Lead_phase_0": lead.get("phase_0", pd.NA),
            "Lead_amp_0": lead.get("amp_0", pd.NA),
            "Lead_pval_1": lead.get("pval_1", pd.NA),
            "Lead_phase_1": lead.get("phase_1", pd.NA),
            "Lead_amp_1": lead.get("amp_1", pd.NA),
            "Lead_phase_diff_1_minus_0": lead.get("phase_diff_1_minus_0", pd.NA),
            "Lead_amp_diff_1_minus_0": lead.get("amp_diff_1_minus_0", pd.NA),
        }
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)

    # Sort genes in the requested order, then by min p-value
    gene_order = {g: i for i, g in enumerate(TARGET_GENES)}
    summary["_gene_order"] = summary["Gene"].map(gene_order)
    summary = summary.sort_values(
        ["_gene_order", "min_pval", "Tissue_or_source"],
        na_position="last"
    ).drop(columns="_gene_order")

    summary_out = root / "rhyQTL_target_summary.csv"
    summary.to_csv(summary_out, index=False)

    # Gene-level overview
    print("\n" + "=" * 100)
    print("GENE-LEVEL OVERVIEW")
    print("=" * 100)

    overview_rows = []
    for gene in TARGET_GENES:
        g = all_hits[all_hits["rhyGene.name"].astype(str) == gene]
        if g.empty:
            overview_rows.append({
                "Gene": gene,
                "Found": "NO",
                "Files/Tissues": 0,
                "Total_rows": 0,
                "Best_pval": pd.NA,
            })
        else:
            overview_rows.append({
                "Gene": gene,
                "Found": "YES",
                "Files/Tissues": g["source_path"].nunique(),
                "Total_rows": len(g),
                "Best_pval": g["pval"].min() if "pval" in g.columns else pd.NA,
            })

    overview = pd.DataFrame(overview_rows)
    print(overview.to_string(index=False))

    # Main requested summary
    display_cols = [
        "Gene",
        "Tissue_or_source",
        "Num_rows",
        "Lead_SNP",
        "min_pval",
        "Lead_Chosen.model",
        "Lead_G.test_pval",
        "Lead_phase_0",
        "Lead_amp_0",
        "Lead_phase_1",
        "Lead_amp_1",
        "Lead_phase_diff_1_minus_0",
        "Lead_amp_diff_1_minus_0",
    ]
    display_cols = [c for c in display_cols if c in summary.columns]

    print("\n" + "=" * 100)
    print("GENE × TISSUE/FILE SUMMARY")
    print("=" * 100)

    with pd.option_context(
        "display.max_rows", None,
        "display.max_columns", None,
        "display.width", 220,
        "display.max_colwidth", 40,
    ):
        print(summary[display_cols].to_string(index=False))

    print("\nSaved:")
    print(f"  Detailed hits : {detailed_out}")
    print(f"  Summary       : {summary_out}")

    if skipped:
        print(f"\nNote: {len(skipped)} files were skipped.")
        skipped_out = root / "rhyQTL_skipped_files.csv"
        pd.DataFrame(skipped, columns=["file", "reason"]).to_csv(skipped_out, index=False)
        print(f"  Skipped log   : {skipped_out}")


if __name__ == "__main__":
    main()
