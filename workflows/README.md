# Snakemake boundary

This directory owns deterministic, file-based scientific computation for node `5`.

- `Snakefile`: workflow entry point.
- `rules/gwas_regional.smk`: C.b regional GWAS computation.
- `rules/colocalization.smk`: D replication and robustness computation.
- `scripts/`: scripts invoked by rules.
- `schemas/`: JSON validation at the Python/Snakemake boundary.
- `profiles/`: execution profiles.

The flow is:

`C.b.2 or D.4 -> ScientificJobSpec -> node 5 -> Snakemake -> ResultManifest -> C.b.3 or D.5`

Snakemake rule state must not be copied into `RhythmResearchState`.

## Protein-rhythm workflow

The validated first-stage protein-rhythm workflow is configured in
`configs/protein_rhythm.yaml` and exposed as the target `protein_rhythm_all`.

It runs:

`raw matrix -> paper-compatible QC -> simplified mixed cosinor -> candidate
validation -> Japanese trajectory plots -> hashed ScientificResultManifest`

Install the workflow dependencies and inspect the DAG:

```bash
python -m pip install -e '.[workflow]'
snakemake \
  --snakefile workflows/Snakefile \
  --configfile configs/protein_rhythm.yaml \
  --profile workflows/profiles/local \
  --dry-run protein_rhythm_all
```

Run it with:

```bash
snakemake \
  --snakefile workflows/Snakefile \
  --configfile configs/protein_rhythm.yaml \
  --profile workflows/profiles/local \
  protein_rhythm_all
```

The default output root is `outputs/snakemake_protein_rhythm`. The validation
stage is intentionally compute-heavy because it performs leave-one-subject-out,
time-origin and linear-trend sensitivity checks plus 200 parametric bootstrap
replicates per screened candidate.
