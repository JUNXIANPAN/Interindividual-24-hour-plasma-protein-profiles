# Interindividual 24-hour plasma protein profiles

Code for screening 24-hour protein profiles and testing whether subjects differ
in their protein-specific circadian curves.

## Analysis overview

The workflow uses LC-MS protein intensities from repeated samples collected over
the 24-hour cycle. It applies the published sample-level quality-control rule,
log2-transforms observed intensities without imputation, excludes timestamped
technical reinjections, and fits mixed-effects cosinor models with plate
indicators.

Two model specifications are provided:

- `first_pass_mixed_cosinor.py`: correlated subject random intercept, cosine,
	and sine effects.
- `simplified_mixed_cosinor.py`: independent subject baseline, cosine, and sine
	variance components. This specification improves model identifiability for
	the available sample size.

Candidate validation includes parametric bootstrap likelihood-ratio tests,
leave-one-subject-out checks, time-origin sensitivity checks, and adjustment for
a linear time trend.

## Repository layout

```text
analysis/   Analysis, validation, plotting, and manifest scripts
configs/    Project configuration files
docs/       Design notes
src/        Reusable package components
tests/      Unit and scientific tests
workflows/  Snakemake workflow definitions
```

## Main scripts

- `analysis/first_pass_mixed_cosinor.py`
- `analysis/simplified_mixed_cosinor.py`
- `analysis/validate_simplified_candidates.py`
- `analysis/plot_candidate_trajectories.py`
- `analysis/plot_ten_candidate_curves.py`
- `analysis/build_protein_rhythm_manifest.py`

## Installation

Python 3.11 or newer is required. The statistical workflow uses the optional
`workflow` dependencies:

```bash
python -m pip install -e '.[workflow]'
```

## Running the analysis

The input matrix is a tab-separated file whose first four columns contain
protein annotations and whose remaining columns contain sample intensities.
The matrix is not distributed in this repository.

```bash
python analysis/first_pass_mixed_cosinor.py \
	--input data/report.pg_matrix.tsv \
	--output-dir outputs/first_pass_mixed_cosinor \
	--bootstrap-top 10 \
	--bootstrap-reps 50

python analysis/simplified_mixed_cosinor.py \
	--input data/report.pg_matrix.tsv \
	--output-dir outputs/simplified_mixed_cosinor \
	--bootstrap-top 10 \
	--bootstrap-reps 50
```

Candidate validation uses the simplified-model result table:

```bash
python analysis/validate_simplified_candidates.py \
	--matrix data/report.pg_matrix.tsv \
	--screen-results outputs/simplified_mixed_cosinor/simplified_mixed_cosinor_results.tsv \
	--output-dir outputs/simplified_candidate_validation \
	--bootstrap-reps 200
```

For confirmatory work, set the bootstrap count before running and record the
seed and software environment with the results.

## Data and reproducibility

Raw protein matrices, tissue gene lists, generated outputs, and local caches are
excluded from the public repository. To reproduce the analysis, obtain the
input matrix independently and place it at the path supplied to `--input` or
`--matrix`. The scripts write QC tables, model results, bootstrap summaries,
and plots to the requested output directory.

The results are intended for methodological validation and candidate
prioritization. Carryover, acquisition order, missingness, and model-boundary
diagnostics should be reviewed before making biological claims.
