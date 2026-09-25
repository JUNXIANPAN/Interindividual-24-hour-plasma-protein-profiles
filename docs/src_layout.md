# Source layout and stable node addresses

Node filenames begin with their design-memo address in an import-safe form:

- `n05_scientific_compute.py` = node `5`
- `r03_map_protein_to_gene.py` = node `R.3`
- `a02_query_store.py` = node `A.2`
- `ca06_freeze_study.py` = node `C.a.6`
- `cb03_validate_manifest.py` = node `C.b.3`
- `d05_validate_manifest.py` = node `D.5`
- `e02_deep_research.py` = node `E.2`
- `s04_guard_causal_language.py` = node `S.4`

Dots are omitted only because ordinary Python imports cannot use dotted filenames. Every
numbered node module also exports its canonical `NODE_ID` value.

Scientific file computation lives outside `src`, under `workflows/`. Python node `5` and
`infrastructure/jobs/snakemake_runner.py` form the control boundary; Snakemake exclusively
owns its internal rule DAG.
