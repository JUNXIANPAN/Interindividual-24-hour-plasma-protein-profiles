# Rhythm Research Agent

This project is a candidate-centric biomedical evidence agent built around LangGraph and a
Snakemake scientific-compute boundary.

## Purpose
Given a Protein + Trait pair, the system builds an auditable evidence matrix from:

- R: deterministic entity resolution
- A: released protein-rhythm evidence
- B: optional rhythm-genetics supporting evidence
- C: GWAS study selection and regional evidence
- D: conditional colocalization evidence
- E: agentic literature and mechanism evidence
- S: evidence synthesis and causal-language safeguards

## Source architecture

- `contracts/`: stable Pydantic contracts shared across all layers
- `orchestration/`: top-level LangGraph state, routing, dispatch, and barrier
- `evidence/`: R, A–E, and S evidence-domain workflows
- `ports/`: interfaces required by the application core
- `infrastructure/`: DuckDB, artifact, registry, GWAS, Snakemake, and literature adapters
- `scientific/`: deterministic statistical kernels independent of LangGraph

Executable node files carry their stable memo address, for example
`a02_query_store.py`, `ca03_rank_studies.py`, and `d02_prerequisite_gate.py`.
See `docs/src_layout.md` for the complete naming rule.

Snakemake workflows live in `workflows/`, outside the Python package. Node `5` submits
`ScientificJobSpec` values through `infrastructure/jobs/snakemake_runner.py`; C.b.3 or D.5
then validates the returned manifest.

## Design principles
- A-D are deterministic evidence workflows without LLM dependency.
- E is the primary LLM/deep-research workflow and cannot rewrite A-D results.
- A and C are parallel input-anchored evidence branches; B is optional and D is conditional.
- LangGraph state contains references and control state, not large scientific datasets.
- Snakemake owns deterministic file computation behind `ScientificJobSpec` and result manifests.
- The project uses a src/ layout.
- Candidate evidence is reported independently; no unsupported composite score is produced.

## Status
The v0.2 package structure and core cross-layer contracts are scaffolded. Scientific providers,
workflow nodes, and statistical implementations remain to be implemented and validated.

## Requirements
- Python 3.11
- Modern Python packaging with src/ layout
