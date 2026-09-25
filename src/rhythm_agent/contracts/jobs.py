"""Contract at the LangGraph-to-Snakemake boundary."""

from enum import StrEnum

from pydantic import BaseModel, Field

from .artifacts import ArtifactRef
from .regions import GenomicRegion


class JobStatus(StrEnum):
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    PARTIAL = "PARTIAL"


class ScientificJobSpec(BaseModel):
    run_id: str
    job_id: str
    workflow_name: str
    protein_id: str | None = None
    gene_id: str | None = None
    phenotype_id: str | None = None
    selected_gwas_study_id: str | None = None
    regions: list[GenomicRegion] = Field(default_factory=list)
    input_artifact_refs: list[ArtifactRef] = Field(default_factory=list)
    protocol_version: str
    config_hash: str
    code_version: str
    random_seed: int


class ScientificResultManifest(BaseModel):
    schema_version: str
    run_id: str
    job_id: str
    workflow_name: str
    status: JobStatus
    inputs: list[ArtifactRef] = Field(default_factory=list)
    outputs: list[ArtifactRef] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    metrics: dict[str, float] = Field(default_factory=dict)
