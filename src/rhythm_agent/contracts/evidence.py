"""Common envelope and references for all evidence domains."""

from typing import Any

from pydantic import BaseModel, Field

from .artifacts import ArtifactRef
from .provenance import Provenance


class EvidenceRef(BaseModel):
    evidence_id: str
    evidence_type: str
    schema_version: str


class CommonEvidenceEnvelope(BaseModel):
    evidence_id: str
    evidence_type: str
    schema_version: str
    subject: dict[str, Any]
    status: str
    quality_grade: str | None = None
    interpretation_code: str | None = None
    limitations: list[str] = Field(default_factory=list)
    provenance: list[Provenance] = Field(default_factory=list)
    artifact_refs: list[ArtifactRef] = Field(default_factory=list)
