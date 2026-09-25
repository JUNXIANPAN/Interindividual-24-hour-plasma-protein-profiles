"""Versioned A evidence record schemas."""

from pydantic import BaseModel, Field


class AEvidenceRecord(BaseModel):
    schema_version: str
    release_id: str
    protein_id: str
    analysis_status: str
    interpretation_code: str
    limitations: list[str] = Field(default_factory=list)
