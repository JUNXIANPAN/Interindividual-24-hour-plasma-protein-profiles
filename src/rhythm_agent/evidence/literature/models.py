"""Structured literature evidence, conflicts, and mechanism hypotheses."""

from enum import StrEnum

from pydantic import BaseModel, Field


class EvidenceStance(StrEnum):
    SUPPORTING = "SUPPORTING"
    CONTRADICTING = "CONTRADICTING"
    CONTEXTUAL = "CONTEXTUAL"
    UNCLEAR = "UNCLEAR"


class LiteratureStatus(StrEnum):
    SUPPORTIVE = "SUPPORTIVE"
    MIXED = "MIXED"
    CONTRADICTORY = "CONTRADICTORY"
    SPARSE = "SPARSE"


class LiteratureEvidenceItem(BaseModel):
    """One claim tied to a retrieved source identifier."""

    citation_id: str
    stance: EvidenceStance
    claim: str
    species: str | None = None
    tissue: str | None = None
    experiment_type: str | None = None
    limitations: list[str] = Field(default_factory=list)


class LiteratureExtraction(BaseModel):
    """Validated output of E.3."""

    overall_status: LiteratureStatus
    evidence: list[LiteratureEvidenceItem] = Field(default_factory=list)
    unresolved_questions: list[str] = Field(default_factory=list)
    mechanism_hypotheses: list[str] = Field(default_factory=list)


class EEvidenceRecord(BaseModel):
    """Final E evidence record assembled without changing A-D."""

    schema_version: str = "1.0"
    status: LiteratureStatus
    evidence: list[LiteratureEvidenceItem] = Field(default_factory=list)
    unresolved_questions: list[str] = Field(default_factory=list)
    mechanism_hypotheses: list[str] = Field(default_factory=list)
    source_ids: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
