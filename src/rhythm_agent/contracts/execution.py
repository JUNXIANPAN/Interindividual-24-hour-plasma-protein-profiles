"""Branch lifecycle contracts for parallel LangGraph execution."""

from datetime import datetime
from enum import StrEnum

from pydantic import BaseModel, Field

from .artifacts import ArtifactRef
from .diagnostics import Diagnostic
from .evidence import EvidenceRef


class BranchStatus(StrEnum):
    NOT_STARTED = "NOT_STARTED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    NO_EVIDENCE = "NO_EVIDENCE"
    INSUFFICIENT_INFORMATION = "INSUFFICIENT_INFORMATION"
    NOT_EXECUTABLE = "NOT_EXECUTABLE"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"


class BranchExecutionState(BaseModel):
    status: BranchStatus = BranchStatus.NOT_STARTED
    reason_code: str | None = None
    started_at: datetime | None = None
    completed_at: datetime | None = None
    evidence_ref: EvidenceRef | None = None
    artifact_refs: list[ArtifactRef] = Field(default_factory=list)
    diagnostics: list[Diagnostic] = Field(default_factory=list)
