"""Global LangGraph state for candidate-centric evidence evaluation."""

from __future__ import annotations

from typing import Annotated, NotRequired, TypeVar

from typing_extensions import TypedDict

from rhythm_agent.contracts.artifacts import ArtifactRef
from rhythm_agent.contracts.diagnostics import Diagnostic
from rhythm_agent.contracts.evidence import EvidenceRef
from rhythm_agent.contracts.execution import BranchExecutionState
from rhythm_agent.contracts.requests import RawResearchRequest, ResolvedResearchRequest

ValueT = TypeVar("ValueT")


def merge_mappings(left: dict[str, ValueT], right: dict[str, ValueT]) -> dict[str, ValueT]:
    """Merge non-overlapping branch updates into one checkpointable mapping."""
    return {**left, **right}


def append_values(left: list[ValueT], right: list[ValueT]) -> list[ValueT]:
    """Append branch-produced records while preserving their execution order."""
    return [*left, *right]


class RhythmResearchInput(TypedDict):
    """External fields accepted by the top-level graph."""

    message: str
    protein: NotRequired[str | None]
    trait: NotRequired[str | None]

    tissue: NotRequired[str | None]
    study_id: NotRequired[str | None]


class RhythmResearchState(RhythmResearchInput, total=False):
    """Small control state shared by top-level nodes and evidence subgraphs."""

    raw_request: RawResearchRequest
    run_id: str
    resolved_request: ResolvedResearchRequest
    branch_statuses: Annotated[dict[str, BranchExecutionState], merge_mappings]
    evidence_refs: Annotated[dict[str, EvidenceRef], merge_mappings]
    scientific_job_refs: Annotated[list[ArtifactRef], append_values]
    diagnostics: Annotated[list[Diagnostic], append_values]
    review_request_ids: Annotated[list[str], append_values]
    config_hash: str
    protocol_versions: Annotated[dict[str, str], merge_mappings]
    audit_event_refs: Annotated[list[str], append_values]


class RhythmResearchOutput(TypedDict, total=False):
    """Public output while the top-level graph is built incrementally."""

    raw_request: RawResearchRequest
    protein: str
    trait: str
