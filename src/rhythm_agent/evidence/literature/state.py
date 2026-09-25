"""Message-oriented LangGraph state for the E subgraph."""

from typing import Annotated, Any

from langchain_core.messages import AnyMessage
from langgraph.graph.message import add_messages
from typing_extensions import TypedDict

from .models import EEvidenceRecord, LiteratureExtraction


class LiteratureResearchState(TypedDict, total=False):
    """E may read deterministic evidence context but cannot mutate A-D records."""

    messages: Annotated[list[AnyMessage], add_messages]
    research_questions: list[str]
    deterministic_evidence_context: str
    retrieved_source_ids: list[str]
    extraction: LiteratureExtraction
    e_evidence: EEvidenceRecord
    iteration_count: int
    max_iterations: int
    warnings: list[str]
    runtime: dict[str, Any]
