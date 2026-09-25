"""Node 1: extract a Protein + Trait request and obtain human approval."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from langgraph.types import interrupt
from pydantic import BaseModel, ConfigDict, Field

from rhythm_agent.contracts.requests import RawResearchRequest
from rhythm_agent.orchestration.state import RhythmResearchState

if TYPE_CHECKING:
    from langchain_core.language_models import BaseChatModel

NODE_ID = "1"
EXTRACTION_PROMPT = (
    "Extract the protein and trait explicitly requested by the user. "
    + "Do not resolve, normalize, or infer biomedical identifiers."
)


class ProteinTraitExtraction(BaseModel):
    """The two required entities extracted from a user's free-form message."""

    protein: str = Field(min_length=1)
    trait: str = Field(min_length=1)

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")


def _review_response(value: Any) -> tuple[bool, str | None]:
    """Normalize boolean and object-shaped interrupt responses."""
    if isinstance(value, bool):
        return value, None
    if isinstance(value, dict):
        approved = bool(value.get("approved", value.get("approve", False)))
        message = value.get("message")
        return approved, message.strip() if isinstance(message, str) and message.strip() else None
    raise TypeError("N01 review response must be a bool or a mapping")


def make_accept_user_input_node(model: BaseChatModel) -> Callable:
    """Build N01 with an injected LLM so the graph remains testable."""
    structured_model = model.with_structured_output(ProteinTraitExtraction)

    def accept_user_input(state: RhythmResearchState) -> dict[str, Any]:
        message = state["message"].strip()
        if not message:
            raise ValueError("message must not be blank")

        while True:
            extracted = structured_model.invoke(
                [
                    ("system", EXTRACTION_PROMPT),
                    ("human", message),
                ]
            )
            request = RawResearchRequest(
                protein=extracted.protein,
                trait=extracted.trait,
                tissue=state.get("tissue"),
                study_id=state.get("study_id"),
            )

            approved, correction = _review_response(
                interrupt(
                    {
                        "kind": "approve_protein_trait",
                        "message": message,
                        "protein": request.protein,
                        "trait": request.trait,
                        "instruction": (
                            "Approve this extraction, or reject it and provide a corrected message."
                        ),
                    }
                )
            )
            if approved:
                return {
                    "message": message,
                    "protein": request.protein,
                    "trait": request.trait,
                    "raw_request": request,
                }
            if correction:
                message = correction
                continue

            correction = interrupt(
                {
                    "kind": "request_correction",
                    "instruction": "Please provide a new message describing the protein and trait.",
                }
            )
            if not isinstance(correction, str) or not correction.strip():
                raise ValueError("A non-empty message is required after rejecting the extraction")
            message = correction.strip()

    return accept_user_input
