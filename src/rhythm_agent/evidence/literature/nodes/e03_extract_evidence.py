"""Node E.3: Extract validated evidence from the complete research trace."""

from collections.abc import Callable

from langchain_core.language_models import BaseChatModel

from rhythm_agent.evidence.literature.models import LiteratureExtraction
from rhythm_agent.evidence.literature.state import LiteratureResearchState

NODE_ID = "E.3"


def make_extract_evidence_node(model: BaseChatModel) -> Callable:
    """Use LangChain structured output to obtain a validated extraction."""
    structured_model = model.with_structured_output(LiteratureExtraction)

    async def extract_evidence(state: LiteratureResearchState) -> dict:
        extraction = await structured_model.ainvoke(state["messages"])
        return {"extraction": extraction}

    return extract_evidence
