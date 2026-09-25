"""Node E.1: Generate focused questions from gaps and conflicts."""

from langchain_core.messages import HumanMessage, SystemMessage

from rhythm_agent.evidence.literature.state import LiteratureResearchState

NODE_ID = "E.1"


def build_questions(state: LiteratureResearchState) -> dict:
    """Seed the message history with bounded questions and immutable A-D context."""
    questions = state.get("research_questions", [])
    if not questions:
        raise ValueError("E.1 requires at least one focused research question.")

    context = state.get("deterministic_evidence_context", "")
    return {
        "messages": [
            SystemMessage(
                content=(
                    "You are a biomedical literature researcher. Use retrieved sources, "
                    "preserve citations, actively seek contradictory evidence, and keep "
                    "association, shared signal, mediation, and causality distinct. "
                    "The supplied A-D statuses are immutable."
                )
            ),
            HumanMessage(
                content=(
                    "Research questions:\n- "
                    + "\n- ".join(questions)
                    + "\n\nDeterministic A-D evidence context:\n"
                    + context
                )
            ),
        ],
        "iteration_count": 0,
        "retrieved_source_ids": [],
    }
