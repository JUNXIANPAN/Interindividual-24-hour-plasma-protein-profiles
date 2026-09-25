"""Node E.4: Build literature evidence without changing A-D statuses."""

from rhythm_agent.evidence.literature.models import EEvidenceRecord
from rhythm_agent.evidence.literature.state import LiteratureResearchState

NODE_ID = "E.4"


def build_record(state: LiteratureResearchState) -> dict:
    """Convert the validated E.3 result into the final E evidence record."""
    extraction = state["extraction"]
    record = EEvidenceRecord(
        status=extraction.overall_status,
        evidence=extraction.evidence,
        unresolved_questions=extraction.unresolved_questions,
        mechanism_hypotheses=extraction.mechanism_hypotheses,
        source_ids=sorted({item.citation_id for item in extraction.evidence}),
        limitations=[
            "Literature evidence does not modify deterministic A-D results.",
            "Publication and retrieval bias may affect the available evidence.",
        ],
    )
    return {"e_evidence": record}
