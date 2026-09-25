"""Evidence-store port."""

from typing import Protocol

from rhythm_agent.contracts.evidence import CommonEvidenceEnvelope


class EvidenceStore(Protocol):
    def get(self, evidence_id: str) -> CommonEvidenceEnvelope | None: ...
    def put(self, evidence: CommonEvidenceEnvelope) -> None: ...
