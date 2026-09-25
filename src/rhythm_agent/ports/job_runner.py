"""External deterministic scientific-job runner port."""

from typing import Protocol

from rhythm_agent.contracts.jobs import ScientificJobSpec, ScientificResultManifest


class ScientificJobRunner(Protocol):
    def run(self, spec: ScientificJobSpec) -> ScientificResultManifest: ...
