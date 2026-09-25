"""Artifact-store port."""

from pathlib import Path
from typing import Protocol

from rhythm_agent.contracts.artifacts import ArtifactRef


class ArtifactStore(Protocol):
    def publish(self, path: Path, *, kind: str, schema_version: str) -> ArtifactRef: ...
