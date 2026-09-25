"""Rhythm-genetics release provider port."""

from typing import Any, Protocol


class RhythmGeneticsProvider(Protocol):
    def find_by_gene(self, gene_id: str) -> list[dict[str, Any]]: ...
