"""GWAS metadata and summary-statistics provider port."""

from typing import Any, Protocol


class GWASProvider(Protocol):
    def search_studies(self, phenotype_id: str) -> list[dict[str, Any]]: ...
