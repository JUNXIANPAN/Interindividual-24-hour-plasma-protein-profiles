"""Literature-search provider port."""

from typing import Any, Protocol


class LiteratureProvider(Protocol):
    async def search(self, query: str) -> list[dict[str, Any]]: ...
