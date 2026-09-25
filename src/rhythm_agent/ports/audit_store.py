"""Append-only audit-event store port."""

from typing import Any, Protocol


class AuditStore(Protocol):
    def append(self, event: dict[str, Any]) -> str: ...
