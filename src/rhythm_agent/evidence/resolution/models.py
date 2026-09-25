"""Resolution candidates, mapping outcomes, and human-review payloads."""

from enum import StrEnum


class ResolutionStatus(StrEnum):
    EXACT = "EXACT"
    AMBIGUOUS = "AMBIGUOUS"
    NOT_FOUND = "NOT_FOUND"
