"""Portable warning and error records."""

from enum import StrEnum

from pydantic import BaseModel


class Severity(StrEnum):
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


class Diagnostic(BaseModel):
    code: str
    message: str
    severity: Severity
    source_node: str
    recoverable: bool = False
