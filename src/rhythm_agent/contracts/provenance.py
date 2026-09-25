"""Scientific provenance contracts."""

from datetime import datetime

from pydantic import BaseModel


class Provenance(BaseModel):
    source: str
    version: str
    retrieved_at: datetime
    dataset_id: str | None = None
    checksum: str | None = None
