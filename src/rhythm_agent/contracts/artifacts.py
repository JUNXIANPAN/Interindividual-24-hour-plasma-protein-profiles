"""References to large files kept outside LangGraph state."""

from pydantic import BaseModel


class ArtifactRef(BaseModel):
    artifact_id: str
    uri: str
    sha256: str
    kind: str
    schema_version: str
