"""Raw and resolved candidate-evaluation requests."""

from pydantic import BaseModel, ConfigDict, Field

from .entities import Gene, Phenotype, Protein
from .provenance import Provenance


class RawResearchRequest(BaseModel):
    """Validated structural input; it contains no scientific interpretation."""

    protein: str = Field(min_length=1)
    trait: str = Field(min_length=1)
    ancestry: str | None = None
    tissue: str | None = None
    study_id: str | None = None

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")


class ResolvedResearchRequest(BaseModel):
    request_id: str
    protein: Protein
    genes: list[Gene]
    phenotype: Phenotype
    user_constraints: dict[str, str] = Field(default_factory=dict)
    resolution_provenance: list[Provenance] = Field(default_factory=list)
