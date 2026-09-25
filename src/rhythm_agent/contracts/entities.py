"""Canonical biomedical entities."""

from pydantic import BaseModel, Field


class Gene(BaseModel):
    symbol: str
    ensembl_id: str | None = None
    hgnc_id: str | None = None
    species: str = "Homo sapiens"


class Protein(BaseModel):
    protein_id: str
    namespace: str
    isoform: str | None = None
    species: str = "Homo sapiens"


class Phenotype(BaseModel):
    name: str
    ontology_id: str | None = None
    definition: str | None = None
    synonyms: list[str] = Field(default_factory=list)


class Variant(BaseModel):
    chromosome: str
    position: int
    ref: str
    alt: str
    genome_build: str
    rsid: str | None = None
