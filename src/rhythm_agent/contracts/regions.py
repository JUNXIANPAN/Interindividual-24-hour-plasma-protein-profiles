"""Genomic region contracts."""

from pydantic import BaseModel


class GenomicRegion(BaseModel):
    chromosome: str
    start: int
    end: int
    genome_build: str
    definition_method: str
