"""Enumerations used by versioned data contracts."""

from enum import StrEnum


class EvidenceType(StrEnum):
    PROTEIN_RHYTHM = "A"
    RHYTHM_GENETICS = "B"
    GWAS_STUDY = "C_STUDY"
    GWAS_REGIONAL = "C_REGIONAL"
    COLOCALIZATION = "D"
    LITERATURE = "E"
