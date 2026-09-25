"""Canonical entity and ontology registry port."""

from typing import Protocol

from rhythm_agent.contracts.entities import Phenotype, Protein


class EntityRegistry(Protocol):
    def resolve_protein(self, value: str) -> list[Protein]: ...
    def resolve_phenotype(self, value: str) -> list[Phenotype]: ...
