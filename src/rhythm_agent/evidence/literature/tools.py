"""Build LangChain tools around provider-neutral literature search ports."""

from langchain_core.tools import BaseTool, StructuredTool

from rhythm_agent.ports.literature_provider import LiteratureProvider


def build_literature_tools(
    *,
    pubmed: LiteratureProvider,
    crossref: LiteratureProvider,
) -> list[BaseTool]:
    """Expose curated literature providers to E.2 without leaking adapters."""

    async def search_pubmed(query: str) -> list[dict]:
        """Search PubMed and return records with stable citation identifiers."""
        return await pubmed.search(query)

    async def search_crossref(query: str) -> list[dict]:
        """Search Crossref for publication metadata and DOI identifiers."""
        return await crossref.search(query)

    return [
        StructuredTool.from_function(
            coroutine=search_pubmed,
            name="search_pubmed",
            description=(
                "Search biomedical publications in PubMed. Use focused queries and retain "
                "PMIDs from every record used as evidence."
            ),
        ),
        StructuredTool.from_function(
            coroutine=search_crossref,
            name="search_crossref",
            description=(
                "Resolve publication metadata and DOI records through Crossref. Use it to "
                "verify citation identity; it is not a substitute for reading evidence."
            ),
        ),
    ]
