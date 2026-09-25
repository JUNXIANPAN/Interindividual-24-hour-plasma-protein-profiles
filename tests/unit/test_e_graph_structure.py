"""Offline structure checks for the native LangGraph E subgraph."""

from unittest import TestCase

from rhythm_agent.evidence.literature.graph import TOOL_NODE_ID, build_graph
from rhythm_agent.evidence.literature.tools import build_literature_tools
from rhythm_agent.infrastructure.llm.azure_chat_model import (
    azure_v1_endpoint,
    build_azure_chat_model,
)
from rhythm_agent.settings import Settings


class _FakeLiteratureProvider:
    async def search(self, query: str) -> list[dict]:
        return [{"id": "test-source", "query": query}]


class ELiteratureGraphStructureTest(TestCase):
    def test_compiles_numbered_nodes_and_tool_loop(self) -> None:
        settings = Settings(
            _env_file=None,
            azure_openai_endpoint="https://example.openai.azure.com",
            azure_openai_api_key="not-a-real-key",
            azure_openai_deployment="test-deployment",
        )
        model = build_azure_chat_model(settings)
        provider = _FakeLiteratureProvider()
        tools = build_literature_tools(pubmed=provider, crossref=provider)

        graph = build_graph(model=model, tools=tools)
        node_names = set(graph.get_graph().nodes)

        self.assertTrue({"E.1", "E.2", TOOL_NODE_ID, "E.3", "E.4"} <= node_names)

    def test_v1_endpoint_is_appended_only_once(self) -> None:
        expected = "https://example.openai.azure.com/openai/v1/"
        self.assertEqual(azure_v1_endpoint("https://example.openai.azure.com"), expected)
        self.assertEqual(azure_v1_endpoint(expected), expected)
