"""Assemble the bounded, tool-using E literature research subgraph."""

from langchain_core.language_models import BaseChatModel
from langchain_core.tools import BaseTool
from langgraph.graph import END, START, StateGraph
from langgraph.prebuilt import ToolNode

from .nodes.e01_build_questions import build_questions
from .nodes.e02_deep_research import make_deep_research_node, route_after_deep_research
from .nodes.e03_extract_evidence import make_extract_evidence_node
from .nodes.e04_build_record import build_record
from .state import LiteratureResearchState

NODE_IDS = tuple(f"E.{index}" for index in range(1, 5))
TOOL_NODE_ID = "E.2.tools"


def build_graph(
    *,
    model: BaseChatModel,
    tools: list[BaseTool],
    checkpointer=None,
):
    """Compile E.1–E.4 while keeping tool execution inside the E.2 loop."""
    if not tools:
        raise ValueError("The E literature graph requires at least one retrieval tool.")

    model_with_tools = model.bind_tools(tools)
    builder = StateGraph(LiteratureResearchState)
    builder.add_node("E.1", build_questions)
    builder.add_node("E.2", make_deep_research_node(model_with_tools))
    builder.add_node(TOOL_NODE_ID, ToolNode(tools))
    builder.add_node("E.3", make_extract_evidence_node(model))
    builder.add_node("E.4", build_record)

    builder.add_edge(START, "E.1")
    builder.add_edge("E.1", "E.2")
    builder.add_conditional_edges(
        "E.2",
        route_after_deep_research,
        {"tools": TOOL_NODE_ID, "extract": "E.3"},
    )
    builder.add_edge(TOOL_NODE_ID, "E.2")
    builder.add_edge("E.3", "E.4")
    builder.add_edge("E.4", END)
    return builder.compile(checkpointer=checkpointer)
