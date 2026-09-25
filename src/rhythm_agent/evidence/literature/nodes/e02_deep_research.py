"""Node E.2: Let a LangChain chat model reason and request literature tools."""

from collections.abc import Callable

from langchain_core.language_models import BaseChatModel

from rhythm_agent.evidence.literature.state import LiteratureResearchState

NODE_ID = "E.2"


def make_deep_research_node(model_with_tools: BaseChatModel) -> Callable:
    """Create the E.2 node around a tool-bound LangChain chat model."""

    async def deep_research(state: LiteratureResearchState) -> dict:
        response = await model_with_tools.ainvoke(state["messages"])
        return {
            "messages": [response],
            "iteration_count": state.get("iteration_count", 0) + 1,
        }

    return deep_research


def route_after_deep_research(state: LiteratureResearchState) -> str:
    """Continue tool execution only while calls exist and the loop is bounded."""
    last_message = state["messages"][-1]
    tool_calls = getattr(last_message, "tool_calls", [])
    if tool_calls and state.get("iteration_count", 0) < state.get("max_iterations", 6):
        return "tools"
    return "extract"
