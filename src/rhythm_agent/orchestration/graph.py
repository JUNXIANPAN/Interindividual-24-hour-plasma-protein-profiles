"""Assembly point for global nodes 1–7."""

from __future__ import annotations

from typing import TYPE_CHECKING

from langgraph.graph import END, START, StateGraph

from rhythm_agent.infrastructure.llm.azure_chat_model import build_azure_chat_model
from rhythm_agent.orchestration.nodes.n01_user_input import make_accept_user_input_node
from rhythm_agent.orchestration.nodes.n02_initialize_run import initialize_run
from rhythm_agent.orchestration.state import (
    RhythmResearchInput,
    RhythmResearchOutput,
    RhythmResearchState,
)
from rhythm_agent.settings import Settings

if TYPE_CHECKING:
    from langchain_core.language_models import BaseChatModel


def build_graph(
    *,
    model: BaseChatModel | None = None,
    settings: Settings | None = None,
    checkpointer=None,
):
    """Compile START → N01 → N02, using the configured Azure LLM by default."""
    if model is None:
        model = build_azure_chat_model(settings or Settings())

    builder = StateGraph(
        RhythmResearchState,
        input_schema=RhythmResearchInput,
        output_schema=RhythmResearchOutput,
    )
    builder.add_node("1", make_accept_user_input_node(model))
    builder.add_node("2", initialize_run)
    builder.add_edge(START, "1")
    builder.add_edge("1", "2")
    builder.add_edge("2", END)
    return builder.compile(checkpointer=checkpointer)
