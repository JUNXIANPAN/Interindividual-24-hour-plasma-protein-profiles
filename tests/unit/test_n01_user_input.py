"""Node 1 extraction and human-review tests."""

from unittest.mock import Mock, patch

from langgraph.checkpoint.memory import InMemorySaver
from langgraph.types import Command

from rhythm_agent.orchestration.graph import build_graph
from rhythm_agent.orchestration.nodes.n01_user_input import (
    ProteinTraitExtraction,
    make_accept_user_input_node,
)


def _model_with(*extractions: ProteinTraitExtraction) -> tuple[Mock, Mock]:
    structured = Mock()
    structured.invoke.side_effect = list(extractions)
    model = Mock()
    model.with_structured_output.return_value = structured
    return model, structured


def test_n01_extracts_and_writes_approved_protein_and_trait() -> None:
    model, structured = _model_with(ProteinTraitExtraction(protein="P01308", trait="T2D"))
    node = make_accept_user_input_node(model)

    with patch(
        "rhythm_agent.orchestration.nodes.n01_user_input.interrupt", return_value=True
    ) as review:
        result = node({"message": "Does insulin affect type 2 diabetes?"})

    assert result["protein"] == "P01308"
    assert result["trait"] == "T2D"
    assert result["raw_request"].protein == "P01308"
    assert structured.invoke.call_count == 1
    assert review.call_args.args[0]["kind"] == "approve_protein_trait"


def test_n01_reparses_message_attached_to_rejection() -> None:
    model, structured = _model_with(
        ProteinTraitExtraction(protein="CRP", trait="diabetes"),
        ProteinTraitExtraction(protein="CRP", trait="coronary artery disease"),
    )
    node = make_accept_user_input_node(model)

    with patch(
        "rhythm_agent.orchestration.nodes.n01_user_input.interrupt",
        side_effect=[
            {"approved": False, "message": "I meant CRP and coronary artery disease"},
            {"approved": True},
        ],
    ):
        result = node({"message": "Study CRP and diabetes"})

    assert result["trait"] == "coronary artery disease"
    assert result["message"] == "I meant CRP and coronary artery disease"
    assert structured.invoke.call_count == 2


def test_graph_interrupts_then_continues_through_n02_after_approval() -> None:
    extraction = ProteinTraitExtraction(protein="P01308", trait="T2D")
    # LangGraph restarts an interrupted node from its beginning on resume.
    model, _ = _model_with(extraction, extraction)
    graph = build_graph(model=model, checkpointer=InMemorySaver())
    config = {"configurable": {"thread_id": "n01-approval"}}

    paused = graph.invoke({"message": "insulin and type 2 diabetes"}, config)
    assert paused["__interrupt__"][0].value["protein"] == "P01308"

    result = graph.invoke(Command(resume={"approved": True}), config)
    assert result["protein"] == "P01308"
    assert result["trait"] == "T2D"


def test_graph_builds_the_configured_azure_model_by_default() -> None:
    model, _ = _model_with(ProteinTraitExtraction(protein="P01308", trait="T2D"))

    with patch(
        "rhythm_agent.orchestration.graph.build_azure_chat_model", return_value=model
    ) as build_model:
        graph = build_graph()

    assert graph is not None
    build_model.assert_called_once()
