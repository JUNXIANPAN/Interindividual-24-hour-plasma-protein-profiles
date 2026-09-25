"""Construct a LangChain chat model backed by Azure OpenAI's v1 API."""

from __future__ import annotations

from typing import TYPE_CHECKING

from rhythm_agent.settings import Settings

if TYPE_CHECKING:
    from langchain_core.language_models import BaseChatModel


class ChatModelConfigurationError(RuntimeError):
    """Required Azure chat-model configuration is missing."""


def azure_v1_endpoint(endpoint: str) -> str:
    """Normalize either an Azure resource root or an existing v1 base URL."""
    root = endpoint.rstrip("/")
    if root.endswith("/openai/v1"):
        return f"{root}/"
    return f"{root}/openai/v1/"


def build_azure_chat_model(settings: Settings) -> BaseChatModel:
    """Create the model used by LangGraph nodes, tools, and structured output."""
    missing = [
        name
        for name, value in (
            ("RHYTHM_AGENT_AZURE_OPENAI_ENDPOINT", settings.azure_openai_endpoint),
            ("RHYTHM_AGENT_AZURE_OPENAI_API_KEY", settings.azure_openai_api_key),
            ("RHYTHM_AGENT_AZURE_OPENAI_DEPLOYMENT", settings.azure_openai_deployment),
        )
        if value is None
    ]
    if missing:
        raise ChatModelConfigurationError(
            f"Missing Azure OpenAI settings: {', '.join(missing)}"
        )

    try:
        from langchain_openai import ChatOpenAI
    except ImportError as exc:
        raise ChatModelConfigurationError(
            'Azure LangChain dependencies are not installed; run pip install -e ".[azure]".'
        ) from exc

    endpoint = azure_v1_endpoint(settings.azure_openai_endpoint)
    return ChatOpenAI(
        model=settings.azure_openai_deployment,
        base_url=endpoint,
        api_key=settings.azure_openai_api_key,
        timeout=settings.azure_openai_timeout_seconds,
        max_retries=settings.azure_openai_max_retries,
        use_responses_api=True,
        stream_usage=True,
    )
