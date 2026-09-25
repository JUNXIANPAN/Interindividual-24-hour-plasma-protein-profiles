"""LangChain chat-model construction for configured model providers."""

from .azure_chat_model import build_azure_chat_model

__all__ = ["build_azure_chat_model"]
