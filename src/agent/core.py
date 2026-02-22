"""LangChain agent creation and invocation.

Creates a tool-calling agent bound to a dataset, using OpenAI's gpt-4o-mini
and the structured plotting/analysis tools.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

from langchain_core.messages import AIMessage, HumanMessage, SystemMessage, ToolMessage
from langchain_openai import ChatOpenAI

from src.agent.prompts import build_system_prompt
from src.agent.tools import bind_dataset, bind_dataset_state, get_all_tools

if TYPE_CHECKING:
    from anndata import AnnData
    from src.types import DatasetState

logger = logging.getLogger(__name__)


@dataclass
class AgentResponse:
    """Structured response from the agent."""

    text: str
    tool_called: bool


def create_agent(
    adata: AnnData,
    api_key: str,
    model: str = "gpt-4o-mini",
    dataset_state: DatasetState | None = None,
):
    """Create a LangChain tool-calling agent bound to the dataset.

    Args:
        adata: The annotated data matrix.
        api_key: OpenAI API key.
        model: Model name to use.
        dataset_state: Optional processing state tracker.
    """
    bind_dataset(adata)
    if dataset_state is not None:
        bind_dataset_state(dataset_state)

    llm = ChatOpenAI(model=model, api_key=api_key, temperature=0)
    tools = get_all_tools()
    llm_with_tools = llm.bind_tools(tools)
    system_prompt = build_system_prompt(adata, dataset_state=dataset_state)

    return AgentRunner(llm_with_tools, tools, system_prompt)


class AgentRunner:
    """Runs the tool-calling agent loop.

    Handles the LLM → tool call → LLM response cycle, supporting
    multiple sequential tool calls in a single turn.
    """

    def __init__(self, llm_with_tools, tools: list, system_prompt: str) -> None:
        self._llm = llm_with_tools
        self._tools = {t.name: t for t in tools}
        self._system_prompt = system_prompt

    def invoke(
        self,
        user_input: str,
        chat_history: list | None = None,
    ) -> AgentResponse:
        """Run the agent on a user query."""
        messages = [SystemMessage(content=self._system_prompt)]

        if chat_history:
            for role, content in chat_history:
                if role == "user":
                    messages.append(HumanMessage(content=content))
                elif role == "assistant":
                    messages.append(AIMessage(content=content))

        messages.append(HumanMessage(content=user_input))

        tool_called = False
        max_iterations = 5

        for _ in range(max_iterations):
            response = self._llm.invoke(messages)
            messages.append(response)

            if not response.tool_calls:
                break

            tool_called = True
            for tool_call in response.tool_calls:
                tool_name = tool_call["name"]
                tool_args = tool_call["args"]
                logger.info("Tool call: %s(%s)", tool_name, tool_args)

                tool_fn = self._tools.get(tool_name)
                if tool_fn is None:
                    result = f"Error: Unknown tool '{tool_name}'."
                else:
                    try:
                        result = tool_fn.invoke(tool_args)
                    except Exception as e:
                        logger.exception("Tool execution failed: %s", tool_name)
                        result = f"Error executing {tool_name}: {e}"

                messages.append(
                    ToolMessage(content=str(result), tool_call_id=tool_call["id"])
                )

        final_text = response.content if response.content else ""
        return AgentResponse(text=final_text, tool_called=tool_called)
