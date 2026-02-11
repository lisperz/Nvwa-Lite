"""LangChain agent creation and invocation.

Creates a tool-calling agent bound to a dataset, using OpenAI's gpt-4o
and the structured plotting tools.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

from langchain_core.messages import AIMessage, HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI

from src.agent.prompts import build_system_prompt
from src.agent.tools import bind_dataset, get_all_tools

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


@dataclass
class AgentResponse:
    """Structured response from the agent."""

    text: str
    tool_called: bool


def create_agent(adata: AnnData, api_key: str, model: str = "gpt-4o"):
    """Create a LangChain tool-calling agent bound to the dataset.

    Args:
        adata: The annotated data matrix.
        api_key: OpenAI API key.
        model: Model name to use.

    Returns:
        A callable agent (invoke function).
    """
    bind_dataset(adata)

    llm = ChatOpenAI(
        model=model,
        api_key=api_key,
        temperature=0,
    )

    tools = get_all_tools()
    llm_with_tools = llm.bind_tools(tools)
    system_prompt = build_system_prompt(adata)

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
        """Run the agent on a user query.

        Args:
            user_input: The user's natural language query.
            chat_history: Optional list of previous (role, content) tuples.

        Returns:
            AgentResponse with the final text and whether a tool was called.
        """
        messages = [SystemMessage(content=self._system_prompt)]

        # Add chat history
        if chat_history:
            for role, content in chat_history:
                if role == "user":
                    messages.append(HumanMessage(content=content))
                elif role == "assistant":
                    messages.append(AIMessage(content=content))

        messages.append(HumanMessage(content=user_input))

        tool_called = False
        max_iterations = 5  # Safety limit

        for _ in range(max_iterations):
            response = self._llm.invoke(messages)
            messages.append(response)

            # If no tool calls, we're done
            if not response.tool_calls:
                break

            tool_called = True

            # Execute each tool call
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

                # Add tool result as a ToolMessage
                from langchain_core.messages import ToolMessage
                messages.append(
                    ToolMessage(content=str(result), tool_call_id=tool_call["id"])
                )

        # Extract final text
        final_text = response.content if response.content else ""
        return AgentResponse(text=final_text, tool_called=tool_called)
