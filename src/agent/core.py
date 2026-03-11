"""LangChain agent creation and invocation.

Creates a tool-calling agent bound to a dataset, using OpenAI's gpt-4o-mini
and the structured plotting/analysis tools.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING

from langchain_core.messages import AIMessage, HumanMessage, SystemMessage, ToolMessage
from langchain_openai import ChatOpenAI

from src.agent.prompts import build_system_prompt
from src.agent.tools import bind_dataset, bind_dataset_state, get_all_tools
from src.logging.service import EventLogger

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
    user_id: str | None = None,
    session_id: str | None = None,
):
    """Create a LangChain tool-calling agent bound to the dataset.

    Args:
        adata: The annotated data matrix.
        api_key: OpenAI API key.
        model: Model name to use.
        dataset_state: Optional processing state tracker.
        user_id: User ID for logging.
        session_id: Session ID for logging.
    """
    bind_dataset(adata)
    if dataset_state is not None:
        bind_dataset_state(dataset_state)

    llm = ChatOpenAI(model=model, api_key=api_key, temperature=0)
    tools = get_all_tools()

    # Log available tools for debugging
    tool_names = [t.name for t in tools]
    logger.info(f"Creating agent with {len(tools)} tools: {', '.join(tool_names)}")

    llm_with_tools = llm.bind_tools(tools)
    system_prompt = build_system_prompt(adata, dataset_state=dataset_state)

    # Initialize event logger if user_id and session_id provided
    event_logger = None
    if user_id and session_id:
        event_logger = EventLogger()

    return AgentRunner(llm_with_tools, tools, system_prompt, event_logger, user_id, session_id)


class AgentRunner:
    """Runs the tool-calling agent loop.

    Handles the LLM → tool call → LLM response cycle, supporting
    multiple sequential tool calls in a single turn.
    """

    def __init__(
        self,
        llm_with_tools,
        tools: list,
        system_prompt: str,
        event_logger: EventLogger | None = None,
        user_id: str | None = None,
        session_id: str | None = None,
    ) -> None:
        self._llm = llm_with_tools
        self._tools = {t.name: t for t in tools}
        self._system_prompt = system_prompt
        self._event_logger = event_logger
        self._user_id = user_id
        self._session_id = session_id

    def invoke(
        self,
        user_input: str,
        chat_history: list | None = None,
    ) -> AgentResponse:
        """Run the agent on a user query."""
        start_time = time.time()

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
        total_tokens = 0

        for _ in range(max_iterations):
            response = self._llm.invoke(messages)
            messages.append(response)

            # Extract token usage if available
            if hasattr(response, "response_metadata"):
                usage = response.response_metadata.get("token_usage", {})
                if usage:
                    total_tokens += usage.get("total_tokens", 0)

            if not response.tool_calls:
                break

            tool_called = True
            for tool_call in response.tool_calls:
                tool_name = tool_call["name"]
                tool_args = tool_call["args"]
                logger.info("Tool call: %s(%s)", tool_name, tool_args)

                tool_start = time.time()
                tool_fn = self._tools.get(tool_name)

                if tool_fn is None:
                    result = f"Error: Unknown tool '{tool_name}'."
                else:
                    try:
                        result = tool_fn.invoke(tool_args)
                        tool_duration = time.time() - tool_start

                        # Log successful tool execution
                        if self._event_logger and self._user_id and self._session_id:
                            self._event_logger.log_tool_execution(
                                user_id=self._user_id,
                                session_id=self._session_id,
                                tool_name=tool_name,
                                args=tool_args,
                                result=str(result),
                                duration_ms=tool_duration * 1000,
                                status="success"
                            )
                    except Exception as e:
                        logger.exception("Tool execution failed: %s", tool_name)
                        result = f"Error executing {tool_name}: {e}"
                        tool_duration = time.time() - tool_start

                        # Log failed tool execution
                        if self._event_logger and self._user_id and self._session_id:
                            self._event_logger.log_tool_execution(
                                user_id=self._user_id,
                                session_id=self._session_id,
                                tool_name=tool_name,
                                args=tool_args,
                                result=str(e),
                                duration_ms=tool_duration * 1000,
                                status="error"
                            )

                messages.append(
                    ToolMessage(content=str(result), tool_call_id=tool_call["id"])
                )

        final_text = response.content if response.content else ""
        response_time = time.time() - start_time

        # Log user message with response time
        if self._event_logger and self._user_id and self._session_id:
            self._event_logger.log_user_message(
                user_id=self._user_id,
                session_id=self._session_id,
                message=user_input,
                response_time_ms=response_time * 1000,
                tool_called=tool_called
            )

        # Log token usage
        if self._event_logger and self._user_id and self._session_id and total_tokens > 0:
            self._event_logger.log_token_usage(
                user_id=self._user_id,
                session_id=self._session_id,
                input_tokens=0,  # Not easily extractable from aggregated usage
                output_tokens=0,
                total_tokens=total_tokens,
                model=self._llm.model_name if hasattr(self._llm, "model_name") else "gpt-4o-mini"
            )

        return AgentResponse(text=final_text, tool_called=tool_called)
