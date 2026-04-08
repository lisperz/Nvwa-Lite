"""LangChain agent creation and invocation.

Creates a tool-calling agent bound to a dataset, using OpenAI's gpt-4o-mini
and the structured plotting/analysis tools.
"""

from __future__ import annotations

import logging
import time
import traceback
import uuid
from dataclasses import dataclass
from typing import TYPE_CHECKING

from langchain_core.messages import AIMessage, HumanMessage, SystemMessage, ToolMessage
from langchain_openai import ChatOpenAI

from src.agent.prompts import build_system_prompt
from src.agent.router import classify_intent
from src.agent.tools import bind_dataset, bind_dataset_state, get_all_tools
from src.agent.viz_state import VisualizationState, bind_viz_state
from src.db.logger import DatabaseLogger
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
    viz_state: VisualizationState | None = None,
    user_id: str | None = None,
    session_id: str | None = None,
):
    """Create a LangChain tool-calling agent bound to the dataset.

    Args:
        adata: The annotated data matrix.
        api_key: OpenAI API key.
        model: Model name to use.
        dataset_state: Optional processing state tracker.
        viz_state: Optional visualization state for multi-turn refinement.
        user_id: User ID for logging.
        session_id: Session ID for logging.
    """
    bind_dataset(adata)
    if dataset_state is not None:
        bind_dataset_state(dataset_state)
    if viz_state is not None:
        bind_viz_state(viz_state)

    llm = ChatOpenAI(model=model, api_key=api_key, temperature=0)
    tools = get_all_tools()

    # Log available tools for debugging
    tool_names = [t.name for t in tools]
    logger.info(f"Creating agent with {len(tools)} tools: {', '.join(tool_names)}")

    llm_with_tools = llm.bind_tools(tools)
    viz_block = viz_state.to_prompt_block() if viz_state else ""
    system_prompt = build_system_prompt(adata, dataset_state=dataset_state, viz_state_block=viz_block)

    # Extract dataset metadata for session-start logging
    clustering_key = next(
        (k for k in adata.obs.columns if "leiden" in k or "louvain" in k), None
    )
    dataset_metadata = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "has_umap": "X_umap" in adata.obsm,
        "has_clustering": clustering_key is not None,
        "clustering_key": clustering_key,
    }

    # Initialize event logger if user_id and session_id provided
    event_logger = None
    db_logger = None
    if user_id and session_id:
        event_logger = EventLogger()
        db_logger = DatabaseLogger()

    return AgentRunner(
        llm_with_tools, tools, system_prompt, event_logger, db_logger,
        user_id, session_id, dataset_metadata,
    )


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
        db_logger: DatabaseLogger | None = None,
        user_id: str | None = None,
        session_id: str | None = None,
        dataset_metadata: dict | None = None,
    ) -> None:
        self._llm = llm_with_tools
        self._tools = {t.name: t for t in tools}
        self._system_prompt = system_prompt
        self._event_logger = event_logger
        self._db_logger = db_logger
        self._user_id = user_id
        self._session_id = session_id
        self._dataset_metadata = dataset_metadata

    def invoke(
        self,
        user_input: str,
        chat_history: list | None = None,
        filename: str = "unknown",
    ) -> AgentResponse:
        """Run the agent on a user query."""
        start_time = time.time()
        turn_id = str(uuid.uuid4())

        # Ensure DB session row exists before logging anything
        if self._db_logger and self._user_id and self._session_id:
            self._db_logger.ensure_session(
                self._user_id, self._session_id, filename, self._dataset_metadata
            )

        messages = [SystemMessage(content=self._system_prompt)]

        if chat_history:
            for role, content in chat_history:
                if role == "user":
                    messages.append(HumanMessage(content=content))
                elif role == "assistant":
                    messages.append(AIMessage(content=content))

        messages.append(HumanMessage(content=user_input))

        # Classify intent — logging only, no behavior change (router v0)
        router_result = classify_intent(user_input)
        logger.info(
            "Router: layer=%s task_type=%s confidence=%s matched_on=%r",
            router_result.layer,
            router_result.task_type,
            router_result.confidence,
            router_result.matched_on,
        )
        if self._event_logger and self._user_id and self._session_id:
            self._event_logger.log_session_event(
                user_id=self._user_id,
                session_id=self._session_id,
                event_type="router_classification",
                metadata={
                    "layer": router_result.layer,
                    "task_type": router_result.task_type,
                    "confidence": router_result.confidence,
                    "matched_on": router_result.matched_on,
                },
            )

        tool_called = False
        call_index = 0
        max_iterations = 5
        hit_max_iterations = True
        total_input_tokens = 0
        total_output_tokens = 0

        try:
            for _ in range(max_iterations):
                response = self._llm.invoke(messages)
                messages.append(response)

                # Extract token usage if available
                if hasattr(response, "response_metadata"):
                    usage = response.response_metadata.get("token_usage", {})
                    if usage:
                        total_input_tokens += usage.get("prompt_tokens", 0)
                        total_output_tokens += usage.get("completion_tokens", 0)

                if not response.tool_calls:
                    hit_max_iterations = False
                    break

                tool_called = True
                for tool_call in response.tool_calls:
                    tool_name = tool_call["name"]
                    tool_args = tool_call["args"]
                    call_index += 1
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
                                    status="success",
                                    turn_id=turn_id,
                                    call_index=call_index,
                                )
                            if self._db_logger and self._user_id and self._session_id:
                                self._db_logger.log_tool_execution(
                                    user_id=self._user_id,
                                    session_id=self._session_id,
                                    tool_name=tool_name,
                                    args=tool_args,
                                    result=str(result),
                                    duration_ms=tool_duration * 1000,
                                    status="success",
                                    turn_id=turn_id,
                                    call_index=call_index,
                                )
                        except Exception as e:
                            logger.exception("Tool execution failed: %s", tool_name)
                            error_tb = traceback.format_exc()
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
                                    status="error",
                                    error_stacktrace=error_tb,
                                    turn_id=turn_id,
                                    call_index=call_index,
                                )
                            if self._db_logger and self._user_id and self._session_id:
                                self._db_logger.log_tool_execution(
                                    user_id=self._user_id,
                                    session_id=self._session_id,
                                    tool_name=tool_name,
                                    args=tool_args,
                                    result=str(e),
                                    duration_ms=tool_duration * 1000,
                                    status="error",
                                    error_stacktrace=error_tb,
                                    turn_id=turn_id,
                                    call_index=call_index,
                                )

                    messages.append(
                        ToolMessage(content=str(result), tool_call_id=tool_call["id"])
                    )

            final_text = response.content if response.content else ""
            response_time = time.time() - start_time
            end_reason = "max_iterations" if hit_max_iterations else "normal"

            # Log user message with response time
            if self._event_logger and self._user_id and self._session_id:
                self._event_logger.log_user_message(
                    user_id=self._user_id,
                    session_id=self._session_id,
                    message=user_input,
                    response_time_ms=response_time * 1000,
                    tool_called=tool_called,
                    turn_id=turn_id,
                )
            if self._db_logger and self._user_id and self._session_id:
                self._db_logger.log_user_message(
                    user_id=self._user_id,
                    session_id=self._session_id,
                    message=user_input,
                    response_time_ms=response_time * 1000,
                    tool_called=tool_called,
                    turn_id=turn_id,
                )

            # Log assistant response (always, even if empty)
            if self._db_logger and self._user_id and self._session_id:
                msg_id = self._db_logger.log_assistant_message(
                    user_id=self._user_id,
                    session_id=self._session_id,
                    message=final_text or "(no text response)",
                )
                # Persist any plots/tables generated this turn
                if msg_id is not None:
                    from src.agent.tools import get_plot_results, get_table_results
                    plots = get_plot_results()
                    tables = get_table_results()
                    logger.info("Artifact logging: msg_id=%s plots=%d tables=%d", msg_id, len(plots), len(tables))
                    self._db_logger.log_artifacts(
                        message_id=msg_id,
                        session_id=self._session_id,
                        user_id=self._user_id,
                        plot_results=plots,
                        table_results=tables,
                    )

            # Log token usage
            total_tokens = total_input_tokens + total_output_tokens
            if self._event_logger and self._user_id and self._session_id and total_tokens > 0:
                self._event_logger.log_token_usage(
                    user_id=self._user_id,
                    session_id=self._session_id,
                    input_tokens=total_input_tokens,
                    output_tokens=total_output_tokens,
                    model=self._llm.model_name if hasattr(self._llm, "model_name") else "gpt-4o-mini"
                )
            if self._db_logger and self._user_id and self._session_id and total_tokens > 0:
                self._db_logger.log_token_usage(
                    user_id=self._user_id,
                    session_id=self._session_id,
                    input_tokens=total_input_tokens,
                    output_tokens=total_output_tokens,
                    model=self._llm.model_name if hasattr(self._llm, "model_name") else "gpt-4o-mini"
                )

            if self._db_logger and self._user_id and self._session_id:
                self._db_logger.end_session(self._session_id, end_reason)

            return AgentResponse(text=final_text, tool_called=tool_called)

        except Exception:
            if self._db_logger and self._user_id and self._session_id:
                self._db_logger.end_session(self._session_id, "error")
            raise
