"""Agent loop with three execution paths keyed by router classification.

  - 2a → spec pipeline (extractor → resolver → validator → dispatch → responder → gatekeeper).
         Legacy fallback: if extractor picks a tool not in the new @register REGISTRY
         (or extract raises), drop to the LangChain tool loop (keeps output_guard).
  - 2b → plain LLM text (no bind_tools). output_guard still applies (text-only LLM can
         still fabricate artifact claims until step 11 hardens the 2b system prompt).
  - ambiguous → LangChain tool loop (unchanged from pre-step-10 behavior).

Shared post-logging (user_message, assistant_message, artifacts, tokens, end_session)
runs once per invoke, after the path-specific code returns a _PathResult.
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

import src.domain  # noqa: F401 — import triggers @register side-effects so REGISTRY is populated
from src.agent import artifacts
from src.agent.extractor import extract
from src.agent.gatekeeper import check as gatekeeper_check
from src.agent.output_guard import (
    REPHRASE_FALLBACK_MESSAGE,
    build_corrective_system_message,
    check_output_honesty,
)
from src.agent.prompts import build_system_prompt
from src.agent.responder import draft_response
from src.agent.router import classify_intent
from src.agent.tools import (
    bind_dataset,
    bind_dataset_state,
    bind_logger,
    get_all_tools,
    get_plot_results,
    get_table_results,
)
from src.agent.viz_state import VisualizationState, bind_viz_state
from src.core.registry import REGISTRY, get_tool
from src.core.results import ArtifactResult, TextResult, ToolExecutionError
from src.domain.resolver.resolver import resolve
from src.platform.infra.db.logger import DatabaseLogger
from src.platform.observability.events import EventLogger
from src.spec_validation.result import Issue
from src.spec_validation.validator import validate

if TYPE_CHECKING:
    from anndata import AnnData
    from src.core.types import DatasetState

logger = logging.getLogger(__name__)


@dataclass
class AgentResponse:
    """Structured response from the agent."""

    text: str
    tool_called: bool


@dataclass
class _PathResult:
    """Per-path output consumed by invoke() for shared post-logging."""

    text: str
    tool_called: bool
    end_reason: str
    input_tokens: int = 0
    output_tokens: int = 0


# Per-reason wording for needs_input rendering (Q6: generic MVP template).
_REASON_DETAIL: dict[str, str] = {
    "missing": "value missing",
    "not_found": "not found in the dataset",
    "ambiguous": "ambiguous match",
    "wrong_field": "looks like it belongs to a different field",
    "wrong_type": "wrong type",
    "unknown_tool": "tool not recognized",
    "unknown_scenario": "scenario not recognized",
}


def _classify_tool_status(result: object) -> tuple[str, str | None]:
    """Infer (status, error_msg) from a tool return value.

    Contract: tool error returns must start with ``Error:`` (capital E).
    See ``src/agent/tools.py`` and related tool modules for canonicalized
    prefixes.
    """
    if not isinstance(result, str):
        return ("success", None)
    first_line = result.lstrip().split("\n", 1)[0][:200]
    if first_line.startswith("Error:"):
        return ("error", first_line)
    return ("success", None)


def _render_needs_input(issues: list[Issue]) -> str:
    """Deterministic template for needs_input responses (one line per issue)."""
    lines: list[str] = []
    for issue in issues:
        field_name = issue.field.removeprefix("params.")
        detail = _REASON_DETAIL.get(issue.reason, issue.reason)
        suggestion_str = ""
        if issue.suggestions:
            suggestion_str = f" Suggestions: {', '.join(str(s) for s in issue.suggestions)}."
        lines.append(f"I need clarification on {field_name}: {detail}.{suggestion_str}")
    return "\n".join(lines)


def create_agent(
    adata: "AnnData",
    api_key: str,
    model: str = "gpt-4o-mini",
    dataset_state: "DatasetState | None" = None,
    viz_state: VisualizationState | None = None,
    user_id: str | None = None,
    session_id: str | None = None,
):
    """Create an AgentRunner bound to the dataset + both LLM flavors (with + without tools)."""
    bind_dataset(adata)
    if dataset_state is not None:
        bind_dataset_state(dataset_state)
    if viz_state is not None:
        bind_viz_state(viz_state)

    llm = ChatOpenAI(model=model, api_key=api_key, temperature=0)
    tools = get_all_tools()

    tool_names = [t.name for t in tools]
    logger.info(f"Creating agent with {len(tools)} tools: {', '.join(tool_names)}")

    llm_with_tools = llm.bind_tools(tools)
    viz_block = viz_state.to_prompt_block() if viz_state else ""
    system_prompt = build_system_prompt(adata, dataset_state=dataset_state, viz_state_block=viz_block)

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

    event_logger = None
    db_logger = None
    if user_id and session_id:
        event_logger = EventLogger()
        db_logger = DatabaseLogger()
    bind_logger(event_logger, user_id, session_id)

    return AgentRunner(
        adata=adata,
        llm_with_tools=llm_with_tools,
        llm_plain=llm,
        tools=tools,
        system_prompt=system_prompt,
        event_logger=event_logger,
        db_logger=db_logger,
        user_id=user_id,
        session_id=session_id,
        dataset_metadata=dataset_metadata,
    )


class AgentRunner:
    """Runs the agent across three paths keyed by router classification."""

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
        adata: "AnnData | None" = None,
        llm_plain=None,
    ) -> None:
        self._adata = adata
        self._llm = llm_with_tools
        self._llm_plain = llm_plain
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
        history = chat_history or []

        if self._db_logger and self._user_id and self._session_id:
            self._db_logger.ensure_session(
                self._user_id, self._session_id, filename, self._dataset_metadata
            )

        router_result = classify_intent(user_input)
        logger.info(
            "Router: layer=%s task_type=%s confidence=%s matched_on=%r",
            router_result.layer,
            router_result.task_type,
            router_result.confidence,
            router_result.matched_on,
        )
        self._log_session_event("router_classification", {
            "layer": router_result.layer,
            "task_type": router_result.task_type,
            "confidence": router_result.confidence,
            "matched_on": router_result.matched_on,
        })

        try:
            if router_result.layer == "2a":
                result = self._run_spec_pipeline(user_input, history, turn_id)
                if result is None:
                    result = self._run_langchain_loop(user_input, history, turn_id)
            elif router_result.layer == "2b":
                result = self._run_plain_llm(user_input, history)
            else:
                result = self._run_langchain_loop(user_input, history, turn_id)

            response_time = time.time() - start_time
            self._log_post_turn(
                user_input=user_input,
                result=result,
                response_time=response_time,
                turn_id=turn_id,
            )
            return AgentResponse(text=result.text, tool_called=result.tool_called)

        except Exception:
            if self._db_logger and self._user_id and self._session_id:
                self._db_logger.end_session(self._session_id, "error")
            raise

    # ------------------------------------------------------------------
    # Path: 2a spec pipeline
    # ------------------------------------------------------------------

    def _run_spec_pipeline(
        self,
        user_input: str,
        chat_history: list,
        turn_id: str,
    ) -> _PathResult | None:
        """Run extractor → resolver → validator → dispatch → responder → gatekeeper.

        Returns None to signal a legacy-tool fallback (extractor picked a tool not in
        the new REGISTRY, or extract raised). Caller then runs the LangChain loop.
        """
        try:
            spec = extract(user_input, chat_history)
        except Exception as e:
            logger.warning("Extractor failed, falling back to LangChain loop: %s", e)
            return None

        if spec.tool_name == "none":
            logger.info("Extractor declined (no REGISTRY tool fits); falling back to LangChain")
            self._log_session_event("extractor_declined", {
                "turn_id": turn_id,
                "scenario_id": spec.scenario_id,
                "pre_canonical_params": spec.pre_canonical_params,
            })
            return None

        if spec.tool_name not in REGISTRY:
            logger.info(
                "Extractor picked unknown tool %r; falling back to LangChain loop",
                spec.tool_name,
            )
            return None

        self._log_session_event("spec_emitted", {
            "turn_id": turn_id,
            "scenario_id": spec.scenario_id,
            "tool_name": spec.tool_name,
            "pre_canonical_params": spec.pre_canonical_params,
        })

        try:
            spec, resolver_issues = resolve(spec, self._adata)
        except Exception as e:
            logger.exception("Resolver failed for tool=%s", spec.tool_name)
            return _PathResult(
                text=f"Error resolving request: {e}",
                tool_called=False,
                end_reason="error",
            )

        self._log_session_event("resolver_completed", {
            "turn_id": turn_id,
            "tool_name": spec.tool_name,
            "resolver_issues": [i.reason for i in resolver_issues],
            "canonicalizations": len(spec.canonicalizations_applied),
        })

        if resolver_issues:
            return _PathResult(
                text=_render_needs_input(resolver_issues),
                tool_called=False,
                end_reason="normal",
            )

        spec, val_result = validate(spec)
        self._log_session_event("validator_completed", {
            "turn_id": turn_id,
            "tool_name": spec.tool_name,
            "status": val_result.status,
            "issues": [i.reason for i in val_result.issues],
        })

        if val_result.status == "needs_input":
            if any(i.reason == "unknown_tool" for i in val_result.issues):
                return None
            return _PathResult(
                text=_render_needs_input(val_result.issues),
                tool_called=False,
                end_reason="normal",
            )

        tool_output = self._dispatch_registered(spec, turn_id)

        rr = draft_response(spec, tool_output, chat_history)
        input_tokens = rr.input_tokens
        output_tokens = rr.output_tokens
        self._log_session_event("responder_completed", {
            "turn_id": turn_id,
            "tool_name": spec.tool_name,
            "tokens": rr.input_tokens + rr.output_tokens,
            "retry": False,
        })

        plots = get_plot_results() + artifacts.get_plot_results()
        tables = get_table_results() + artifacts.get_table_results()
        g = gatekeeper_check(spec, tool_output, rr.text, plots, tables)

        final_text = rr.text
        if g.status == "block":
            if g.retry_hint:
                logger.info("Gatekeeper block (%s) — retrying responder with hint", g.reason)
                rr2 = draft_response(
                    spec, tool_output, chat_history, retry_hint=g.retry_hint,
                )
                input_tokens += rr2.input_tokens
                output_tokens += rr2.output_tokens
                self._log_session_event("responder_completed", {
                    "turn_id": turn_id,
                    "tool_name": spec.tool_name,
                    "tokens": rr2.input_tokens + rr2.output_tokens,
                    "retry": True,
                })
                g2 = gatekeeper_check(spec, tool_output, rr2.text, plots, tables)
                if g2.status == "block":
                    logger.warning(
                        "Gatekeeper exhausted after retry (reason=%s); surfacing explanation",
                        g2.reason,
                    )
                    final_text = g2.explanation
                else:
                    final_text = rr2.text
            else:
                logger.warning(
                    "Gatekeeper block (%s) not retry-able; surfacing explanation", g.reason,
                )
                final_text = g.explanation

        return _PathResult(
            text=final_text,
            tool_called=True,
            end_reason="normal",
            input_tokens=input_tokens,
            output_tokens=output_tokens,
        )

    def _dispatch_registered(self, spec, turn_id: str) -> str:
        """Dispatch a @register tool with adata injected positionally. Log execution.

        Contract (L1): @register tools return ToolResult (TextResult | ArtifactResult).
        ArtifactResult bytes are bridged into the legacy artifact channel via
        artifacts.consume_artifact_result so UI consumption stays unchanged during
        migration. _classify_tool_status is NOT called here — it's the legacy
        LangChain-loop convention; the new path uses raw.status directly.
        """
        entry = get_tool(spec.tool_name)
        assert entry is not None  # validator guarantees it's in REGISTRY

        tool_start = time.time()
        error_stacktrace = None
        try:
            raw = entry.callable(self._adata, **spec.params)
            assert isinstance(raw, (TextResult, ArtifactResult)), (
                f"Tool {spec.tool_name!r} returned {type(raw).__name__}, expected ToolResult"
            )
            if isinstance(raw, ArtifactResult):
                artifacts.consume_artifact_result(raw)
            tool_output = raw.text
            status = raw.status
            error_msg = raw.error_message
            tool_duration = time.time() - tool_start
        except ToolExecutionError as e:
            tool_duration = time.time() - tool_start
            tool_output = f"Error: {e}"
            status = "error"
            error_msg = str(e)[:200]
        except Exception as e:
            logger.exception("Registered-tool dispatch failed: %s", spec.tool_name)
            tool_duration = time.time() - tool_start
            tool_output = f"Error executing {spec.tool_name}: {e}"
            status = "error"
            error_msg = tool_output[:200]
            error_stacktrace = traceback.format_exc()

        if self._event_logger and self._user_id and self._session_id:
            self._event_logger.log_tool_execution(
                user_id=self._user_id,
                session_id=self._session_id,
                tool_name=spec.tool_name,
                args=spec.params,
                result=tool_output,
                duration_ms=tool_duration * 1000,
                status=status,
                turn_id=turn_id,
                call_index=0,
                error=error_msg,
                **({"error_stacktrace": error_stacktrace} if error_stacktrace else {}),
            )
        if self._db_logger and self._user_id and self._session_id:
            self._db_logger.log_tool_execution(
                user_id=self._user_id,
                session_id=self._session_id,
                tool_name=spec.tool_name,
                args=spec.params,
                result=tool_output,
                duration_ms=tool_duration * 1000,
                status=status,
                turn_id=turn_id,
                call_index=0,
                **({"error_stacktrace": error_stacktrace} if error_stacktrace else {}),
            )
        return tool_output

    # ------------------------------------------------------------------
    # Path: 2b plain LLM
    # ------------------------------------------------------------------

    def _run_plain_llm(
        self,
        user_input: str,
        chat_history: list,
    ) -> _PathResult:
        """Plain LLM call (no tools) with output_guard one-shot retry.

        output_guard still applies because a text-only LLM can still fabricate
        artifact-creation claims. Will be redundant once step 11 hardens the
        2b-specific system prompt.
        """
        messages = self._build_messages(user_input, chat_history)
        input_tokens = 0
        output_tokens = 0

        response = self._llm_plain.invoke(messages)
        input_tokens, output_tokens = self._accumulate_tokens(
            response, input_tokens, output_tokens,
        )
        text = response.content or ""
        guard = check_output_honesty(text, [])
        if not guard.triggered:
            return _PathResult(
                text=text, tool_called=False, end_reason="normal",
                input_tokens=input_tokens, output_tokens=output_tokens,
            )

        # One-shot corrective retry
        logger.warning(
            "2b output_guard triggered (attempt 1): matched=%r", guard.matched_phrase,
        )
        self._log_session_event("output_guard_triggered", {
            "attempt": 1,
            "matched_phrase": guard.matched_phrase,
            "assistant_text_excerpt": text[:500],
            "path": "2b",
        })
        messages.append(response)
        messages.append(SystemMessage(
            content=build_corrective_system_message(guard.matched_phrase or ""),
        ))
        response2 = self._llm_plain.invoke(messages)
        input_tokens, output_tokens = self._accumulate_tokens(
            response2, input_tokens, output_tokens,
        )
        text2 = response2.content or ""
        guard2 = check_output_honesty(text2, [])
        if not guard2.triggered:
            return _PathResult(
                text=text2, tool_called=False, end_reason="normal",
                input_tokens=input_tokens, output_tokens=output_tokens,
            )

        logger.warning(
            "2b output_guard exhausted (attempt 2): matched=%r", guard2.matched_phrase,
        )
        self._log_session_event("output_guard_exhausted", {
            "attempt": 2,
            "matched_phrase": guard2.matched_phrase,
            "first_assistant_text_excerpt": text[:500],
            "second_assistant_text_excerpt": text2[:500],
            "path": "2b",
        })
        return _PathResult(
            text=REPHRASE_FALLBACK_MESSAGE,
            tool_called=False,
            end_reason="normal",
            input_tokens=input_tokens,
            output_tokens=output_tokens,
        )

    # ------------------------------------------------------------------
    # Path: legacy LangChain loop (ambiguous + 2a-fallback)
    # ------------------------------------------------------------------

    def _run_langchain_loop(
        self,
        user_input: str,
        chat_history: list,
        turn_id: str,
    ) -> _PathResult:
        """Legacy LangChain tool-calling loop with output_guard. Behavior preserved.

        Used for: router=ambiguous, AND 2a fallback when the extractor picked a
        legacy @tool not in the new REGISTRY.
        """
        messages = self._build_messages(user_input, chat_history)

        tool_called = False
        call_index = -1
        max_iterations = 5
        hit_max_iterations = True
        input_tokens = 0
        output_tokens = 0

        guard_attempts = 0
        first_hallucinated_excerpt: str | None = None
        guard_exhausted = False
        response = None

        for _ in range(max_iterations):
            response = self._llm.invoke(messages)
            messages.append(response)
            input_tokens, output_tokens = self._accumulate_tokens(
                response, input_tokens, output_tokens,
            )

            if not response.tool_calls:
                if tool_called:
                    hit_max_iterations = False
                    break

                guard = check_output_honesty(
                    response.content or "", response.tool_calls,
                )
                if not guard.triggered:
                    hit_max_iterations = False
                    break

                excerpt = (response.content or "")[:500]
                if guard_attempts == 0:
                    logger.warning(
                        "Output guard triggered (attempt 1): matched=%r",
                        guard.matched_phrase,
                    )
                    self._log_session_event("output_guard_triggered", {
                        "attempt": 1,
                        "matched_phrase": guard.matched_phrase,
                        "assistant_text_excerpt": excerpt,
                        "turn_id": turn_id,
                    })
                    first_hallucinated_excerpt = excerpt
                    guard_attempts += 1
                    messages.append(SystemMessage(
                        content=build_corrective_system_message(
                            guard.matched_phrase or "",
                        ),
                    ))
                    continue

                logger.warning(
                    "Output guard exhausted (attempt 2): matched=%r",
                    guard.matched_phrase,
                )
                self._log_session_event("output_guard_exhausted", {
                    "attempt": 2,
                    "matched_phrase": guard.matched_phrase,
                    "first_assistant_text_excerpt": first_hallucinated_excerpt,
                    "second_assistant_text_excerpt": excerpt,
                    "turn_id": turn_id,
                })
                guard_exhausted = True
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
                        status, error_msg = _classify_tool_status(result)

                        if self._event_logger and self._user_id and self._session_id:
                            self._event_logger.log_tool_execution(
                                user_id=self._user_id,
                                session_id=self._session_id,
                                tool_name=tool_name,
                                args=tool_args,
                                result=str(result),
                                duration_ms=tool_duration * 1000,
                                status=status,
                                turn_id=turn_id,
                                call_index=call_index,
                                error=error_msg,
                            )
                        if self._db_logger and self._user_id and self._session_id:
                            self._db_logger.log_tool_execution(
                                user_id=self._user_id,
                                session_id=self._session_id,
                                tool_name=tool_name,
                                args=tool_args,
                                result=str(result),
                                duration_ms=tool_duration * 1000,
                                status=status,
                                turn_id=turn_id,
                                call_index=call_index,
                            )
                    except Exception as e:
                        logger.exception("Tool execution failed: %s", tool_name)
                        error_tb = traceback.format_exc()
                        result = f"Error executing {tool_name}: {e}"
                        tool_duration = time.time() - tool_start

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
                    ToolMessage(content=str(result), tool_call_id=tool_call["id"]),
                )

        if guard_exhausted:
            final_text = REPHRASE_FALLBACK_MESSAGE
        else:
            final_text = response.content if response and response.content else ""
        end_reason = "max_iterations" if hit_max_iterations else "normal"
        return _PathResult(
            text=final_text,
            tool_called=tool_called,
            end_reason=end_reason,
            input_tokens=input_tokens,
            output_tokens=output_tokens,
        )

    # ------------------------------------------------------------------
    # Shared helpers
    # ------------------------------------------------------------------

    def _build_messages(
        self, user_input: str, chat_history: list,
    ) -> list:
        """Construct LangChain message list from history tuples + current user input."""
        messages = [SystemMessage(content=self._system_prompt)]
        for role, content in chat_history:
            if role == "user":
                messages.append(HumanMessage(content=content))
            elif role == "assistant":
                messages.append(AIMessage(content=content))
        messages.append(HumanMessage(content=user_input))
        return messages

    @staticmethod
    def _accumulate_tokens(
        response, input_tokens: int, output_tokens: int,
    ) -> tuple[int, int]:
        """Pull prompt/completion tokens off a LangChain response metadata blob."""
        if hasattr(response, "response_metadata"):
            usage = response.response_metadata.get("token_usage", {})
            if usage:
                input_tokens += usage.get("prompt_tokens", 0)
                output_tokens += usage.get("completion_tokens", 0)
        return input_tokens, output_tokens

    def _log_session_event(self, event_type: str, metadata: dict) -> None:
        """Fire-and-forget EventLogger session event; no-op when logger unbound."""
        if self._event_logger and self._user_id and self._session_id:
            self._event_logger.log_session_event(
                user_id=self._user_id,
                session_id=self._session_id,
                event_type=event_type,
                metadata=metadata,
            )

    def _log_post_turn(
        self,
        user_input: str,
        result: _PathResult,
        response_time: float,
        turn_id: str,
    ) -> None:
        """Shared logging suffix: user_message, assistant_message, artifacts, tokens, end_session."""
        if self._event_logger and self._user_id and self._session_id:
            self._event_logger.log_user_message(
                user_id=self._user_id,
                session_id=self._session_id,
                message=user_input,
                response_time_ms=response_time * 1000,
                tool_called=result.tool_called,
                turn_id=turn_id,
            )
        if self._db_logger and self._user_id and self._session_id:
            self._db_logger.log_user_message(
                user_id=self._user_id,
                session_id=self._session_id,
                message=user_input,
                response_time_ms=response_time * 1000,
                tool_called=result.tool_called,
                turn_id=turn_id,
            )

        if self._db_logger and self._user_id and self._session_id:
            msg_id = self._db_logger.log_assistant_message(
                user_id=self._user_id,
                session_id=self._session_id,
                message=result.text or "(no text response)",
                turn_id=turn_id,
            )
            if msg_id is not None:
                plots = get_plot_results() + artifacts.get_plot_results()
                tables = get_table_results() + artifacts.get_table_results()
                logger.info(
                    "Artifact logging: msg_id=%s plots=%d tables=%d",
                    msg_id, len(plots), len(tables),
                )
                self._db_logger.log_artifacts(
                    message_id=msg_id,
                    session_id=self._session_id,
                    user_id=self._user_id,
                    plot_results=plots,
                    table_results=tables,
                )

        total_tokens = result.input_tokens + result.output_tokens
        model_name = getattr(self._llm, "model_name", "gpt-4o-mini")
        if self._event_logger and self._user_id and self._session_id and total_tokens > 0:
            self._event_logger.log_token_usage(
                user_id=self._user_id,
                session_id=self._session_id,
                input_tokens=result.input_tokens,
                output_tokens=result.output_tokens,
                model=model_name,
            )
        if self._db_logger and self._user_id and self._session_id and total_tokens > 0:
            self._db_logger.log_token_usage(
                user_id=self._user_id,
                session_id=self._session_id,
                input_tokens=result.input_tokens,
                output_tokens=result.output_tokens,
                model=model_name,
            )

        if self._db_logger and self._user_id and self._session_id:
            self._db_logger.end_session(self._session_id, result.end_reason)
