-- Migration 003: logging improvements
-- Applied: 2026-04-07

-- 1. Rename result_preview → result (full content, no truncation)
ALTER TABLE tool_executions RENAME COLUMN result_preview TO result;

-- 2. Add error_stacktrace to tool_executions
ALTER TABLE tool_executions ADD COLUMN error_stacktrace TEXT;

-- 3. Add turn_id and call_index to tool_executions
ALTER TABLE tool_executions ADD COLUMN turn_id TEXT;
ALTER TABLE tool_executions ADD COLUMN call_index INTEGER;

-- 3b. Add turn_id to chat_messages (links user message to its tool calls)
ALTER TABLE chat_messages ADD COLUMN turn_id TEXT;

-- 4. Add dataset_metadata to analysis_sessions
ALTER TABLE analysis_sessions ADD COLUMN dataset_metadata JSONB;

-- 5. Add end_reason to analysis_sessions
ALTER TABLE analysis_sessions ADD COLUMN end_reason TEXT;
