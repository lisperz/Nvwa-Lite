-- Migration: Add feedback_responses table
-- Created: 2026-03-28
-- Purpose: Store user feedback after analysis sessions

CREATE TABLE IF NOT EXISTS feedback_responses (
    response_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id TEXT NOT NULL,
    user_id TEXT NOT NULL,
    q1_score INTEGER NOT NULL CHECK (q1_score BETWEEN 1 AND 5),
    q2_time_saved VARCHAR(50) NOT NULL,
    q3_open_text TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- Index for querying by session and user
CREATE INDEX IF NOT EXISTS idx_feedback_session ON feedback_responses(session_id);
CREATE INDEX IF NOT EXISTS idx_feedback_user ON feedback_responses(user_id);
CREATE INDEX IF NOT EXISTS idx_feedback_created ON feedback_responses(created_at);
