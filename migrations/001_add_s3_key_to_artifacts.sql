-- Add s3_key column to message_artifacts table for S3 storage integration
-- Migration: 001_add_s3_key_to_artifacts
-- Date: 2026-03-25

ALTER TABLE message_artifacts
ADD COLUMN IF NOT EXISTS s3_key TEXT;

COMMENT ON COLUMN message_artifacts.s3_key IS 'S3 object key for the artifact (plot PNG or table CSV)';
