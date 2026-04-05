# S3 Integration Guide

## Overview

Nvwa-Lite now uses AWS S3 for storing:
- User-uploaded .h5ad dataset files
- Generated plot images (PNG)
- Generated table data (CSV)

## Configuration

### Environment Variables

Add to `.env`:
```bash
AWS_REGION=us-east-2
S3_BUCKET_NAME=nvwa-mvp-pilot
```

### IAM Permissions

The EC2 instance role needs:
- `s3:PutObject` - Upload files
- `s3:GetObject` - Download files
- `s3:DeleteObject` - Cleanup (optional)

## S3 Key Structure

```
users/{user_id}/sessions/{session_id}/
  ├── uploads/{filename}.h5ad
  └── results/
      ├── plot/plot_{message_id}_{idx}.png
      └── table/table_{message_id}_{idx}.csv
```

## Database Schema

The `message_artifacts` table includes an `s3_key` column:
```sql
ALTER TABLE message_artifacts ADD COLUMN s3_key TEXT;
```

## Migration

Run the migration to add S3 support:
```bash
bash scripts/migrate_db.sh
```

## Fallback Behavior

- If S3 is unavailable, files are stored locally in `data/uploads/`
- Plots/tables are always stored as base64 in the database as backup
- S3 failures are logged but don't block operations

## Testing

Local testing works without S3 - the system gracefully falls back to local storage.
