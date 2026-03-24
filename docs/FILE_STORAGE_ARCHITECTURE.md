# File Upload Storage Architecture

## Current Implementation

### Storage Location
**Local filesystem storage** - Files are currently stored locally, NOT in S3.

**Upload Directory:**
```
/Users/zhuchen/Downloads/Nvwa Bio Technical Challange/nvwa-lite/data/uploads/
```

**Docker Volume Mount:**
```yaml
volumes:
  - ./data/uploads:/app/data/uploads
```

This means uploaded files are stored on the host machine at `data/uploads/` and mounted into the Docker container at `/app/data/uploads/`.

### Current Files in Upload Directory
```
data/uploads/
├── pbmc_test.h5ad (24.7 MB)
├── Sub UMap of the forebrain.h5ad (42.6 MB)
├── Test.h5ad (71.7 MB)
└── Turtle.h5ad (65.0 MB)

Total: ~204 MB
```

### How Upload Works

**File Upload Flow:**
1. User uploads file via Streamlit file uploader (`src/ui/components.py:18-42`)
2. File is read into memory as bytes
3. File size is validated (max 2000 MB)
4. File is written to local disk: `data/uploads/{filename}`
5. Path is returned to the app for processing

**Code Location:**
```python
# src/ui/components.py
UPLOAD_DIR = Path("data/uploads")
MAX_UPLOAD_MB = 2000

def file_upload_widget() -> Path | None:
    uploaded = st.file_uploader("Upload a .h5ad file", type=["h5ad"])
    if uploaded is not None:
        file_bytes = uploaded.getvalue()
        UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
        dest = UPLOAD_DIR / uploaded.name
        dest.write_bytes(file_bytes)
        return dest
```

## S3 Storage Service (Available but Not Used)

### Implementation Status
✓ **Code exists** - `src/storage/service.py` contains full S3 implementation
✗ **Not configured** - No S3 bucket name in `.env`
✗ **Not integrated** - App uses local storage, not S3 service

### S3 Service Features (Ready to Use)
The `S3StorageService` class provides:
- Pre-signed URLs for direct browser uploads
- Per-user namespace isolation: `users/{user_id}/sessions/{session_id}/uploads/{filename}`
- Automatic TTL-based cleanup (14 days default)
- Upload/download methods
- Result storage for plots, tables, logs

### Fallback: LocalStorageService
The code also includes `LocalStorageService` class that mimics S3 API but stores files locally at `data/storage/` (different from current `data/uploads/`).

## Migration Path to S3

### Prerequisites
1. Create S3 bucket in AWS
2. Configure AWS credentials
3. Set environment variables

### Required Environment Variables
```bash
S3_BUCKET_NAME=your-bucket-name
AWS_ACCESS_KEY_ID=your-access-key
AWS_SECRET_ACCESS_KEY=your-secret-key
AWS_REGION=us-west-2  # or your preferred region
```

### Code Changes Needed
1. **Update `src/ui/components.py`:**
   - Replace direct file write with S3StorageService
   - Use pre-signed URLs or direct upload

2. **Update `src/ui/app.py`:**
   - Initialize S3StorageService
   - Update file loading to download from S3
   - Update session tracking to use S3 keys

3. **Update Docker configuration:**
   - Remove `./data/uploads` volume mount (no longer needed)
   - Add AWS credentials to environment

### Benefits of S3 Migration
- **Scalability:** No local disk space limits
- **Durability:** AWS handles backups and redundancy
- **Multi-instance:** Multiple containers can share same storage
- **Automatic cleanup:** TTL-based deletion of old files
- **Security:** Pre-signed URLs with expiration
- **Cost:** Pay only for what you use

### Risks of Current Local Storage
- **Data loss:** Files deleted if container is removed
- **Disk space:** Limited by host machine
- **Single instance:** Can't scale horizontally
- **No backup:** Files not backed up automatically
- **Manual cleanup:** Old files accumulate

## Recommendations

### For Development/Testing
✓ Current local storage is fine
- Simple and fast
- No AWS costs
- Easy to debug

### For Production/Pilot
✗ Must migrate to S3
- Required for scalability
- Required for multi-instance deployment
- Required for data durability
- Required for automatic cleanup

### Migration Priority
**Medium-High Priority**
- Not blocking for single-user testing
- Critical before multi-user pilot
- Critical before production deployment

## Next Steps

1. **Create S3 bucket** in AWS console
2. **Configure IAM permissions** for bucket access
3. **Add environment variables** to `.env`
4. **Test S3StorageService** with existing code
5. **Update UI components** to use S3
6. **Migrate existing files** from `data/uploads/` to S3
7. **Remove local volume mount** from docker-compose.yml
8. **Test end-to-end** upload/download flow
9. **Set lifecycle policy** for automatic cleanup

## Cost Estimate (AWS S3)

**Assumptions:**
- 10 users
- 100 MB average file size per user
- 5 sessions per user per month
- 14-day retention

**Monthly Storage:**
- Total data: 10 users × 100 MB × 5 sessions = 5 GB
- S3 Standard: $0.023/GB = ~$0.12/month

**Monthly Requests:**
- Uploads: 50 PUT requests = $0.005/1000 = ~$0.00025
- Downloads: 200 GET requests = $0.0004/1000 = ~$0.00008

**Total: ~$0.12/month** (negligible cost)
