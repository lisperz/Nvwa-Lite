# Large File Upload Architecture Proposal

## Current Problem
- Cloudflare has 100 MB upload limit (Free/Pro plans)
- Nginx/Streamlit proxy adds latency and consumes server bandwidth
- Not scalable for many concurrent users uploading large files

## Proposed Solution: Direct S3 Uploads with Presigned URLs

### Architecture Flow

```
┌─────────────┐
│   Browser   │
│  (Streamlit)│
└──────┬──────┘
       │ 1. Request upload URL
       ▼
┌─────────────────┐
│  Streamlit App  │
│   (Backend)     │
└──────┬──────────┘
       │ 2. Generate presigned URL
       │    (boto3.generate_presigned_post)
       ▼
┌─────────────┐
│   Browser   │
└──────┬──────┘
       │ 3. Upload directly to S3
       │    (bypasses nginx/Cloudflare)
       ▼
┌─────────────┐
│   AWS S3    │
│   Bucket    │
└──────┬──────┘
       │ 4. Upload complete
       ▼
┌─────────────────┐
│  Streamlit App  │
│  (Load from S3) │
└─────────────────┘
```

### Implementation Steps

#### 1. Backend: Generate Presigned URL

```python
# src/storage/presigned.py
import boto3
from datetime import timedelta

def generate_upload_url(
    user_id: str,
    session_id: str,
    filename: str,
    file_size: int,
    expiration: int = 3600
) -> dict:
    """Generate presigned POST URL for direct S3 upload.

    Returns:
        {
            'url': 'https://s3.amazonaws.com/bucket',
            'fields': {'key': '...', 'policy': '...', ...}
        }
    """
    s3_client = boto3.client('s3', region_name='us-east-2')

    s3_key = f"users/{user_id}/sessions/{session_id}/uploads/{filename}"

    # Generate presigned POST (allows browser to upload directly)
    presigned_post = s3_client.generate_presigned_post(
        Bucket='nvwa-mvp-pilot',
        Key=s3_key,
        Fields={
            'Content-Type': 'application/octet-stream',
            'x-amz-meta-user-id': user_id,
            'x-amz-meta-session-id': session_id,
        },
        Conditions=[
            ['content-length-range', 1, 5368709120],  # 1 byte to 5 GB
        ],
        ExpiresIn=expiration
    )

    return {
        'url': presigned_post['url'],
        'fields': presigned_post['fields'],
        's3_key': s3_key,
        'expires_in': expiration
    }
```

#### 2. Frontend: Upload with JavaScript

```python
# src/ui/components.py - Modified file_upload_widget

import streamlit as st
import streamlit.components.v1 as components

def file_upload_widget_v2(user_id: str, session_id: str):
    """File upload widget using direct S3 upload."""

    # Custom HTML/JS component for direct S3 upload
    upload_html = """
    <div id="upload-container">
        <input type="file" id="file-input" accept=".h5ad" />
        <button id="upload-btn">Upload to S3</button>
        <div id="progress"></div>
    </div>

    <script>
    document.getElementById('upload-btn').onclick = async () => {
        const file = document.getElementById('file-input').files[0];
        if (!file) return;

        // 1. Request presigned URL from backend
        const response = await fetch('/api/upload-url', {
            method: 'POST',
            body: JSON.stringify({
                filename: file.name,
                filesize: file.size
            })
        });
        const {url, fields, s3_key} = await response.json();

        // 2. Upload directly to S3
        const formData = new FormData();
        Object.entries(fields).forEach(([k, v]) => formData.append(k, v));
        formData.append('file', file);

        const xhr = new XMLHttpRequest();
        xhr.upload.onprogress = (e) => {
            const percent = (e.loaded / e.total * 100).toFixed(1);
            document.getElementById('progress').innerText = `${percent}%`;
        };

        xhr.onload = () => {
            if (xhr.status === 204) {
                // 3. Notify Streamlit that upload is complete
                window.parent.postMessage({
                    type: 'upload_complete',
                    s3_key: s3_key
                }, '*');
            }
        };

        xhr.open('POST', url);
        xhr.send(formData);
    };
    </script>
    """

    components.html(upload_html, height=200)
```

#### 3. Backend: Load from S3

```python
# src/ui/app.py - Modified to load from S3

def load_dataset_from_s3(s3_key: str) -> ad.AnnData:
    """Download .h5ad from S3 and load into memory."""
    import boto3
    import tempfile

    s3_client = boto3.client('s3', region_name='us-east-2')

    with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp:
        s3_client.download_file('nvwa-mvp-pilot', s3_key, tmp.name)
        adata = ad.read_h5ad(tmp.name)

    return adata
```

### Benefits Summary

| Aspect | Current (Proxy Upload) | Proposed (Direct S3) |
|--------|----------------------|---------------------|
| Max file size | 100 MB (Cloudflare) | 5 TB (S3 limit) |
| Upload speed | Slow (proxied) | Fast (direct) |
| Server bandwidth | High | Zero |
| Cloudflare cost | $200/mo for 500 MB | $0 (bypassed) |
| Scalability | Limited by EC2 | Unlimited (S3) |
| Resumable | No | Yes (multipart) |

### Migration Plan

1. **Phase 1**: Implement presigned URL generation (backend)
2. **Phase 2**: Create custom upload component (frontend)
3. **Phase 3**: Test with pilot users
4. **Phase 4**: Switch all users to new upload flow
5. **Phase 5**: Remove old upload code

### Additional Enhancements

#### A. Multipart Upload for Files > 100 MB
- Split large files into 5 MB chunks
- Upload chunks in parallel
- Resume failed uploads

#### B. Upload Progress Tracking
- Store upload state in Redis
- Show real-time progress to user
- Handle browser refresh gracefully

#### C. Background Processing
- After upload completes, process file asynchronously
- Send notification when ready
- Allows user to continue browsing

#### D. File Validation
- Validate .h5ad format before processing
- Check file integrity (checksums)
- Reject corrupted files early

## Cost Analysis

### Current Approach (Cloudflare Business)
- Cloudflare Business: $200/month
- 500 MB upload limit
- Still consumes EC2 bandwidth

### Proposed Approach (Direct S3)
- Cloudflare: $0 (stay on Free plan)
- S3 storage: ~$0.023/GB/month
- S3 upload: Free (ingress)
- S3 download: $0.09/GB (to EC2 in same region)

**Example**: 100 users × 500 MB/month = 50 GB
- Storage: 50 GB × $0.023 = $1.15/month
- Download: 50 GB × $0.09 = $4.50/month
- **Total: ~$6/month** vs $200/month for Cloudflare Business

## Timeline Estimate

- **Week 1**: Backend presigned URL generation + S3 download
- **Week 2**: Frontend custom upload component
- **Week 3**: Testing + bug fixes
- **Week 4**: Production deployment + monitoring

## Risks & Mitigations

| Risk | Mitigation |
|------|-----------|
| Browser compatibility | Test on Chrome, Firefox, Safari |
| CORS issues | Configure S3 bucket CORS policy |
| Upload failures | Implement retry logic + multipart |
| Security | Presigned URLs expire after 1 hour |

## Conclusion

Direct S3 uploads are the industry standard for handling large files at scale. This approach will:
- Eliminate all current upload size limitations
- Reduce infrastructure costs significantly
- Scale to thousands of concurrent users
- Improve upload speed and reliability

**Recommendation**: Implement this architecture before scaling to more users.
