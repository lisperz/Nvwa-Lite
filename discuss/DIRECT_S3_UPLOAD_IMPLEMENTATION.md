# Direct S3 Upload Implementation Plan

## Current State Analysis

### What You Already Have ✅
1. **`S3StorageService.generate_upload_url()`** - Creates presigned URLs (storage/service.py:48-80)
2. **S3 bucket configured** - `nvwa-mvp-pilot` with IAM role
3. **User/session management** - Already tracking user_id and session_id

### What's Wrong ❌
The `file_upload_widget()` in `components.py` uses Streamlit's `st.file_uploader()` which:
- Receives file through Cloudflare → nginx → Streamlit (proxy upload)
- Hits Cloudflare's 100 MB limit
- Consumes server bandwidth
- Then uploads to S3 from backend (double upload!)

## Recommended Fix: Minimal Changes Approach

### Option 1: Custom HTML Upload Component (Recommended)

Replace `st.file_uploader()` with a custom HTML/JS component that:
1. Requests presigned URL from backend
2. Uploads directly to S3 from browser
3. Notifies Streamlit when complete

**Pros:**
- Bypasses Cloudflare completely
- No size limits
- Fast uploads
- Uses your existing `generate_upload_url()` method

**Cons:**
- Requires custom JavaScript component
- More complex than Streamlit widget

### Option 2: Turn Off Cloudflare + Keep Current Code (Quick Fix)

Turn off Cloudflare proxy for `nvwa.bio`:
- URL stays the same
- Current code works as-is
- Uploads up to 2 GB work

**Pros:**
- Zero code changes
- Immediate fix
- Simple

**Cons:**
- Still uses proxy upload (slower)
- Consumes server bandwidth
- Loses DDoS protection

## Implementation: Option 1 (Direct S3 Upload)

### Step 1: Create Custom Upload Component

```python
# src/ui/s3_uploader.py
import streamlit as st
import streamlit.components.v1 as components
from src.storage.service import S3StorageService
import os

def s3_file_uploader(user_id: str, session_id: str) -> str | None:
    """Custom file uploader that uploads directly to S3.

    Returns:
        S3 key of uploaded file, or None if no upload.
    """

    # Generate presigned URL
    if 'upload_url' not in st.session_state:
        s3_service = S3StorageService(
            bucket_name=os.getenv("S3_BUCKET_NAME"),
            region=os.getenv("AWS_REGION", "us-east-2")
        )
        # This will be called when user selects a file
        st.session_state.upload_url = None
        st.session_state.s3_key = None

    # Custom HTML/JS component
    upload_html = f"""
    <div style="padding: 20px; border: 2px dashed #ccc; border-radius: 8px;">
        <input type="file" id="file-input" accept=".h5ad"
               style="margin-bottom: 10px;" />
        <button id="upload-btn"
                style="padding: 8px 16px; background: #ff4b4b; color: white;
                       border: none; border-radius: 4px; cursor: pointer;">
            Upload to S3
        </button>
        <div id="progress" style="margin-top: 10px;"></div>
        <div id="status" style="margin-top: 10px; color: green;"></div>
    </div>

    <script>
    const uploadBtn = document.getElementById('upload-btn');
    const fileInput = document.getElementById('file-input');
    const progress = document.getElementById('progress');
    const status = document.getElementById('status');

    uploadBtn.onclick = async () => {{
        const file = fileInput.files[0];
        if (!file) {{
            alert('Please select a file first');
            return;
        }}

        uploadBtn.disabled = true;
        uploadBtn.innerText = 'Uploading...';

        try {{
            // 1. Request presigned URL from Streamlit backend
            const response = await fetch(window.location.origin + '/api/upload-url', {{
                method: 'POST',
                headers: {{'Content-Type': 'application/json'}},
                body: JSON.stringify({{
                    user_id: '{user_id}',
                    session_id: '{session_id}',
                    filename: file.name
                }})
            }});

            const {{url, s3_key}} = await response.json();

            // 2. Upload directly to S3 using presigned URL
            const xhr = new XMLHttpRequest();

            xhr.upload.onprogress = (e) => {{
                if (e.lengthComputable) {{
                    const percent = (e.loaded / e.total * 100).toFixed(1);
                    progress.innerHTML = `<div style="width: ${{percent}}%; height: 20px;
                                          background: #ff4b4b; border-radius: 4px;"></div>
                                          <div>${{percent}}% uploaded</div>`;
                }}
            }};

            xhr.onload = () => {{
                if (xhr.status === 200) {{
                    status.innerText = '✅ Upload complete!';
                    // Notify Streamlit
                    window.parent.postMessage({{
                        type: 'streamlit:setComponentValue',
                        value: s3_key
                    }}, '*');
                }} else {{
                    status.innerText = '❌ Upload failed: ' + xhr.statusText;
                    status.style.color = 'red';
                }}
                uploadBtn.disabled = false;
                uploadBtn.innerText = 'Upload to S3';
            }};

            xhr.onerror = () => {{
                status.innerText = '❌ Upload failed';
                status.style.color = 'red';
                uploadBtn.disabled = false;
                uploadBtn.innerText = 'Upload to S3';
            }};

            xhr.open('PUT', url);
            xhr.setRequestHeader('Content-Type', 'application/octet-stream');
            xhr.send(file);

        }} catch (error) {{
            status.innerText = '❌ Error: ' + error.message;
            status.style.color = 'red';
            uploadBtn.disabled = false;
            uploadBtn.innerText = 'Upload to S3';
        }}
    }};
    </script>
    """

    # Render component and get S3 key when upload completes
    s3_key = components.html(upload_html, height=200)

    return s3_key
```

### Step 2: Add API Endpoint for Presigned URLs

Since Streamlit doesn't have built-in API routes, we need to add a simple endpoint:

```python
# src/api/upload.py
from fastapi import FastAPI, Request
from src.storage.service import S3StorageService
import os

app = FastAPI()

@app.post("/api/upload-url")
async def get_upload_url(request: Request):
    """Generate presigned URL for direct S3 upload."""
    data = await request.json()

    s3_service = S3StorageService(
        bucket_name=os.getenv("S3_BUCKET_NAME"),
        region=os.getenv("AWS_REGION", "us-east-2")
    )

    url, s3_key = s3_service.generate_upload_url(
        user_id=data['user_id'],
        session_id=data['session_id'],
        filename=data['filename']
    )

    return {"url": url, "s3_key": s3_key}
```

### Step 3: Update app.py to Use New Uploader

```python
# src/ui/app.py - Replace file_upload_widget with s3_file_uploader

from src.ui.s3_uploader import s3_file_uploader

# In main():
s3_key = s3_file_uploader(user_id, session_id)

if s3_key:
    # Download from S3 and load
    s3_service = S3StorageService(...)
    file_bytes = s3_service.download_file(s3_key)

    # Save locally for processing
    dest = UPLOAD_DIR / Path(s3_key).name
    dest.write_bytes(file_bytes)

    # Load into AnnData
    adata = ad.read_h5ad(dest)
```

## Simpler Alternative: Streamlit + S3 Presigned POST

Actually, there's an even simpler approach using Streamlit's existing uploader but with presigned POST:

```python
# src/ui/components.py - Modified version

def file_upload_widget_v2(user_id: str, session_id: str):
    """File upload with S3 presigned POST."""

    with st.sidebar:
        st.subheader("Upload Data")

        # Show upload instructions
        if st.button("Get Upload Link"):
            s3_service = S3StorageService(...)
            url, s3_key = s3_service.generate_upload_url(
                user_id, session_id, "dataset.h5ad"
            )

            st.info(f"""
            **Direct S3 Upload Link (bypasses Cloudflare):**

            Use this curl command to upload your file:
            ```bash
            curl -X PUT "{url}" \\
                 -H "Content-Type: application/octet-stream" \\
                 --data-binary @your_file.h5ad
            ```

            Or use this Python script:
            ```python
            import requests
            with open('your_file.h5ad', 'rb') as f:
                requests.put('{url}', data=f)
            ```

            After upload, click "Load from S3" below.
            """)

            st.session_state.pending_s3_key = s3_key

        # Load from S3
        if st.button("Load from S3") and 'pending_s3_key' in st.session_state:
            s3_key = st.session_state.pending_s3_key
            s3_service = S3StorageService(...)
            file_bytes = s3_service.download_file(s3_key)

            # Process file...
            return file_bytes, s3_key

    return None, None
```

## Recommendation

**For immediate fix (today):**
- Turn off Cloudflare proxy for `nvwa.bio`
- No code changes needed
- Uploads up to 2 GB work immediately

**For long-term scalability (next sprint):**
- Implement custom S3 upload component (Option 1)
- Or use the simpler "Get Upload Link" approach above
- This gives you unlimited file size and scales infinitely

## Cost Comparison

| Approach | Setup Time | File Size Limit | Monthly Cost | Scalability |
|----------|-----------|-----------------|--------------|-------------|
| Current (broken) | 0 | 100 MB | $0 | Poor |
| Turn off Cloudflare | 5 min | 2 GB | $0 | Medium |
| Direct S3 upload | 1-2 weeks | 5 TB | ~$6 | Unlimited |
| Cloudflare Business | 5 min | 500 MB | $200 | Medium |

## Next Steps

1. **Immediate**: Turn off Cloudflare proxy (5 minutes)
2. **This week**: Test with pilot users
3. **Next sprint**: Implement direct S3 upload component
4. **Before scaling**: Deploy direct S3 upload to production
