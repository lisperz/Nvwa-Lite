"""Browser-to-S3 direct multipart upload component.

Key design constraint: st.file_uploader loses the file DOM reference on every
Streamlit rerun. Therefore, presigned URLs must be generated and passed to JS
in the SAME render that shows the file uploader — no rerun between file
selection and JS upload start.

State machine:
  idle      → user hasn't selected a file
  uploading → file selected, URLs generated, JS is uploading (no rerun yet)
  completing → JS finished, Python calls complete_multipart_upload
  completed → temp file ready, return result
"""

from __future__ import annotations

import json
import logging
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import streamlit as st
import streamlit.components.v1 as components

if TYPE_CHECKING:
    from src.storage.service import S3StorageService

logger = logging.getLogger(__name__)

_CHUNK_SIZE = 5 * 1024 * 1024  # 5 MB — S3 minimum part size
_COMPLETION_KEY = "s3_upload_completion"


def s3_direct_upload_widget(
    user_id: str,
    session_id: str,
    s3_service: S3StorageService,
    max_size_mb: int = 2000,
) -> tuple[str | None, str | None, Path | None]:
    """Browser-to-S3 direct multipart upload widget.

    Returns:
        (filename, s3_key, temp_path) on success, (None, None, None) otherwise.
    """
    if "s3_upload_state" not in st.session_state:
        _reset_state()

    state = st.session_state.s3_upload_state

    st.subheader("Upload Data")
    st.caption(f"Max file size: {max_size_mb} MB")

    # Phase: completing — JS finished, finalize on S3 and download
    if state["status"] == "completing":
        filename = state["filename"]
        s3_key = state["s3_key"]
        parts = state["parts"]

        with st.spinner("Finalizing and loading dataset..."):
            try:
                s3_service.complete_multipart_upload(s3_key, state["upload_id"], parts)
                logger.info(f"Multipart upload complete: {s3_key}")

                temp_path = Path(
                    tempfile.mktemp(suffix=".h5ad", prefix=f"nvwa_{session_id}_")
                )
                s3_service.download_file_to_path(s3_key, temp_path)
                logger.info(f"Downloaded to temp: {temp_path}")

                state["temp_path"] = str(temp_path)
                state["status"] = "completed"
                st.rerun()
            except Exception as e:
                st.error(f"Failed to finalize upload: {e}")
                logger.exception("Failed to complete multipart upload")
                _abort_and_reset(s3_service, state)
                return None, None, None

    # Phase: completed
    elif state["status"] == "completed":
        st.success(f"✓ {state['filename']} ready")
        return state["filename"], state["s3_key"], Path(state["temp_path"])

    # Phase: idle or uploading — always show file uploader so file reference stays alive
    else:
        uploaded_file = st.file_uploader(
            "Upload a .h5ad file",
            type=["h5ad"],
            key="s3_file_ref",
            help="File goes directly to cloud storage — server never handles file bytes.",
        )

        if uploaded_file is not None:
            filename = uploaded_file.name
            file_size = uploaded_file.size
            size_mb = file_size / (1024 * 1024)

            if size_mb > max_size_mb:
                st.error(f"File too large ({size_mb:.1f} MB). Maximum: {max_size_mb} MB.")
                return None, None, None

            # Generate presigned URLs if this is a new file
            if state["filename"] != filename or state["status"] == "idle":
                with st.spinner("Preparing upload..."):
                    try:
                        from src.storage.service import prepare_multipart_upload
                        upload_info = prepare_multipart_upload(
                            s3_service, user_id, session_id, filename,
                            file_size, _CHUNK_SIZE,
                        )
                        state["filename"] = filename
                        state["file_size"] = file_size
                        state["upload_id"] = upload_info["upload_id"]
                        state["s3_key"] = upload_info["s3_key"]
                        state["upload_info"] = upload_info
                        state["status"] = "uploading"
                        logger.info(
                            f"Prepared {len(upload_info['part_urls'])} part URLs for {filename}"
                        )
                    except Exception as e:
                        st.error(f"Failed to prepare upload: {e}")
                        logger.exception("Failed to prepare multipart upload")
                        _reset_state()
                        return None, None, None

            # Render progress bar and JS uploader
            # The file input stays visible so JS can access the file reference
            st.caption(f"📁 {filename} ({file_size / (1024*1024):.1f} MB)")

            # Hidden text input — JS sets this to signal completion
            completion_json = st.text_input(
                "upload_completion",
                key=_COMPLETION_KEY,
                label_visibility="collapsed",
            )

            # Check if JS signaled completion from a previous interaction
            if completion_json:
                try:
                    completion = json.loads(completion_json)
                    if completion.get("action") == "upload_complete":
                        state["parts"] = completion["parts"]
                        state["status"] = "completing"
                        st.rerun()
                    elif completion.get("action") == "upload_error":
                        err = completion.get("error", "Unknown error")
                        st.error(f"Upload failed: {err}")
                        _abort_and_reset(s3_service, state)
                        return None, None, None
                except json.JSONDecodeError:
                    pass

            # Render JS uploader — file input is still in the DOM above
            _render_js_uploader(state["upload_info"])

    return None, None, None


def _reset_state() -> None:
    st.session_state.s3_upload_state = {
        "status": "idle",
        "filename": None,
        "file_size": None,
        "upload_id": None,
        "s3_key": None,
        "upload_info": None,
        "parts": None,
        "temp_path": None,
    }


def _abort_and_reset(s3_service: S3StorageService, state: dict) -> None:
    try:
        if state.get("s3_key") and state.get("upload_id"):
            s3_service.abort_multipart_upload(state["s3_key"], state["upload_id"])
    except Exception:
        pass
    _reset_state()


def _render_js_uploader(upload_info: dict) -> None:
    """Render JS multipart uploader. On completion, sets the hidden text input."""
    upload_info_json = json.dumps(upload_info)
    html = f"""
    <style>
      body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; padding: 8px; }}
      .progress-wrap {{ background: #e0e0e0; border-radius: 10px; height: 26px; overflow: hidden; margin: 6px 0; }}
      .progress-bar {{
        height: 100%; background: linear-gradient(90deg, #4CAF50, #45a049);
        width: 0%; transition: width 0.2s; display: flex; align-items: center;
        justify-content: center; color: white; font-size: 13px; font-weight: bold;
        min-width: 36px;
      }}
      .status {{ font-size: 13px; color: #555; margin-top: 4px; }}
    </style>
    <div class="progress-wrap"><div class="progress-bar" id="bar">0%</div></div>
    <div class="status" id="status">Starting upload...</div>

    <script>
    (async () => {{
      const info = {upload_info_json};
      const {{ upload_id, s3_key, chunk_size, part_urls }} = info;
      const bar = document.getElementById('bar');
      const status = document.getElementById('status');

      // Find the .h5ad file from Streamlit's file uploader in the parent frame
      // It is still in the DOM because we keep the file uploader rendered above
      function findFile() {{
        const inputs = window.parent.document.querySelectorAll('input[type=file]');
        for (const inp of inputs) {{
          if (inp.files && inp.files.length > 0 && inp.files[0].name.endsWith('.h5ad')) {{
            return inp.files[0];
          }}
        }}
        return null;
      }}

      // Wait up to 3 seconds for the file input to be available
      let file = null;
      for (let i = 0; i < 30; i++) {{
        file = findFile();
        if (file) break;
        await new Promise(r => setTimeout(r, 100));
      }}

      if (!file) {{
        status.textContent = '⚠ File reference not found. Please re-select the file.';
        return;
      }}

      const totalBytes = file.size;
      let uploadedBytes = 0;
      const results = [];
      let uploadError = null;

      function updateProgress(added) {{
        uploadedBytes += added;
        const pct = Math.min(99, Math.round((uploadedBytes / totalBytes) * 100));
        bar.style.width = pct + '%';
        bar.textContent = pct + '%';
        const mb = (uploadedBytes / 1024 / 1024).toFixed(1);
        const total = (totalBytes / 1024 / 1024).toFixed(1);
        status.textContent = `Uploading ${{mb}} / ${{total}} MB (${{pct}}%)`;
      }}

      async function uploadPart(partInfo) {{
        const {{ part_number, url }} = partInfo;
        const start = (part_number - 1) * chunk_size;
        const end = Math.min(start + chunk_size, file.size);
        const chunk = file.slice(start, end);

        return new Promise((resolve, reject) => {{
          const xhr = new XMLHttpRequest();
          let lastLoaded = 0;
          xhr.upload.onprogress = e => {{
            updateProgress(e.loaded - lastLoaded);
            lastLoaded = e.loaded;
          }};
          xhr.onload = () => {{
            if (xhr.status >= 200 && xhr.status < 300) {{
              resolve({{ PartNumber: part_number, ETag: xhr.getResponseHeader('ETag') }});
            }} else {{
              reject(new Error(`Part ${{part_number}} failed: HTTP ${{xhr.status}}`));
            }}
          }};
          xhr.onerror = () => reject(new Error(`Part ${{part_number}} network error`));
          xhr.open('PUT', url);
          xhr.setRequestHeader('Content-Type', 'application/octet-stream');
          xhr.send(chunk);
        }});
      }}

      // Upload with concurrency of 3
      const queue = [...part_urls];
      const active = new Set();

      async function runNext() {{
        if (queue.length === 0 || uploadError) return;
        const partInfo = queue.shift();
        active.add(partInfo.part_number);
        try {{
          results.push(await uploadPart(partInfo));
        }} catch (e) {{
          uploadError = e.message;
        }} finally {{
          active.delete(partInfo.part_number);
          await runNext();
        }}
      }}

      await Promise.all(
        Array.from({{ length: Math.min(3, part_urls.length) }}, () => runNext())
      );
      while (active.size > 0) await new Promise(r => setTimeout(r, 100));

      if (uploadError) {{
        bar.style.background = '#f44336';
        status.textContent = '✗ ' + uploadError;
        signalStreamlit(JSON.stringify({{ action: 'upload_error', error: uploadError }}));
        return;
      }}

      results.sort((a, b) => a.PartNumber - b.PartNumber);
      bar.style.width = '100%';
      bar.textContent = '100%';
      status.textContent = '✓ Upload complete — finalizing...';
      signalStreamlit(JSON.stringify({{ action: 'upload_complete', s3_key, parts: results }}));
    }})();

    function signalStreamlit(value) {{
      // Set the hidden text input in the parent Streamlit frame
      const inputs = window.parent.document.querySelectorAll('input[type=text]');
      for (const inp of inputs) {{
        // The hidden input has an empty value — find it by checking for empty
        if (inp.value === '') {{
          const setter = Object.getOwnPropertyDescriptor(
            window.parent.HTMLInputElement.prototype, 'value'
          ).set;
          setter.call(inp, value);
          inp.dispatchEvent(new Event('input', {{ bubbles: true }}));
          return;
        }}
      }}
    }}
    </script>
    """
    components.html(html, height=80, scrolling=False)
