"""S3 storage service for dataset uploads and result storage.

Manages file uploads/downloads using S3 with pre-signed URLs.
Implements per-user namespace isolation and TTL-based cleanup.
"""

from __future__ import annotations

import os
from pathlib import Path

try:
    import boto3
    from botocore.config import Config
    from botocore.exceptions import ClientError
    BOTO3_AVAILABLE = True
except ImportError:
    BOTO3_AVAILABLE = False


class S3StorageService:
    """Manage dataset uploads and result storage in S3.

    For the pilot, this provides:
    - Pre-signed URLs for direct browser uploads (avoids server bottleneck)
    - Per-user namespace isolation (users/<user_id>/sessions/<session_id>/...)
    - Automatic TTL-based cleanup (14 days default)
    """

    def __init__(self, bucket_name: str | None = None, region: str = "us-west-2"):
        """Initialize the S3 storage service.

        Args:
            bucket_name: S3 bucket name. If None, reads from environment.
            region: AWS region.
        """
        if not BOTO3_AVAILABLE:
            raise ImportError("boto3 is required for S3 storage. Install with: uv add boto3")

        self.bucket_name = bucket_name or os.getenv("S3_BUCKET_NAME")
        if not self.bucket_name:
            raise ValueError("S3_BUCKET_NAME must be set in environment or passed to constructor")

        self.region = region
        self.s3_client = boto3.client(
            "s3",
            region_name=region,
            config=Config(s3={"addressing_style": "virtual"}),
        )

    def generate_upload_url(
        self,
        user_id: str,
        session_id: str,
        filename: str,
        expires_in: int = 3600,
    ) -> tuple[str, str]:
        """Generate pre-signed URL for dataset upload.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            filename: Original filename.
            expires_in: URL expiration time in seconds (default 1 hour).

        Returns:
            Tuple of (presigned_url, s3_key).
        """
        key = f"users/{user_id}/sessions/{session_id}/uploads/{filename}"

        try:
            url = self.s3_client.generate_presigned_url(
                "put_object",
                Params={
                    "Bucket": self.bucket_name,
                    "Key": key,
                    "ContentType": "application/octet-stream",
                },
                ExpiresIn=expires_in,
            )
            return url, key
        except ClientError as e:
            raise RuntimeError(f"Failed to generate upload URL: {e}")

    def upload_file(
        self,
        user_id: str,
        session_id: str,
        file_data: bytes,
        filename: str,
        file_type: str = "upload",
    ) -> str:
        """Upload a file directly to S3.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            file_data: File content as bytes.
            filename: Filename.
            file_type: Type of file ("upload", "result", "log").

        Returns:
            S3 key of the uploaded file.
        """
        key = f"users/{user_id}/sessions/{session_id}/{file_type}s/{filename}"

        try:
            self.s3_client.put_object(
                Bucket=self.bucket_name,
                Key=key,
                Body=file_data,
                ContentType=self._get_content_type(filename),
            )
            return key
        except ClientError as e:
            raise RuntimeError(f"Failed to upload file: {e}")

    def upload_result(
        self,
        user_id: str,
        session_id: str,
        result_type: str,
        data: bytes,
        filename: str,
    ) -> str:
        """Upload analysis result to S3.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            result_type: Type of result ("plot", "table", "log").
            data: Result data as bytes.
            filename: Filename.

        Returns:
            S3 key of the uploaded result.
        """
        key = f"users/{user_id}/sessions/{session_id}/results/{result_type}/{filename}"

        try:
            self.s3_client.put_object(
                Bucket=self.bucket_name,
                Key=key,
                Body=data,
                ContentType=self._get_content_type(filename),
            )
            return key
        except ClientError as e:
            raise RuntimeError(f"Failed to upload result: {e}")

    def generate_download_url(
        self,
        key: str,
        expires_in: int = 3600,
    ) -> str:
        """Generate pre-signed URL for file download.

        Args:
            key: S3 object key.
            expires_in: URL expiration time in seconds (default 1 hour).

        Returns:
            Pre-signed download URL.
        """
        try:
            url = self.s3_client.generate_presigned_url(
                "get_object",
                Params={"Bucket": self.bucket_name, "Key": key},
                ExpiresIn=expires_in,
            )
            return url
        except ClientError as e:
            raise RuntimeError(f"Failed to generate download URL: {e}")

    def download_file(self, key: str) -> bytes:
        """Download a file from S3.

        Args:
            key: S3 object key.

        Returns:
            File content as bytes.
        """
        try:
            response = self.s3_client.get_object(Bucket=self.bucket_name, Key=key)
            return response["Body"].read()
        except ClientError as e:
            raise RuntimeError(f"Failed to download file: {e}")

    def download_file_to_path(self, key: str, dest: Path) -> None:
        """Download a file from S3 directly to a local path (streaming, low memory).

        Args:
            key: S3 object key.
            dest: Local destination path.
        """
        try:
            self.s3_client.download_file(self.bucket_name, key, str(dest))
        except ClientError as e:
            raise RuntimeError(f"Failed to download file to path: {e}")

    def initiate_multipart_upload(
        self,
        user_id: str,
        session_id: str,
        filename: str,
    ) -> tuple[str, str]:
        """Start a multipart upload session.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            filename: Original filename.

        Returns:
            Tuple of (upload_id, s3_key).
        """
        key = f"users/{user_id}/sessions/{session_id}/uploads/{filename}"
        try:
            response = self.s3_client.create_multipart_upload(
                Bucket=self.bucket_name,
                Key=key,
                ContentType="application/octet-stream",
            )
            return response["UploadId"], key
        except ClientError as e:
            raise RuntimeError(f"Failed to initiate multipart upload: {e}")

    def generate_part_upload_url(
        self,
        s3_key: str,
        upload_id: str,
        part_number: int,
        expires_in: int = 7200,
    ) -> str:
        """Generate a presigned URL for uploading one part.

        Args:
            s3_key: S3 object key.
            upload_id: Multipart upload ID.
            part_number: Part number (1-indexed).
            expires_in: URL expiration in seconds (default 2 hours).

        Returns:
            Presigned URL string.
        """
        try:
            return self.s3_client.generate_presigned_url(
                "upload_part",
                Params={
                    "Bucket": self.bucket_name,
                    "Key": s3_key,
                    "UploadId": upload_id,
                    "PartNumber": part_number,
                },
                ExpiresIn=expires_in,
            )
        except ClientError as e:
            raise RuntimeError(f"Failed to generate part upload URL: {e}")

    def complete_multipart_upload(
        self,
        s3_key: str,
        upload_id: str,
        parts: list[dict],
    ) -> None:
        """Complete a multipart upload.

        Args:
            s3_key: S3 object key.
            upload_id: Multipart upload ID.
            parts: List of {"PartNumber": int, "ETag": str} dicts.
        """
        try:
            self.s3_client.complete_multipart_upload(
                Bucket=self.bucket_name,
                Key=s3_key,
                UploadId=upload_id,
                MultipartUpload={"Parts": parts},
            )
        except ClientError as e:
            raise RuntimeError(f"Failed to complete multipart upload: {e}")

    def abort_multipart_upload(self, s3_key: str, upload_id: str) -> None:
        """Abort and clean up a multipart upload.

        Args:
            s3_key: S3 object key.
            upload_id: Multipart upload ID.
        """
        try:
            self.s3_client.abort_multipart_upload(
                Bucket=self.bucket_name,
                Key=s3_key,
                UploadId=upload_id,
            )
        except ClientError as e:
            raise RuntimeError(f"Failed to abort multipart upload: {e}")

    def configure_cors(self, allowed_origins: list[str]) -> None:
        """Set CORS policy on the bucket for browser-direct uploads.

        Args:
            allowed_origins: List of allowed origins (e.g. ["https://nvwa.bio"]).
        """
        try:
            self.s3_client.put_bucket_cors(
                Bucket=self.bucket_name,
                CORSConfiguration={
                    "CORSRules": [{
                        "AllowedOrigins": allowed_origins,
                        "AllowedMethods": ["PUT"],
                        "AllowedHeaders": ["Content-Type", "Content-Length"],
                        "ExposeHeaders": ["ETag"],
                        "MaxAgeSeconds": 3600,
                    }]
                },
            )
        except ClientError as e:
            raise RuntimeError(f"Failed to configure CORS: {e}")

    def set_lifecycle_policy(self, days: int = 14) -> None:
        """Set TTL for automatic deletion of old data.

        Args:
            days: Number of days before deletion (default 14).
        """
        lifecycle_config = {
            "Rules": [
                {
                    "Id": "DeleteOldPilotData",
                    "Status": "Enabled",
                    "Prefix": "users/",
                    "Expiration": {"Days": days},
                }
            ]
        }

        try:
            self.s3_client.put_bucket_lifecycle_configuration(
                Bucket=self.bucket_name,
                LifecycleConfiguration=lifecycle_config,
            )
        except ClientError as e:
            raise RuntimeError(f"Failed to set lifecycle policy: {e}")

    def _get_content_type(self, filename: str) -> str:
        """Determine content type from filename.

        Args:
            filename: Filename with extension.

        Returns:
            MIME type string.
        """
        ext = Path(filename).suffix.lower()
        content_types = {
            ".h5ad": "application/octet-stream",
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".csv": "text/csv",
            ".json": "application/json",
            ".log": "text/plain",
        }
        return content_types.get(ext, "application/octet-stream")


def prepare_multipart_upload(
    s3_service: S3StorageService,
    user_id: str,
    session_id: str,
    filename: str,
    file_size: int,
    chunk_size: int = 5 * 1024 * 1024,
) -> dict:
    """Pre-generate all presigned part URLs for a multipart upload.

    For a 1GB file with 5MB chunks, this generates ~200 URLs (~100KB payload).
    All URLs are generated upfront and returned to the JS component in one round-trip.

    Args:
        s3_service: Initialized S3StorageService instance.
        user_id: User identifier.
        session_id: Session identifier.
        filename: Original filename.
        file_size: File size in bytes.
        chunk_size: Size of each chunk in bytes (min 5MB for S3 multipart).

    Returns:
        Dict with upload_id, s3_key, chunk_size, and list of {part_number, url}.
    """
    import math

    upload_id, s3_key = s3_service.initiate_multipart_upload(user_id, session_id, filename)
    num_parts = math.ceil(file_size / chunk_size)

    part_urls = [
        {
            "part_number": i + 1,
            "url": s3_service.generate_part_upload_url(s3_key, upload_id, i + 1),
        }
        for i in range(num_parts)
    ]

    return {
        "upload_id": upload_id,
        "s3_key": s3_key,
        "chunk_size": chunk_size,
        "part_urls": part_urls,
    }


class LocalStorageService:
    """Fallback local storage for development/testing without S3.

    Mimics S3 API but stores files locally.
    """

    def __init__(self, base_dir: Path = Path("data/storage")):
        """Initialize local storage.

        Args:
            base_dir: Base directory for file storage.
        """
        self.base_dir = base_dir
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def upload_file(
        self,
        user_id: str,
        session_id: str,
        file_data: bytes,
        filename: str,
        file_type: str = "upload",
    ) -> str:
        """Upload file to local storage."""
        key = f"users/{user_id}/sessions/{session_id}/{file_type}s/{filename}"
        file_path = self.base_dir / key
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_bytes(file_data)
        return key

    def download_file(self, key: str) -> bytes:
        """Download file from local storage."""
        file_path = self.base_dir / key
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {key}")
        return file_path.read_bytes()

    def generate_download_url(self, key: str, expires_in: int = 3600) -> str:
        """Generate local file path (no actual URL)."""
        return str(self.base_dir / key)
