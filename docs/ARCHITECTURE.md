# Nvwa MVP Architecture - Phase 1 & Phase 2

## Executive Summary

This document describes the architecture transformation from the current monolithic Streamlit app to a scalable MVP system that supports 5-10 pilot users (Phase 1) and can scale to 50+ users (Phase 2).

**Key Insight:** Phase 1 implements clean service boundaries within the Streamlit app. Phase 2 extracts these services into a distributed system. This approach ships fast while maintaining a clear path to scale.

## Current State (Before Phase 1)

```
┌─────────────────────────────────────┐
│      Streamlit Monolith             │
│                                     │
│  ┌──────────────────────────────┐  │
│  │  UI + Agent + Plotting       │  │
│  │  (all in one process)        │  │
│  └──────────────────────────────┘  │
│                                     │
│  - No authentication               │
│  - No structured logging           │
│  - Local file uploads              │
│  - No session management           │
│  - No monitoring                   │
└─────────────────────────────────────┘
```

**Limitations:**
- No user tracking
- No concurrency control
- No cost monitoring
- Not scalable beyond single instance

## Phase 1 Architecture (Pilot: 5-10 Users)

```
┌──────────────────────────────────────────────────────────────┐
│                 Streamlit App (Port 8501)                     │
│                                                               │
│  ┌─────────────┐  ┌──────────────┐  ┌──────────────────┐   │
│  │   Auth      │  │   Session    │  │   Analysis       │   │
│  │   Service   │  │   Manager    │  │   Engine         │   │
│  │             │  │              │  │   (existing)     │   │
│  │ - Token     │  │ - Redis      │  │                  │   │
│  │   validate  │  │ - Concurrency│  │ - Agent          │   │
│  │ - User      │  │   limits     │  │ - Tools          │   │
│  │   tracking  │  │ - Timeouts   │  │ - Plotting       │   │
│  └──────┬──────┘  └──────┬───────┘  └────────┬─────────┘   │
│         │                │                    │              │
│  ┌──────▼────────────────▼────────────────────▼──────────┐  │
│  │           Event Logger (structured JSON)              │  │
│  │  - tool_execution.log                                 │  │
│  │  - user_interaction.log                               │  │
│  │  - system_metrics.log                                 │  │
│  └──────────────────────┬────────────────────────────────┘  │
│                         │                                    │
│  ┌──────────────────────▼────────────────────────────────┐  │
│  │           Storage Service                             │  │
│  │  - S3 (primary)                                       │  │
│  │  - Local (fallback)                                   │  │
│  └───────────────────────────────────────────────────────┘  │
└───────────────────────────────────────────────────────────────┘
                         │
            ┌────────────▼────────────┐
            │  Dashboard (Port 8502)  │
            │  - Active users         │
            │  - Token usage          │
            │  - Tool success rates   │
            └─────────────────────────┘
```

### Key Components

#### 1. Authentication Service
- **Purpose:** Token-based user authentication
- **Implementation:** Pre-generated tokens for pilot users
- **Storage:** Environment variables or JSON file
- **Scalability:** Can migrate to JWT/OAuth in Phase 2

#### 2. Session Manager
- **Purpose:** Track user sessions and enforce concurrency limits
- **Implementation:** Redis-backed with in-memory fallback
- **Limits:** Max 1 session per user, max 2 concurrent system-wide
- **Timeout:** 30 minutes of inactivity

#### 3. Event Logger
- **Purpose:** Structured logging for monitoring and analytics
- **Format:** JSON with user_id, session_id, timestamp, event, payload
- **Files:** Separate logs for tools, interactions, metrics
- **Retention:** Indefinite (for analytics)

#### 4. Storage Service
- **Purpose:** Manage dataset uploads and result storage
- **Primary:** S3 with pre-signed URLs
- **Fallback:** Local filesystem for development
- **Structure:** `users/<user_id>/sessions/<session_id>/...`
- **TTL:** 14 days automatic cleanup

#### 5. Analysis Engine (Existing)
- **Purpose:** Core scRNA-seq analysis
- **Components:** Agent, Tools, Plotting
- **No changes:** Existing code preserved
- **Integration:** Wrapped with logging and session tracking

### Data Flow

```
1. User accesses http://3.150.203.87/?token=ABC123
2. AuthService validates token → user_id
3. SessionManager creates session → session_id
4. User uploads file → S3StorageService → S3 key
5. User sends message → Agent processes
6. Tools execute → EventLogger logs
7. Results generated → S3StorageService stores
8. Response returned → EventLogger logs
9. Dashboard reads logs → displays metrics
```

### Concurrency Control

```
User A uploads dataset
  ↓
SessionManager checks:
  - User A has no active session? ✓
  - System has < 2 active sessions? ✓
  ↓
Session created, analysis proceeds

User B uploads dataset (while A active)
  ↓
SessionManager checks:
  - User B has no active session? ✓
  - System has < 2 active sessions? ✓
  ↓
Session created, analysis proceeds

User C uploads dataset (while A & B active)
  ↓
SessionManager checks:
  - User C has no active session? ✓
  - System has < 2 active sessions? ✗
  ↓
Session rejected, user sees "System busy" message
```

### Monitoring & Analytics

**Real-time Metrics (CEO Dashboard):**
- Active users (30min window)
- Token consumption (24h)
- Estimated cost
- Tool success rates
- Session depth (avg messages per session)
- Queue patience (users waiting)

**Log Analysis:**
- Intent routing accuracy
- Semantic mapping hit rate
- Task cycle time
- Artifact download rate
- Self-correction traces
- Parameter iteration patterns

## Phase 2 Architecture (Scale: 50+ Users)

```
┌──────────────────────────────────────────────────────────────┐
│                   Control Plane (EC2)                         │
│                                                               │
│  ┌──────────────┐         ┌─────────────────┐               │
│  │  Streamlit   │────────▶│   FastAPI       │               │
│  │  Frontend    │◀────────│   Backend       │               │
│  │  (Port 8501) │         │   (Port 8000)   │               │
│  └──────────────┘         │                 │               │
│                           │  - Auth API     │               │
│                           │  - Sessions API │               │
│                           │  - Jobs API     │               │
│                           │  - S3 URLs API  │               │
│                           │  - WebSocket    │               │
│                           └────┬────────────┘               │
│                                │                             │
│                           ┌────▼────────┐                   │
│                           │  PostgreSQL │                   │
│                           │  - Users    │                   │
│                           │  - Sessions │                   │
│                           │  - Jobs     │                   │
│                           │  - Analytics│                   │
│                           └────┬────────┘                   │
│                                │                             │
│                           ┌────▼────────┐                   │
│                           │ Redis Queue │                   │
│                           │  (Bull/BullMQ)│                 │
│                           └────┬────────┘                   │
│                                │                             │
│                           ┌────▼────────┐                   │
│                           │ Job Worker  │                   │
│                           │  (max N     │                   │
│                           │   concurrent)│                  │
│                           └────┬────────┘                   │
│                                │                             │
│                           ┌────▼────────┐                   │
│                           │  Docker     │                   │
│                           │  Engine     │                   │
│                           │  (per-session│                  │
│                           │   containers)│                  │
│                           └────┬────────┘                   │
└────────────────────────────────┼──────────────────────────────┘
                                 │
                        ┌────────▼────────┐
                        │    AWS S3       │
                        │  - Datasets     │
                        │  - Results      │
                        │  - Logs         │
                        └─────────────────┘
```

### Key Changes from Phase 1

#### 1. FastAPI Backend Extraction
- **Auth:** JWT tokens, OAuth integration
- **Sessions:** RESTful API for session management
- **Jobs:** Job submission and status tracking
- **WebSocket:** Real-time status updates

#### 2. PostgreSQL Database
- **Users:** Persistent user accounts
- **Sessions:** Session history and analytics
- **Jobs:** Job queue and execution history
- **Analytics:** Aggregated metrics

#### 3. Job Queue (Redis + Workers)
- **Queue:** Bull/BullMQ for job management
- **Workers:** Multiple worker processes
- **Concurrency:** Configurable per-worker
- **Priorities:** User-based or task-based

#### 4. Docker Container Isolation
- **Per-session containers:** Each session gets dedicated container
- **Resource limits:** CPU, memory, timeout enforcement
- **Isolation:** No cross-user interference
- **Cleanup:** Automatic container destruction

#### 5. Scalability Improvements
- **Horizontal scaling:** Multiple worker nodes
- **Load balancing:** Nginx or ALB
- **Auto-scaling:** Based on queue depth
- **Monitoring:** CloudWatch, Prometheus

### Migration Path (Phase 1 → Phase 2)

**Step 1: Extract API Layer**
- Wrap existing services with FastAPI endpoints
- Streamlit calls API instead of direct service calls
- No logic changes, just interface extraction

**Step 2: Add PostgreSQL**
- Migrate session state from Redis to PostgreSQL
- Keep Redis for job queue only
- Add user account management

**Step 3: Add Job Queue**
- Implement job submission API
- Create worker process
- Migrate analysis execution to workers

**Step 4: Add Container Isolation**
- Package analysis engine as Docker image
- Worker launches containers per job
- Implement resource limits and cleanup

**Step 5: Scale Infrastructure**
- Add more worker nodes
- Implement load balancing
- Add auto-scaling policies

## Technology Stack

### Phase 1
- **UI:** Streamlit 1.40+
- **Agent:** LangChain + OpenAI gpt-4o-mini
- **Analysis:** Scanpy, Matplotlib
- **Auth:** Token-based (pre-shared secrets)
- **Session:** Redis (with in-memory fallback)
- **Storage:** AWS S3 (with local fallback)
- **Logging:** Python logging (JSON format)
- **Monitoring:** Streamlit dashboard

### Phase 2 (Additional)
- **Backend:** FastAPI
- **Database:** PostgreSQL
- **Queue:** Redis + Bull/BullMQ
- **Containers:** Docker Engine
- **Orchestration:** Docker Compose → ECS/Fargate
- **Monitoring:** CloudWatch, Prometheus, Grafana
- **Auth:** JWT, OAuth 2.0

## Cost Analysis

### Phase 1 (5-10 Users)

**AWS Infrastructure:**
- EC2 t3.xlarge: $120/month (on-demand) or $75/month (reserved)
- S3 storage: $5-10/month (100GB, 14-day TTL)
- Data transfer: $5-10/month
- **Total:** $130-145/month (on-demand) or $90-100/month (reserved)

**OpenAI API:**
- gpt-4o-mini: $0.15/1M input, $0.60/1M output
- Estimated: $20-50/month for 5-10 active users

**Total Phase 1:** $150-195/month (on-demand) or $110-150/month (reserved)

### Phase 2 (50+ Users)

**AWS Infrastructure:**
- EC2 t3.2xlarge (control plane): $255/month (on-demand) or $160/month (reserved)
- ECS Fargate (workers): $130-260/month (2-4 tasks)
- RDS PostgreSQL: $50-100/month (db.t3.medium)
- ElastiCache Redis: $15-30/month
- S3 storage: $20-50/month (500GB)
- Data transfer: $20-50/month
- **Total:** $490-745/month (on-demand) or $395-650/month (reserved)

**OpenAI API:**
- Estimated: $200-500/month for 50 active users

**Total Phase 2:** $690-1,245/month (on-demand) or $595-1,150/month (reserved)

## Security Considerations

### Phase 1
- **Authentication:** Pre-shared tokens (treat as passwords)
- **Authorization:** Token = full access (no role-based access)
- **Data isolation:** Per-user S3 namespaces
- **Encryption:** HTTPS for transport, S3 encryption at rest
- **Secrets:** Environment variables, not in code
- **Logs:** Contain user data, restrict access

### Phase 2 (Additional)
- **Authentication:** JWT tokens with expiration
- **Authorization:** Role-based access control (RBAC)
- **API security:** Rate limiting, CORS, input validation
- **Container security:** Resource limits, network isolation
- **Database security:** Encrypted connections, least privilege
- **Audit logs:** Immutable audit trail

## Operational Considerations

### Phase 1
- **Deployment:** Docker Compose on single EC2
- **Monitoring:** Log files + Streamlit dashboard
- **Backups:** S3 versioning, no database backups needed
- **Updates:** Rolling restart (brief downtime acceptable)
- **Debugging:** SSH to EC2, check logs, inspect Redis

### Phase 2
- **Deployment:** ECS/Fargate with blue-green deployment
- **Monitoring:** CloudWatch, Prometheus, Grafana
- **Backups:** RDS automated backups, S3 versioning
- **Updates:** Zero-downtime rolling updates
- **Debugging:** Centralized logging, distributed tracing

## Success Metrics

### Phase 1 KPIs
- **Reliability:** >95% tool execution success rate
- **Performance:** <5s average response time
- **Availability:** >99% uptime (pilot acceptable)
- **Cost:** <$200/month total
- **User satisfaction:** Qualitative feedback from pilot users

### Phase 2 KPIs
- **Reliability:** >99% tool execution success rate
- **Performance:** <3s average response time
- **Availability:** >99.9% uptime
- **Scalability:** Support 50+ concurrent users
- **Cost efficiency:** <$20/user/month

## Conclusion

**Phase 1** provides a solid foundation for the pilot with:
- Clean service boundaries
- Comprehensive monitoring
- Graceful degradation
- Clear migration path

**Phase 2** scales the system when needed with:
- Distributed architecture
- Container isolation
- Horizontal scalability
- Production-grade reliability

The key insight is that Phase 1 is not a throwaway prototype - it's a well-architected system that can evolve into Phase 2 without rewrites. The service interfaces defined in Phase 1 become the API contracts in Phase 2.

**Recommendation:** Implement Phase 1 now (2-3 days), launch pilot, gather feedback, then decide if/when Phase 2 is needed based on actual usage patterns and scaling requirements.
