# Feedback Widget Implementation

**Status:** ✅ Implemented (2026-03-28)
**Deadline:** March 28, 2026

## Overview

In-app feedback collection widget that appears after 10 minutes of user inactivity when at least one plot has been generated. Collects 3 key metrics for Demo Day (May 15).

## Implementation

### 1. Database Schema

**Table:** `feedback_responses`

```sql
CREATE TABLE feedback_responses (
    response_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id TEXT NOT NULL,
    user_id TEXT NOT NULL,
    q1_score INTEGER NOT NULL CHECK (q1_score BETWEEN 1 AND 5),
    q2_time_saved VARCHAR(50) NOT NULL,
    q3_open_text TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);
```

**Migration:** `migrations/002_add_feedback_responses.sql`

### 2. Components

**Files Created:**
- `src/ui/feedback_widget.py` - HTML/JS widget component
- `migrations/002_add_feedback_responses.sql` - Database schema

**Files Modified:**
- `src/ui/app.py` - Widget integration and feedback handling
- `src/agent/tools.py` - Plot generation tracking
- `src/db/logger.py` - Feedback logging method

### 3. Widget Behavior

**Trigger Conditions:**
- User inactive for 10 minutes (no mouse/keyboard/scroll)
- At least one plot generated in session
- Feedback not yet submitted for this session

**UI Features:**
- Slide-in from bottom-right corner (320px width)
- Non-blocking overlay (z-index: 9999)
- Auto-minimize to icon after 30 seconds
- Show once per session
- Skip button to close immediately

**Questions:**
1. **Q1:** Star rating (1-5) - AI interpretation accuracy
2. **Q2:** Single select - Time saved (4 options)
3. **Q3:** Optional text - Issues or surprises

### 4. Data Flow

```
User inactive 10min + plot exists
    ↓
Widget appears (slide-in animation)
    ↓
User fills Q1 (required) + Q2 (required) + Q3 (optional)
    ↓
Submit → DatabaseLogger.log_feedback()
    ↓
Save to RDS feedback_responses table
    ↓
Thank you message → Close after 3s
```

### 5. Testing Checklist

- [ ] Widget does NOT appear if no plot generated
- [ ] Widget appears after 10min inactivity with plot
- [ ] Skip button closes widget immediately
- [ ] Q1 star rating saves as integer 1-5
- [ ] Q2 selection saves as string enum
- [ ] Q3 text saves correctly (NULL if blank)
- [ ] Widget minimizes after 30s no interaction
- [ ] Widget shows once per session only
- [ ] Thank you message displays after submit
- [ ] Test with PBMC 3k dataset

## Deployment

**Local Testing:**
```bash
cd nvwa-lite
docker-compose restart nvwa-lite
# Open http://localhost:8501/app
```

**Production Deployment:**
```bash
# Run migration on EC2
ssh -i nvwa-key.pem ubuntu@3.150.203.87
cd Nvwa-Lite
bash scripts/migrate_db.sh

# Deploy code
bash scripts/deploy-to-ec2.sh
```

## Demo Day Metrics

This widget collects the data needed for:
- **Quality metric:** "X% of analyses rated 4+ by expert PIs"
- **Value prop:** "Average time saved: Y hours per analysis"
- **Product feedback:** Qualitative insights for roadmap

Every response = one data point toward funding story.
