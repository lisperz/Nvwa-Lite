# Fix Cloudflare Error 521 - Web Server Down

## Problem
Your site shows "Web server is down (Error 521)" because Cloudflare cannot connect to your EC2 instance.

## Root Cause
AWS Security Group is blocking Cloudflare's origin server IPs.

## Solution

### Step 1: Update EC2 Security Group

1. Go to AWS Console → EC2 → Security Groups
2. Find your instance's security group (check instance `3.150.203.87`)
3. Edit **Inbound Rules**
4. Add these rules:

**For HTTP (port 80):**
- Type: HTTP
- Protocol: TCP
- Port: 80
- Source: Custom → `0.0.0.0/0` (or Cloudflare IP ranges)

**For HTTPS (port 443) - if using:**
- Type: HTTPS
- Protocol: TCP
- Port: 443
- Source: Custom → `0.0.0.0/0`

### Step 2: Verify Cloudflare Settings

1. Go to Cloudflare Dashboard → nvwa.bio
2. Check **SSL/TLS** settings:
   - Set to "Flexible" (Cloudflare → HTTP → Origin)
   - Or "Full" if you have SSL on EC2

3. Check **DNS** settings:
   - Ensure A record points to `3.150.203.87`
   - Orange cloud (proxied) should be enabled

### Step 3: Test Connection

```bash
# From your local machine
curl -I http://3.150.203.87

# Should return HTTP 200 or 403 (not connection refused)
```

## Quick Fix (Temporary)

If you need immediate access, bypass Cloudflare:
1. Go to Cloudflare DNS settings
2. Click the orange cloud next to your A record to turn it gray (DNS only)
3. Wait 1-2 minutes for DNS propagation
4. Access site directly via `http://3.150.203.87`

## Cloudflare IP Ranges (Optional - More Secure)

Instead of `0.0.0.0/0`, use Cloudflare's IP ranges:
https://www.cloudflare.com/ips/

Add each range as a separate inbound rule.
