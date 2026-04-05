# Fix: Change Cloudflare SSL Mode

## Current Issue
- Cloudflare SSL/TLS: **Full (strict)** ❌
- Your EC2: Only HTTP (port 80), no HTTPS ❌
- Result: Cloudflare can't connect → Error 521

## Solution

### Step 1: Change SSL/TLS Mode in Cloudflare

1. Go to Cloudflare Dashboard → nvwa.bio
2. Click **SSL/TLS** in left sidebar
3. Under "SSL/TLS encryption mode", select **Flexible**
4. Wait 1-2 minutes for changes to propagate

### What Each Mode Means

**Flexible** (Use this):
- Browser → Cloudflare: HTTPS ✅
- Cloudflare → Your server: HTTP ✅
- Works with your current setup

**Full**:
- Browser → Cloudflare: HTTPS
- Cloudflare → Your server: HTTPS (requires port 443)
- You don't have HTTPS on EC2 ❌

**Full (strict)** (Current - causing error):
- Same as Full, but requires valid SSL certificate
- You don't have SSL certificate on EC2 ❌

### Step 2: Verify

After changing to Flexible, test:
```bash
curl -I https://nvwa.bio
```

Should return HTTP 200 or 301/302 redirect.

## Alternative: Add HTTPS to EC2 (More Secure)

If you want to use "Full" mode later:
1. Install SSL certificate on EC2 (Let's Encrypt)
2. Configure nginx to listen on port 443
3. Open port 443 in security group
4. Change Cloudflare to "Full" mode
