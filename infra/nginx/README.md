# Production Nginx Configuration

This directory contains the production nginx configuration for the nvwa.bio host.

## File

- `nvwa.bio.conf` — Production nginx config for EC2 host (3.150.203.87)

## Deployment to EC2

The nginx config must be deployed to the EC2 host and loaded by the host nginx service (not Docker).

### Manual Deployment Steps

```bash
# 1. Copy config to EC2
scp -i /path/to/nvwa-key.pem infra/nginx/nvwa.bio.conf ubuntu@3.150.203.87:/tmp/

# 2. SSH into EC2
ssh -i /path/to/nvwa-key.pem ubuntu@3.150.203.87

# 3. Install config
sudo cp /tmp/nvwa.bio.conf /etc/nginx/sites-available/nvwa.bio
sudo ln -sf /etc/nginx/sites-available/nvwa.bio /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default

# 4. Validate config
sudo nginx -t

# 5. Reload nginx (zero-downtime)
sudo systemctl reload nginx
```

### Automated Deployment

The following scripts automate nginx deployment:

- `scripts/deploy_ssl.sh` — Initial SSL setup + nginx config deployment
- `scripts/update_nginx.sh` — Update nginx config only

Both scripts copy `infra/nginx/nvwa.bio.conf` to the EC2 host.

## Architecture Notes

- Host nginx runs on EC2 host (not in Docker)
- Proxies to Docker containers on localhost ports:
  - `/app` → `localhost:8501` (nvwa-lite Streamlit)
  - `/admin` → `localhost:8502` (dashboard Streamlit)
  - `/api/cleanup` → `localhost:8503` (cleanup endpoint, localhost-only binding)
- Serves static landing page from `/var/www/nvwa.bio/`
- Handles SSL termination (Let's Encrypt certs)
- Handles WebSocket upgrades for Streamlit

See `product/infra_architecture.md` for full architecture diagram.
