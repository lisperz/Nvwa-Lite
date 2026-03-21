#!/bin/bash
# Generate pilot user tokens for MVP deployment

set -e

echo "🔐 Generating pilot user tokens..."

python3 << 'EOF'
import json
import secrets
from datetime import datetime

# Define pilot users
pilot_users = [
    ("beta_user_001", "user1@example.com"),
    ("beta_user_002", "user2@example.com"),
    ("beta_user_003", "user3@example.com"),
    ("beta_user_004", "user4@example.com"),
    ("beta_user_005", "user5@example.com"),
    ("beta_user_006", "user6@example.com"),
    ("beta_user_007", "user7@example.com"),
    ("beta_user_008", "user8@example.com"),
    ("beta_user_009", "user9@example.com"),
    ("beta_user_010", "user10@example.com"),
]

tokens = {}
env_lines = []

for user_id, email in pilot_users:
    token = secrets.token_urlsafe(32)
    tokens[user_id] = {
        "email": email,
        "token": token,
        "access_url": f"http://3.150.203.87/?token={token}",
        "created_at": datetime.utcnow().isoformat()
    }

    # Generate .env format
    env_lines.append(f'PILOT_TOKEN_{user_id.split("_")[-1]}={user_id}:{email}:{token}')

# Save to JSON file
with open("pilot_tokens.json", "w") as f:
    json.dump(tokens, f, indent=2)

# Save to .env format
with open("pilot_tokens.env", "w") as f:
    f.write("# Pilot user tokens - add these to your .env file\n")
    f.write("# Format: PILOT_TOKEN_XXX=user_id:email:token\n\n")
    f.write("\n".join(env_lines))

print("✅ Tokens generated successfully!")
print("\nFiles created:")
print("  - pilot_tokens.json (full details)")
print("  - pilot_tokens.env (for .env file)")
print("\n📧 Send access URLs to users:")
for user_id, info in tokens.items():
    print(f"  {info['email']}: {info['access_url']}")

EOF

echo ""
echo "✅ Done! Next steps:"
echo "1. Review pilot_tokens.json"
echo "2. Add pilot_tokens.env contents to your .env file"
echo "3. Send access URLs to pilot users via email"
