Project: Nvwa-Bio MVP (scRNA visualization) — Regression Test v0
Branch: feat/regression-test-v0

Goal (this branch):
- Implement a minimal automated regression test harness v0:
  - tests.yaml with ~20 cases (happy path + input errors + out-of-scope + scale)
  - scripts/run_tests.py to execute the cases against local/dev API
  - output a markdown report with pass/fail and failure reasons

Constraints:
- Do NOT refactor architecture.
- Keep PR small and focused.
- Prefer deterministic checks (artifacts exist, status codes, key error strings).
- Do not add heavy UI testing (no Playwright/Selenium in v0).

Repo conventions:
- Put test cases in tests/tests.yaml
- Put runner in scripts/run_tests.py
- If API endpoints are missing, add minimal internal endpoints in src/api (or equivalent).
- Each test case should define: dataset_id, prompt, expected_status, expected_artifacts, must_contain.