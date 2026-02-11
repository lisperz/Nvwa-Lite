"""Nvwa-Lite: scRNA-seq Visualization Agent.

Minimal entry point — launches the Streamlit app.
"""

import subprocess
import sys


def main() -> None:
    subprocess.run(
        [
            sys.executable, "-m", "streamlit", "run",
            "src/ui/app.py",
            "--server.port=8501",
            "--server.address=0.0.0.0",
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
