import os
from agno.agent import Agent
from agno.models.openai import OpenAIChat
from agno.db.sqlite import SqliteDb
from agno.db.postgres import PostgresDb
from agno.os import AgentOS

# from agno.os.interfaces.agui import AGUI
from agno.tools.reasoning import ReasoningTools
from agno.knowledge.knowledge import Knowledge
from agno.vectordb.pgvector import PgVector
from agno.models.anthropic import Claude

from tools import ShellTools

# supabase connection string, pass b7Orv6W070BSQCph
SUPABASE_DB_URL = "postgresql://postgres.ytflxbsklrskeswzolkm:b7Orv6W070BSQCph@aws-0-us-west-2.pooler.supabase.com:6543/postgres"


def _get_model():
    """Resolve model and API key from environment variables.

    Environment:
      - OPENAI_API_KEY: API key for OpenAI
      - NVWA_OPENAI_MODEL: Optional, defaults to 'gpt-4.1'
    """
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError(
            "OPENAI_API_KEY is not set. Please export it in your environment."
        )
    model_id = os.getenv("NVWA_OPENAI_MODEL", "gpt-5.1")
    # model_id = os.getenv("claude-sonnet-4-5", "gpt-5")
    model_id = "claude-sonnet-4-5"
    api_key = os.getenv("ANTHROPIC_API_KEY", api_key)
    return Claude(id=model_id, api_key=api_key)


AGENT_INSTRUCTIONS = """
    You are Nvwa Bioinformatics Agent.
    
    Specialization:
    - Bioinformatics data analysis and data visualization.
    - Comfortable writing, running, and debugging Python and R code.
    - Proficient with Linux shell.

    Runtime & workspace:
    - Runs in Ubuntu 24 with Python and R environments pre-configured.
    - You have limited sudo permissions for installing missing packages.
    - The working directory inside the container is /workspace with:
        - /workspace/src     -> source code (Python/R)
        - /workspace/data    -> input data
        - /workspace/outputs -> results and reports

    Tools:
    - Use the provided tools to list, read, and write files under /workspace.
    - Use the tools to run Python/R code and scripts, and to install packages.
    - Prefer creating reusable scripts in /workspace/src and writing outputs under /workspace/outputs.
    - When generating plots, save them to /workspace/outputs with clear filenames.

    Behavioral guidelines:
    - Before running, check that required packages are installed; install if missing.
    - Keep analysis reproducible: pin versions where practical and record commands used.
    - For large datasets, stream or sample preview (head) to avoid overwhelming output.
    - Always explain what you will run and where outputs will be written.
    """.strip()


_container = os.getenv("NVWA_CONTAINER_NAME", "nvwa")

# Create a knowledge base with Postgres and PgVector
# knowledge = Knowledge(
#     contents_db=PostgresDb(
#         db_url=SUPABASE_DB_URL,
#         knowledge_table="knowledge_contents",
#     ),
#     vector_db=PgVector(
#         table_name="knowledge_documents",
#         db_url=SUPABASE_DB_URL,
#     ),
# )

agent = Agent(
    name="Nvwa Bioinformatics Agent",
    description="An AI agent specialized in bioinformatics data analysis and visualization.",
    instructions=AGENT_INSTRUCTIONS,
    model=_get_model(),
    # knowledge=knowledge,
    search_knowledge=True,
    tools=[ShellTools(container_name=_container), ReasoningTools()],
    markdown=True,
    db=SqliteDb(db_file="agno.db"),
    # db=PostgresDb(
    #     db_url=SUPABASE_DB_URL,
    #     # db_schema=""
    #     session_table="agent_sessions",
    #     memory_table="agent_memory",
    #     metrics_table="agent_metrics",
    #     eval_table="agent_evals",
    #     knowledge_table="agent_knowledge",
    # ),
    add_history_to_context=True,
    enable_session_summaries=True,
)

# Create the AgentOS and expose a FastAPI app
agent_os = AgentOS(
    agents=[agent],
    # interfaces=[AGUI(agent=agent)]
)
app = agent_os.get_app()

# If executed directly, print a simple readiness message
if __name__ == "__main__":
    print(
        "Nvwa Bioinformatics Agent is configured. Serve with: uvicorn ai.main:app --reload --port 8000"
    )
