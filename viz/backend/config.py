import os

sql_key = "YOUR_DATABASE_PASSWORD"
database_name = "YOUR_DATABASE_NAME"
openai_key = "YOUR_OPENAI_KEY"

try:
    sql_key = os.getenv("POSTGRES_PASSWORD")
    database_name = os.getenv("POSTGRES_DB")
    openai_key = os.getenv("OPENAI_KEY")
except:
    pass