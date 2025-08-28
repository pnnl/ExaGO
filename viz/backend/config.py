import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Database Configuration
sql_key = os.getenv('SQL_KEY', 'your_database_password')
database_name = os.getenv('DATABASE_NAME', 'exago_viz')

# OpenAI Configuration
openai_key = os.getenv('OPENAI_API_KEY')

# development/production
environment = os.getenv('ENVIRONMENT', 'dev')

# Host Configuration
# westmap-database is docker container name
host = 'westmap-database' if environment == 'prod' else 'localhost'