# Backend Environment Setup

## Environment Configuration

This backend requires environment variables to be configured for proper operation.

### Setup Instructions

1. **Copy the template file:**
   ```bash
   cp env.template env
   ```

2. **Edit the `env` file with your actual credentials:**
   - Replace `your_openai_api_key_here` with your OpenAI API key
   - Replace `your_database_password_here` with your PostgreSQL password
   - Replace Firebase configuration values with your actual Firebase project details
   - Replace Mapbox token if using map features

### Required Services

- **OpenAI API**: For ChatGrid natural language processing
- **PostgreSQL**: For data storage and querying
- **Firebase**: For authentication and real-time features
- **Mapbox** (optional): For enhanced map visualizations

### Security Notes

- The `env` file is excluded from version control via `.gitignore`
- Never commit API keys or sensitive credentials to the repository
- Keep your API keys secure and rotate them regularly

### File Structure

- `env.template` - Template with placeholder values (safe to commit)
- `env` - Your actual environment file (excluded from git)
- `server.py` - Main Flask application
- `requirements.txt` - Python dependencies
