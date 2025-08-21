# ChatGrid Configuration Setup & Documentation

## 📋 Overview
ChatGrid is the AI-powered chat interface in the ExaGO visualization that allows users to query power grid data using natural language. It uses OpenAI's GPT models through LangChain to process queries and interact with a PostgreSQL database.

## 🏗️ Architecture
```
Frontend (React) → Backend (Flask) → LangChain → OpenAI GPT → PostgreSQL
```

**How ChatGrid Works:**
1. **Frontend**: Chat widget (`react-chat-widget`) sends user queries via POST to `http://localhost:5000/data`
2. **Backend**: Flask server (`server.py`) receives requests and processes them
3. **AI Processing**: LangChain converts natural language to SQL using OpenAI GPT-3.5-turbo
4. **Database**: PostgreSQL with PostGIS executes queries on power grid data
5. **Response**: Returns structured data + text summary that updates the visualization layers

## ✅ Configuration Changes Made

### 1. Environment Variables Setup
- ✅ Created `viz/.env` with OpenAI API key and database credentials
- ✅ Created `viz/backend/.env` for backend-specific variables
- ✅ Updated `backend/config.py` to use `python-dotenv` instead of hardcoded keys
- ✅ Added `python-dotenv==1.0.0` to `requirements.txt`

### 2. Security Improvements
- ✅ Removed hardcoded OpenAI API key from `config.py`
- ✅ All sensitive credentials now loaded from environment variables
- ✅ Environment files added to `.gitignore` (should be added)

### 3. Database Setup Preparation
- ✅ Created `backend/setup_database.py` automated setup script
- ✅ Configured default database name as `exago_viz`
- ✅ Set default database password as `exago2024`

## 📁 Files Created/Modified

### Configuration Files:
```
viz/.env                     # Frontend environment variables
viz/backend/.env            # Backend environment variables  
viz/backend/config.py       # Updated to use python-dotenv
viz/backend/requirements.txt # Added python-dotenv dependency
viz/backend/setup_database.py # Database setup automation
```

### Available Data Files:
```
viz/data/
├── bus.csv              # Bus/substation data
├── generation.csv       # Generator data  
├── transmission_line.csv # Transmission line data
├── counties.csv         # County boundaries
├── us states.csv        # State boundaries
├── case_ACTIVSg200.json # 200-bus sample network
├── case_ACTIVSg500.json # 500-bus sample network
└── case_ACTIVSg10k.json # 10k-bus sample network
```

## 🚀 Setup Instructions

### Prerequisites
- ✅ Node.js v24.4.1 (installed)
- ✅ npm v11.4.2 (installed)
- ⚠️ Python 3.7+ (needs verification)
- ⚠️ PostgreSQL 12+ (needs installation)

### Step 1: Backend Dependencies
```bash
cd viz/backend
pip install -r requirements.txt
```

### Step 2: PostgreSQL Database Setup

#### Option A: Automated Setup (Recommended)
```bash
cd viz/backend
python setup_database.py
```

#### Option B: Manual Setup
1. **Install PostgreSQL**:
   - Download from: https://www.postgresql.org/download/
   - During installation, set password to `exago2024` (or update `.env`)
   - Remember the username (usually `postgres`)

2. **Create Database**:
   ```sql
   CREATE DATABASE exago_viz;
   \c exago_viz;
   CREATE EXTENSION IF NOT EXISTS postgis;
   ```

3. **Import Data**:
   Use pgAdmin or command line to import CSV files:
   - `generation.csv` → `generation` table
   - `bus.csv` → `bus` table  
   - `transmission_line.csv` → `transmission_line` table
   - `counties.csv` → `counties` table
   - `us states.csv` → `us_states` table

### Step 3: Start Backend Server
```bash
cd viz/backend
python server.py
```
Server will start on: http://localhost:5000

### Step 4: Test ChatGrid
1. Frontend should still be running on http://localhost:8080
2. Click the ChatGrid widget (bottom-left blue chat icon)
3. Try queries like:
   - "Show me all generators with capacity higher than 100"
   - "List transmission lines in the network" 
   - "How many buses are in the system?"

## 🔧 Environment Variables Reference

### Frontend (.env)
```bash
OPENAI_API_KEY=openai_api_key
SQL_KEY=sql_key
DATABASE_NAME=exago_viz
DB_HOST=localhost
DB_PORT=5432
DB_USER=postgres
NODE_ENV=development
```

### Backend (.env)
```bash
OPENAI_API_KEY=openai_api_key
SQL_KEY=sql_key
DATABASE_NAME=exago_viz
DB_HOST=localhost
DB_PORT=5432
DB_USER=postgres
FLASK_ENV=development
FLASK_DEBUG=True
```

## 🧪 Testing & Validation

### Backend Health Check
```bash
curl -X POST http://localhost:5000/data \
  -H "Content-Type: application/json" \
  -d '{"inputText": "How many generators are there?"}'
```

### Database Connection Test
```python
from backend.config import sql_key, database_name, openai_key
print(f"DB: {database_name}, API Key: {openai_key[:10]}...")
```

## 🐛 Troubleshooting

### Common Issues:
1. **ImportError: No module named 'dotenv'**
   - Solution: `pip install python-dotenv`

2. **Database connection refused**
   - Check PostgreSQL is running: `pg_ctl status`
   - Verify credentials in `.env` files
   - Test connection: `psql -U postgres -d exago_viz`

3. **OpenAI API errors**
   - Verify API key has sufficient credits
   - Check key format in `.env` files

4. **ChatGrid not responding**
   - Ensure backend server is running on port 5000
   - Check browser network tab for CORS errors
   - Verify Flask-CORS is installed

## 📊 Current Status

### ✅ Completed:
- [x] Frontend visualization running
- [x] OpenAI API key configured
- [x] Environment variables setup
- [x] Backend configuration updated
- [x] Security improvements implemented
- [x] Database setup scripts created
- [x] Documentation completed

### ⚠️ Pending:
- [ ] PostgreSQL installation
- [ ] Backend dependencies installation
- [ ] Database data import
- [ ] Backend server startup
- [ ] End-to-end ChatGrid testing

## 🔒 Security Notes
- API keys stored in `.env` files (should be git-ignored)
- Never commit `.env` files to version control
- Database credentials use environment variables
- Default password should be changed in production

## 📖 References
- [Original ExaGO README](../README.md)
- [PostgreSQL Installation Guide](https://www.postgresql.org/download/)
- [LangChain Documentation](https://python.langchain.com/docs/get_started/introduction.html)
- [OpenAI API Documentation](https://platform.openai.com/docs/models/overview)

---

**Last Updated**: August 2, 2025  
**Configuration By**: GitHub Copilot Assistant  
**Project**: ExaGO Visualization ChatGrid
