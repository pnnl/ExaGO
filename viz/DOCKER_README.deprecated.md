# Westmap Docker Setup Guide

This guide explains how to run the Westmap visualization platform using Docker containers.

## Architecture Overview

The Westmap platform consists of three main components:
- **Frontend**: React-based visualization interface (Port 8080)
- **Backend**: Flask API with ChatGrid AI functionality (Port 5000)
- **Database**: PostgreSQL with PostGIS for spatial queries (Port 5432)

## Prerequisites

- Docker installed and running
- Docker network support
- `.env` files configured (see Configuration section)

## Quick Start

### 1. Build All Images
```bash
# From the viz directory
make -f Makefile.root build-all
```

### 2. Run Complete Stack
```bash
# Run database, backend, and frontend
make -f Makefile.root run-all
```

### 3. Access the Application
- **Frontend**: http://localhost:8080
- **Backend API**: http://localhost:5000
- **Database**: localhost:5432 (PostgreSQL)
- **ChatGrid**: Click the chat widget in the frontend

## Individual Component Management

### Database
```bash
cd viz/database
make build    # Build PostgreSQL image
make run      # Run database container
make clean    # Stop and remove database container
make logs     # View database logs
make setup    # Initialize database with sample data
```

### Frontend
```bash
cd viz
make build    # Build frontend image
make run      # Run frontend container
make clean    # Stop and remove frontend container
make logs     # View frontend logs
```

### Backend
```bash
cd viz/backend
make build    # Build backend image
make run      # Run backend container
make clean    # Stop and remove backend container
make logs     # View backend logs
```

## Configuration

### Environment Variables
The application uses `.env` files for configuration:

#### Frontend `.env` (viz/.env)
```bash
OPENAI_API_KEY=your_openai_api_key_here
SQL_KEY=your_database_password
DATABASE_NAME=exago_viz
DB_HOST=localhost
DB_PORT=5432
DB_USER=postgres
NODE_ENV=development
```

#### Backend `.env` (viz/backend/.env)
```bash
OPENAI_API_KEY=your_openai_api_key_here
SQL_KEY=your_database_password
DATABASE_NAME=exago_viz
DB_HOST=localhost
DB_PORT=5432
DB_USER=postgres
FLASK_ENV=development
FLASK_DEBUG=True
```

## ✅ Testing and Verification

### Complete Stack Verification
Test all three tiers are working correctly:

```bash
# 1. Test Frontend (Should return HTML page)
curl http://localhost:8080

# 2. Test Backend API (Should return JSON response)
curl -X POST -H "Content-Type: application/json" \
  -d '{"inputText":"How many power plants are there?"}' \
  http://localhost:5000/data

# 3. Test Database Connection (Should return 849)
docker exec -it westmap-database psql -U postgres -d exago_viz \
  -c "SELECT COUNT(*) FROM generation;"
```

### Expected Results
- **Frontend Test**: HTML page with "Westmap - Weather and Energy System..." title
- **Backend Test**: `{"result_list": [{"count": 849}], "text": "There are 849 power plants."}`
- **Database Test**: `849` (total generation records)

### Admin Dashboard Testing (New in V1.0.2)
Test the new admin authentication system:

```bash
# Access admin dashboard directly
curl http://localhost:8080/?admin=true

# Note: Admin functionality requires Firebase setup and authentication
# See ADMIN_SETUP.md for Firebase configuration instructions
```

### Admin Features Available:
- **User Registration & Approval**: Complete user management workflow
- **Real-time Dashboard**: Live updates of user registrations
- **Advanced User Management**: Approve, reject, promote, delete users
- **Filtering & Sorting**: Organize users by status, date, email
- **Role-based Access**: Admin-only dashboard access
- **Security Confirmations**: Safe deletion and promotion workflows

### Sample ChatGrid AI Queries
Test the natural language AI functionality:

```bash
# Power generation queries
curl -X POST -H "Content-Type: application/json" \
  -d '{"inputText":"What types of power generation are available?"}' \
  http://localhost:5000/data

curl -X POST -H "Content-Type: application/json" \
  -d '{"inputText":"How many wind power plants are there?"}' \
  http://localhost:5000/data

curl -X POST -H "Content-Type: application/json" \
  -d '{"inputText":"Show me the top 5 highest capacity power plants"}' \
  http://localhost:5000/data
```

### Database Statistics
The database contains the following data:
- **Generation**: 849 power plants (wind, solar, hydro, natural gas, coal, nuclear)
- **Bus**: 4,762 electrical bus records
- **Transmission Lines**: 12,706 transmission line records  
- **Counties**: 3,221 US county records with spatial data
- **States**: 52 US state records with spatial boundaries

### Container Status Check
```bash
# View all running containers
docker ps --filter "network=gradscaler-network"

# Check container logs
docker logs westmap-frontend --tail 20
docker logs westmap-backend --tail 20  
docker logs westmap-database --tail 20
```
make -f Makefile.root run-all          # Run complete stack (database, backend, frontend)
make -f Makefile.root stop-all         # Stop all containers
make -f Makefile.root clean-all        # Clean all containers and volumes
make -f Makefile.root push-all         # Push all images to registry
make -f Makefile.root logs-all         # Show all container logs
make -f Makefile.root status           # Show container status
make -f Makefile.root dev-restart      # Rebuild and restart everything
make -f Makefile.root setup-data       # Load sample data into database
```

### Frontend Commands (from viz directory)
```bash
make build    # Build frontend Docker image
make run      # Run frontend container
make clean    # Stop and remove container
make push     # Push image to registry
make logs     # Show container logs
make shell    # Shell into container
```

### Backend Commands (from viz/backend directory)
```bash
make build              # Build backend Docker image
make run                # Run backend container
make run-with-setup     # Run with database setup
make clean              # Stop and remove container
make push               # Push image to registry
make logs               # Show container logs
make shell              # Shell into container
make setup-db           # Setup database only
```

## Docker Images

The setup creates the following Docker images:
- `gradscaler/westmap-frontend:V1.0.2` - React frontend with admin authentication system
- `gradscaler/westmap-backend:V1.0.2` - Flask backend with Python 3.9
- `gradscaler/westmap-database:V1.0.2` - PostgreSQL 15 with PostGIS extension

## Network Configuration

All containers run on a custom Docker network: `gradscaler-network`
This allows containers to communicate with each other by container name.

## Data Volume Management

The database uses persistent volumes for data storage:
- `westmap-postgres-data` - PostgreSQL data directory
- CSV data files are automatically loaded during database initialization

## Data Files

The backend includes the following data files for ChatGrid functionality:
- `bus.csv` - Power system bus data
- `generation.csv` - Generator data
- `transmission_line.csv` - Transmission line data
- `counties.csv` - County boundary data
- `us states.csv` - US state data
- Sample JSON files for testing

## 🎉 Deployment Success

✅ **Complete 3-Tier Docker Stack Successfully Deployed**

**What's Running:**
- **PostgreSQL Database** (westmap-database) with PostGIS spatial extensions
- **Flask Backend API** (westmap-backend) with ChatGrid AI functionality
- **React Frontend** (westmap-frontend) with admin authentication and user management

**Data Loaded:**
- 849 power generation facilities (wind, solar, hydro, gas, coal, nuclear)
- 4,762 electrical bus records with transmission data
- 12,706 transmission line connections
- 3,221 US counties with spatial boundaries  
- 52 US states with geographic data
- **Total: 18,590+ records** ready for ChatGrid AI queries

**Verified Features:**
- ✅ Natural language queries to ChatGrid AI
- ✅ PostgreSQL spatial queries with PostGIS
- ✅ Real-time data visualization
- ✅ Cross-container networking and communication
- ✅ Persistent data storage with Docker volumes
- ✅ Firebase authentication system
- ✅ Admin user management dashboard
- ✅ User approval workflow

**Example Working Queries:**
- "How many power plants are there?" → 849 power plants
- "What types of power generation are available?" → wind, solar, hydro, ng, coal, nuclear  
- "How many wind power plants are there?" → 139 wind power plants

## Next Steps

1. **Open the Web Interface**: http://localhost:8080
2. **Setup Admin Authentication**: Follow ADMIN_SETUP.md for Firebase configuration
3. **Access Admin Dashboard**: http://localhost:8080/?admin=true (requires admin authentication)
4. **Click the ChatGrid widget** to start natural language queries
5. **Explore the data** using the visualization tools
6. **Query the API directly** for custom integrations

The complete Westmap platform with ChatGrid AI and Admin Authentication is now ready for use!
3. Ensure PostGIS extension is installed
4. Run database setup script

## Development Workflow

For development with hot reloading:

1. **Setup once**:
   ```bash
   make -f Makefile.root build-all
   make -f Makefile.root setup-db
   ```

2. **Daily development**:
   ```bash
   make -f Makefile.root dev-restart
   ```

3. **View logs**:
   ```bash
   make -f Makefile.root logs-all
   ```

4. **Check status**:
   ```bash
   make -f Makefile.root status
   ```

## Port Configuration

- **Frontend**: 8080 (webpack-dev-server)
- **Backend**: 5000 (Flask)
- **PostgreSQL**: 5432 (host system)

## Security Notes

- API keys are loaded from `.env` files
- Database credentials are configurable
- Production deployments should use secure passwords
- Consider using Docker secrets for production

## Customization

To customize for your setup:

1. Update the `USERNAME` variable in Makefiles
2. Modify `VERSION` for different image tags
3. Change port mappings as needed
4. Update database credentials in `.env` files

---

**Last Updated**: August 3, 2025  
**Docker Configuration**: Alpine-based images with minimal dependencies
