# WestMap Production Deployment Guide (CSV-based Architecture)

## Overview
This guide documents the production deployment of the WestMap application - a simplified 2-tier Docker-based system for power grid visualization with AI-powered ChatGrid functionality using CSV data sources.

## Architecture
- **Frontend**: React application with ChatGrid AI integration
- **Backend**: Flask API with CSV-based data processing and AI chat capabilities  
- **Data Storage**: CSV files for power grid data (no database required)
- **Infrastructure**: Can be deployed on any Docker-compatible server

## Production Environment Details

### Server Directory Structure
```
/home/user/
├── westmap/                    # Main westmap project directory
│   ├── frontend/               # Frontend configuration
│   │   └── .env               # Frontend environment variables
│   ├── backend/               # Backend configuration  
│   │   └── .env               # Backend environment variables
│   └── data/                  # CSV data files
│       ├── bus.csv
│       ├── generation.csv
│       ├── transmission_line.csv
│       ├── counties.csv
│       ├── us states.csv
│       ├── WECC_BA_*.csv
│       └── *.json files
```

### Docker Images (Docker Hub: gradscaler organization)
1. **Backend**: `gradscaler/westmap-backend:V1.0.7`
   - Flask API with ChatGrid AI integration
   - CSV data processing with Pandas and Agno framework
   - OpenAI integration for natural language queries

2. **Frontend**: `gradscaler/westmap-frontend:V1.0.7`
   - React application with admin authentication system
   - Firebase integration for user management
   - ChatGrid widget for AI interactions

### Port Configuration
- **Frontend**: Port 8080 (public access)
- **Backend**: Port 5000 (public access for API)

**Note**: For production environments, consider using a reverse proxy (nginx) for HTTPS and routing.

## Complete Deployment Process

### 1. Server Access
```bash
# Connect to production server
ssh -i {path_to_your_ssh_key} {user}@{server_ip}
```

### 2. Initial Server Setup
```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install Docker
sudo apt install docker.io -y
sudo systemctl start docker
sudo systemctl enable docker
sudo usermod -aG docker $USER

# Logout and login again for group changes to take effect
```

### 3. Environment Configuration
Create the required directory structure and environment files:

```bash
# Create directory structure
mkdir -p ~/westmap/frontend
mkdir -p ~/westmap/backend
mkdir -p ~/westmap/data

# Download or copy CSV data files to ~/westmap/data/
# Ensure the following files are present:
# - bus.csv
# - generation.csv
# - transmission_line.csv
# - counties.csv
# - us states.csv
# - WECC_BA_*.csv files
```

#### Frontend Environment (.env)
```bash
# ~/westmap/frontend/.env
OPENAI_API_KEY=your_openai_key_here
NODE_ENV=production
REACT_APP_BACKEND_URL=http://localhost:5000
FIREBASE_API_KEY=your_firebase_api_key
FIREBASE_AUTH_DOMAIN=your_firebase_domain
FIREBASE_PROJECT_ID=your_firebase_project_id
FIREBASE_STORAGE_BUCKET=your_firebase_bucket
FIREBASE_MESSAGING_SENDER_ID=your_sender_id
FIREBASE_APP_ID=your_app_id
FIREBASE_MEASUREMENT_ID=your_measurement_id
```

#### Backend Environment (.env)
```bash
# ~/westmap/backend/.env
OPENAI_API_KEY=your_openai_key_here
FLASK_ENV=production
FLASK_DEBUG=False
```

### 4. Deploy Backend Container
```bash
# Go to backend directory and use .env file
cd ~/westmap/backend

# Create Docker network
docker network create gradscaler-network || true

# Pull and run backend container
docker pull gradscaler/westmap-backend:V1.0.7
docker run -d \
  --name westmap-backend \
  --network gradscaler-network \
  --env-file .env \
  -v ~/westmap/data:/app/data \
  -p 5000:5000 \
  gradscaler/westmap-backend:V1.0.7
```

### 5. Deploy Frontend Container
```bash
# Go to frontend directory and use .env file
cd ~/westmap/frontend

# Pull and run frontend container
docker pull gradscaler/westmap-frontend:V1.0.7
docker run -d \
  --name westmap-frontend \
  --network gradscaler-network \
  --env-file .env \
  -p 8080:8080 \
  gradscaler/westmap-frontend:V1.0.7
```

### 6. Verify Deployment
```bash
# Check container status
docker ps

# Test backend health
curl http://localhost:5000/health

# Test frontend access
curl http://localhost:8080

# Test ChatGrid API
curl -X POST -H "Content-Type: application/json" \
  -d '{"inputText":"How many power plants are there?"}' \
  http://localhost:5000/data
```

## Optional: Nginx Reverse Proxy Setup

For production environments with HTTPS and better routing:

### Install and Configure Nginx
```bash
# Install nginx
sudo apt install nginx -y
sudo systemctl enable nginx

# Create nginx configuration
sudo nano /etc/nginx/sites-available/westmap
```

### Nginx Configuration
```nginx
server {
    listen 80;
    server_name your-domain.com;

    # Frontend routing
    location / {
        proxy_pass http://localhost:8080;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Backend API routing
    location /api/ {
        proxy_pass http://localhost:5000/;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # CORS headers
        add_header Access-Control-Allow-Origin *;
        add_header Access-Control-Allow-Methods "GET, POST, OPTIONS";
        add_header Access-Control-Allow-Headers "Content-Type, Authorization";
    }

    # Health check
    location /health {
        proxy_pass http://localhost:5000/health;
    }
}
```

### Enable Nginx Configuration
```bash
# Enable the configuration
sudo ln -s /etc/nginx/sites-available/westmap /etc/nginx/sites-enabled/

# Remove default configuration
sudo rm /etc/nginx/sites-enabled/default

# Test configuration
sudo nginx -t

# Restart nginx
sudo systemctl restart nginx
```

## Container Update Process

### Update Backend Container
```bash
# Stop and remove old container
docker stop westmap-backend && docker rm westmap-backend

# Go to backend directory
cd ~/westmap/backend

# Pull updated image
docker pull gradscaler/westmap-backend:V1.0.7

# Start new container
docker run -d \
  --name westmap-backend \
  --network gradscaler-network \
  --env-file .env \
  -v ~/westmap/data:/app/data \
  -p 5000:5000 \
  gradscaler/westmap-backend:V1.0.7
```

### Update Frontend Container
```bash
# Stop and remove old container
docker stop westmap-frontend && docker rm westmap-frontend

# Go to frontend directory
cd ~/westmap/frontend

# Pull updated image
docker pull gradscaler/westmap-frontend:V1.0.7

# Start new container
docker run -d \
  --name westmap-frontend \
  --network gradscaler-network \
  --env-file .env \
  -p 8080:8080 \
  gradscaler/westmap-frontend:V1.0.7
```

## Service URLs & Access Points

### Public URLs (adjust domain as needed)
- **Main Application**: http://your-server:8080
- **Admin Dashboard**: http://your-server:8080/?admin=true
- **ChatGrid AI Widget**: Available in the main application interface
- **Backend API**: http://your-server:5000
- **Health Check**: http://your-server:5000/health

### Admin Features
- **User Registration Approval**: All new users require admin approval
- **Comprehensive User Management**: Approve, reject, promote, and delete users
- **Advanced Filtering**: Filter users by status (Pending, Approved, Admin)
- **Sorting Capabilities**: Sort by name, email, status, or registration date
- **Real-time Updates**: Dashboard updates automatically with Firestore
- **Action Confirmations**: Safety dialogs for destructive actions

## Data Information

### CSV Data Sources
The application processes the following CSV files:
- **bus.csv**: Electrical bus/substation records with voltage levels and coordinates
- **generation.csv**: Power generation facilities with capacity, type, and location data
- **transmission_line.csv**: Transmission line data with capacity and power flow information
- **counties.csv**: US county boundaries and demographic information
- **us states.csv**: US state boundaries and identifiers
- **WECC_BA_*.csv**: Western Electricity Coordinating Council balancing authority data

### Data Updates
To update data in production:
1. Upload new CSV files to `~/westmap/data/`
2. Restart the backend container:
   ```bash
   docker restart westmap-backend
   ```
3. New data is immediately available for ChatGrid queries

## Monitoring & Troubleshooting

### Check Container Status
```bash
# View all containers
docker ps

# Check specific container logs
docker logs westmap-frontend --tail 50
docker logs westmap-backend --tail 50
```

### Test API Endpoints
```bash
# Test frontend (should return HTML)
curl http://localhost:8080

# Test backend health
curl http://localhost:5000/health

# Test ChatGrid API
curl -X POST -H "Content-Type: application/json" \
  -d '{"inputText":"How many power plants are there?"}' \
  http://localhost:5000/data
```

### Common Issues & Solutions

#### 1. ChatGrid Not Responding
- **Symptom**: ChatGrid widget shows errors or no response
- **Solution**: Check backend logs and ensure OpenAI API key is valid
- **Check**: `docker logs westmap-backend --tail 50`

#### 2. Frontend Connection Issues
- **Symptom**: Frontend cannot connect to backend API
- **Solution**: Verify both containers are on the same network
- **Check**: `docker network inspect gradscaler-network`

#### 3. CSV Data Not Loading
- **Symptom**: Backend errors about missing CSV files
- **Solution**: Ensure CSV files are properly mounted in backend container
- **Check**: `docker exec westmap-backend ls -la /app/data`

#### 4. Authentication Issues
- **Symptom**: Admin dashboard not accessible or Firebase errors
- **Solution**: Verify Firebase configuration in frontend .env file
- **Check**: Ensure all Firebase environment variables are set correctly

## Emergency Procedures

### Restart All Services
```bash
# Restart both containers
docker restart westmap-backend westmap-frontend

# If using nginx
sudo systemctl restart nginx
```

### Complete System Reset
```bash
# Stop all containers
docker stop westmap-frontend westmap-backend
docker rm westmap-frontend westmap-backend

# Redeploy (follow deployment steps 4-5)
```

### Data Backup
```bash
# Backup CSV data
tar -czf westmap_data_backup_$(date +%Y%m%d).tar.gz ~/westmap/data/

# Backup environment files
tar -czf westmap_config_backup_$(date +%Y%m%d).tar.gz ~/westmap/frontend/.env ~/westmap/backend/.env
```

## Performance Optimization

### For Large CSV Files
If working with large CSV datasets:

1. **Memory Optimization**: Increase Docker container memory limits
2. **Data Chunking**: Process large files in chunks using Pandas
3. **Caching**: Implement data caching in the backend for frequently accessed queries
4. **Compression**: Use compressed CSV files (gzip) to reduce disk space

### For High Traffic
For production environments with high traffic:

1. **Load Balancing**: Deploy multiple backend containers behind a load balancer
2. **Caching**: Implement Redis for API response caching
3. **CDN**: Use a CDN for static frontend assets
4. **Health Monitoring**: Implement health checks and monitoring

## Version History & Changes

### V1.0.7 (CSV-based Architecture)
- **Removed Database Dependency**: Complete migration from PostgreSQL to CSV files
- **Simplified Deployment**: 2-tier architecture (frontend + backend)
- **Enhanced Data Processing**: Agno framework for AI-powered CSV analysis
- **Improved Performance**: Faster startup and reduced resource requirements
- **Streamlined Maintenance**: No database migrations or complex setup required

### Previous Versions (Deprecated)
- V1.0.2: 3-tier architecture with PostgreSQL database
- V1.0.1: Initial production deployment with database

---

**Last Updated**: August 28, 2025  
**Architecture**: CSV-based 2-tier system with simplified deployment
**Deployment Status**: Production-ready with minimal dependencies
