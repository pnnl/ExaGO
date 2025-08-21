# WestMap Production Deployment Guide

## Overview
This guide documents the complete production deployment of the WestMap application - a 3-tier Docker-based system for power grid visualization with AI-powered ChatGrid functionality, deployed with nginx reverse proxy to resolve CORS issues.

## Architecture
- **Frontend**: React application with ChatGrid AI integration
- **Backend**: Flask API with PostgreSQL integration and AI chat capabilities  
- **Database**: PostgreSQL with PostGIS extensions for spatial data
- **Reverse Proxy**: Nginx for routing and CORS handling
- **Infrastructure**: Ubuntu 22.04.5 LTS server on Azure

## Production Environment Details

### Server Directory Structure
```
/home/azureuser/
├── westmap/                    # Main westmap project directory
│   ├── frontend/               # Frontend configuration
│   │   └── .env               # Frontend environment variables
│   ├── backend/               # Backend configuration  
│   │   └── .env               # Backend environment variables
│   └── database/              # Database configuration
│       └── .env               # Database environment variables
├── data/                      # Data files for ChatGrid
│   ├── bus.csv
│   ├── generation.csv
│   ├── transmission_line.csv
│   ├── counties.csv
│   ├── us states.csv
│   └── *.json files
```

### Docker Images (Docker Hub: gradscaler organization)
1. **Database**: `gradscaler/westmap-database:V1.0.2`
   - PostgreSQL 15.4 with PostGIS extensions
   - Preloaded with power grid data (849 generators, 4762 buses, 12706 transmission lines)

2. **Backend**: `gradscaler/westmap-backend:V1.0.2`
   - Flask API with ChatGrid AI integration
   - Database connectivity and spatial queries

3. **Frontend**: `gradscaler/westmap-frontend:V1.0.2`
   - React application with admin authentication system
   - Firebase integration for user management
   - Enhanced admin dashboard with user approval workflow
   - Modified to use `/api/data` instead of direct backend URLs

### Port Configuration
- **Nginx**: Port 80 (public access, reverse proxy)
- **Frontend**: Port 8080 (internal, proxied via nginx)
- **Backend**: Port 5000 (internal, proxied via nginx)
- **Database**: Port 5432 (internal only)

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

# Install Nginx
sudo apt install nginx -y
sudo systemctl enable nginx
```

### 2.1. Environment Configuration
Each component requires proper .env files in their respective directories:

**Database (.env)**:
```bash
# ~/westmap/database/.env
POSTGRES_DB=exago_viz
POSTGRES_USER=postgres
POSTGRES_PASSWORD=exago2024
PGDATA=/var/lib/postgresql/data
```

**Backend (.env)**:
```bash
# ~/westmap/backend/.env
OPENAI_API_KEY=your_openai_key_here
SQL_KEY=exago2024
DATABASE_NAME=exago_viz
DB_HOST=localhost
DB_PORT=5432
DB_USER=postgres
FLASK_ENV=development
FLASK_DEBUG=True
```

**Frontend (.env)**:
```bash
# ~/westmap/frontend/.env
OPENAI_API_KEY=your_openai_key_here
SQL_KEY=exago2024
DATABASE_NAME=exago_viz
DB_HOST=localhost
DB_PORT=5432
DB_USER=postgres
FIREBASE_API_KEY=your_firebase_api_key
FIREBASE_AUTH_DOMAIN=your_firebase_domain
FIREBASE_PROJECT_ID=your_firebase_project_id
FIREBASE_STORAGE_BUCKET=your_firebase_bucket
FIREBASE_MESSAGING_SENDER_ID=your_sender_id
FIREBASE_APP_ID=your_app_id
FIREBASE_MEASUREMENT_ID=your_measurement_id
REACT_APP_BACKEND_URL=http://localhost:5000
NODE_ENV=development
```

### 3. Deploy Database Container
```bash
# Go to database directory and use .env file
cd ~/westmap/database
docker pull gradscaler/westmap-database:V1.0.1
docker volume create westmap-postgres-data
docker network create gradscaler-network || true
docker run -d \
  --name westmap-database \
  --network gradscaler-network \
  -p 5432:5432 \
  -v westmap-postgres-data:/var/lib/postgresql/data \
  --env-file .env \
  --restart unless-stopped \
  gradscaler/westmap-database:V1.0.1
```

### 4. Deploy Backend Container
```bash
# Go to backend directory and use .env file
cd ~/westmap/backend
docker pull gradscaler/westmap-backend:V1.0.2
docker run -d \
  --name westmap-backend \
  --network gradscaler-network \
  -p 5000:5000 \
  -v ~/data:/app/data \
  --env-file .env \
  --restart unless-stopped \
  gradscaler/westmap-backend:V1.0.2
```

### 5. Deploy Frontend Container
```bash
# Go to frontend directory and use .env file
cd ~/westmap/frontend
docker pull gradscaler/westmap-frontend:V1.0.2
docker run -d \
  --name westmap-frontend \
  --network gradscaler-network \
  -p 8080:8080 \
  --env-file .env \
  --restart unless-stopped \
  gradscaler/westmap-frontend:V1.0.2
```

### 6. Configure Nginx Reverse Proxy

#### Create Nginx Configuration
```bash
sudo nano /etc/nginx/sites-available/westmap
```

#### Nginx Configuration Content
```nginx
server {
    listen 80;
    server_name 20.198.225.186;

    # Frontend - serve React app
    location / {
        proxy_pass http://localhost:8080;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_cache_bypass $http_upgrade;
        
        # CORS headers
        add_header Access-Control-Allow-Origin *;
        add_header Access-Control-Allow-Methods "GET, POST, OPTIONS";
        add_header Access-Control-Allow-Headers "Content-Type, Authorization";
    }

    # Backend API - proxy to Flask backend
    location /api/ {
        # Remove /api prefix when forwarding to backend
        rewrite ^/api/(.*)$ /$1 break;
        
        proxy_pass http://localhost:5000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_cache_bypass $http_upgrade;
        
        # CORS headers
        add_header Access-Control-Allow-Origin *;
        add_header Access-Control-Allow-Methods "GET, POST, OPTIONS";
        add_header Access-Control-Allow-Headers "Content-Type, Authorization";
    }

    # Health check endpoint
    location /health {
        proxy_pass http://localhost:5000/health;
    }
}
```

#### Enable Nginx Configuration
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

## Frontend Code Changes for Nginx Compatibility

### Critical Frontend Modification
The key change that resolved CORS issues was updating `app.js` to use relative API paths:

**File**: `viz/app.js` (around line 696)

**BEFORE** (causing CORS errors):
```javascript
const backendUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:5000';
fetch(`${backendUrl}/data`, {
    method: 'POST',
    body: JSON.stringify(postData),
})
```

**AFTER** (working with nginx reverse proxy):
```javascript
// Use relative API path through nginx reverse proxy
const apiPath = '/api/data';

fetch(apiPath, {
    method: 'POST',
    body: JSON.stringify(postData),
})
```

This change allows the frontend to make requests to `/api/data` which nginx routes to the backend at `localhost:5000/data`, eliminating CORS issues by serving both frontend and backend from the same origin.

## Container Update Process

### Container Update Process

### Update Frontend Container
```bash
# Stop and remove old container
docker stop westmap-frontend && docker rm westmap-frontend

# Go to frontend directory
cd ~/westmap/frontend

# Pull updated image
docker pull gradscaler/westmap-frontend:V1.0.2

# Start new container with proper configuration
docker run -d \
  --name westmap-frontend \
  --network gradscaler-network \
  -p 8080:8080 \
  --env-file .env \
  --restart unless-stopped \
  gradscaler/westmap-frontend:V1.0.2
```

### Update Backend Container
```bash
# Stop and remove old container
docker stop westmap-backend && docker rm westmap-backend

# Go to backend directory
cd ~/westmap/backend

# Pull updated image  
docker pull gradscaler/westmap-backend:V1.0.2

# Start new container with proper configuration
docker run -d \
  --name westmap-backend \
  --network gradscaler-network \
  -p 5000:5000 \
  -v ~/data:/app/data \
  --env-file .env \
  --restart unless-stopped \
  gradscaler/westmap-backend:V1.0.2
```

## Service URLs & Access Points

### Public URLs
- **Main Application**: http://20.198.225.186
- **Admin Dashboard**: http://20.198.225.186/?admin=true
- **ChatGrid AI Widget**: Available in the main application interface
- **Health Check**: http://20.198.225.186/health

### Admin Features (New in V1.0.2)
- **User Registration Approval**: All new users require admin approval
- **Comprehensive User Management**: Approve, reject, promote, and delete users
- **Advanced Filtering**: Filter users by status (Pending, Approved, Admin)
- **Sorting Capabilities**: Sort by name, email, status, or registration date
- **Real-time Updates**: Dashboard updates automatically with Firestore
- **Action Confirmations**: Safety dialogs for destructive actions

### Direct Access (for troubleshooting)
- **Frontend Direct**: http://20.198.225.186:8080
- **Backend API Direct**: http://20.198.225.186:5000

## Database Information

### Connection Details
- **Host**: 20.198.225.186:5432
- **Database**: westmap
- **User**: westmap_user
- **Password**: westmap_password

### Data Statistics
- **Generators**: 849 power generation facilities
- **Buses**: 4,762 electrical connection points
- **Transmission Lines**: 12,706 power transmission lines
- **Total Records**: 18,317+ spatial and operational records

## Monitoring & Troubleshooting

### Check Container Status
```bash
# View all containers
docker ps

# Check specific container logs
docker logs westmap-frontend --tail 50
docker logs westmap-backend --tail 50
docker logs westmap-database --tail 50
```

### Check Nginx Status
```bash
# Test nginx configuration
sudo nginx -t

# Check nginx status
sudo systemctl status nginx

# View nginx logs
sudo tail -f /var/log/nginx/error.log
sudo tail -f /var/log/nginx/access.log
```

### Test API Endpoints
```bash
# Test frontend (should return HTML)
curl http://20.198.225.186

# Test backend API through nginx proxy
curl -X POST -H "Content-Type: application/json" \
  -d '{"inputText":"How many power plants are there?"}' \
  http://20.198.225.186/api/data

# Test backend health check
curl http://20.198.225.186/health
```

### Common Issues & Solutions

#### 1. CORS Errors
- **Symptom**: Browser console shows CORS policy errors
- **Solution**: Ensure nginx is running and frontend uses `/api/` paths
- **Check**: Verify frontend container uses V1.0.2 with updated API paths

#### 2. ChatGrid Widget Not Working
- **Symptom**: "net::ERR_BLOCKED_BY_CLIENT" or connection refused
- **Solution**: Ensure requests go through nginx proxy at `/api/data`
- **Check**: Browser dev tools should show requests to same origin

#### 3. Database Connection Issues
- **Symptom**: Backend logs show database connection errors
- **Solution**: Verify database container is running and accessible
- **Check**: `docker exec westmap-database pg_isready -U westmap_user`

#### 4. Nginx Configuration Errors
- **Symptom**: 502 Bad Gateway or service unavailable
- **Solution**: Test nginx config and restart service
- **Check**: `sudo nginx -t && sudo systemctl restart nginx`

## Emergency Procedures

### Restart All Services
```bash
# Restart all containers
docker restart westmap-database westmap-backend westmap-frontend

# Restart nginx
sudo systemctl restart nginx
```

### Complete System Reset
```bash
# Stop all containers
docker stop westmap-frontend westmap-backend westmap-database
docker rm westmap-frontend westmap-backend westmap-database

# Redeploy (follow deployment steps 3-5)
# Database will preserve data in Docker volumes
```

### Database Backup
```bash
# Create backup
docker exec westmap-database pg_dump -U westmap_user westmap > westmap_backup_$(date +%Y%m%d).sql

# Restore backup (if needed)
docker exec -i westmap-database psql -U westmap_user westmap < westmap_backup_YYYYMMDD.sql
```

## Version History & Changes

### V1.0.2 (Admin Authentication & User Management)
- **Firebase Integration**: Complete user authentication system
- **Admin Approval Workflow**: All users require admin approval before access
- **Advanced Admin Dashboard**: Comprehensive user management interface
- **User Management Features**:
  - Approve/reject new user registrations
  - Promote users to admin status
  - Delete user accounts with confirmation
  - Real-time user status updates
  - Advanced filtering and sorting capabilities
- **Enhanced Security**: Role-based access control
- **Improved UI**: Material-UI components with professional styling

### V1.0.1 (Initial Production Deployment)
- Basic 3-tier Docker stack
- Direct frontend-to-backend communication
- CORS issues encountered

### V1.0.2 (Nginx Reverse Proxy Solution)
- Added nginx reverse proxy configuration
- Updated frontend to use relative API paths (`/api/data`)
- Resolved CORS issues by serving both frontend and backend from same origin
- Improved production architecture with single entry point

## Security Considerations

1. **Firewall**: Only ports 80, 22, and 8080 should be publicly accessible
2. **Database**: Port 5432 should only be accessible internally
3. **SSH Key**: Secure `gradscalerprod_key.pem` with proper permissions (600)
4. **SSL/TLS**: Consider adding HTTPS certificate for production
5. **Environment Variables**: Database credentials should be properly secured

## Performance Optimization

1. **Nginx Caching**: Static assets cached automatically
2. **Database Indexing**: Spatial indexes optimize geographic queries
3. **Container Resources**: Monitor CPU/memory usage and adjust limits
4. **Network Latency**: All containers on same server minimize latency

## Future Enhancements

1. **Load Balancing**: Scale backend with multiple instances
2. **HTTPS Setup**: Add SSL certificate and redirect HTTP to HTTPS
3. **Monitoring**: Implement logging aggregation and health monitoring
4. **Backup Automation**: Schedule automated database backups
5. **Container Orchestration**: Consider Docker Compose or Kubernetes

---

## Deployment Success Confirmation

✅ **Complete Production System Successfully Deployed**

**Architecture**: 3-tier containerized application with nginx reverse proxy
**Services**: PostgreSQL database, Flask backend, React frontend
**Data**: 18,317+ records of power grid infrastructure
**AI Integration**: ChatGrid natural language query interface
**CORS Resolution**: Nginx reverse proxy serving all services from single origin
**High Availability**: Auto-restart policies and persistent data storage

**Deployment Date**: August 9, 2025
**Production URL**: http://20.198.225.186
**Admin Dashboard**: http://20.198.225.186/?admin=true
**Status**: Fully operational with ChatGrid AI and Admin Authentication System

---

*This deployment guide documents a production-ready system that successfully resolved CORS issues through nginx reverse proxy architecture, enabling seamless frontend-backend communication for the ChatGrid AI functionality.*
