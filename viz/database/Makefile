# Westmap PostgreSQL Database Makefile
# This Makefile manages the PostgreSQL database container for ChatGrid functionality

# Configuration
USERNAME = gradscaler
IMAGE_NAME = westmap-database
VERSION = V1.0.2
TAG = $(USERNAME)/$(IMAGE_NAME):$(VERSION)
CONTAINER_NAME = westmap-database
NETWORK = gradscaler-network
DATA_VOLUME = westmap-postgres-data

# Create network if it doesn't exist
.PHONY: network
network:
	@docker network ls | grep -q $(NETWORK) || docker network create $(NETWORK)
	@echo "✓ Network $(NETWORK) is ready"

# Build the database image
.PHONY: build
build: network
	@echo "Building PostgreSQL database image..."
	@cp -r ../data ./data
	@docker build -t $(TAG) .
	@rm -rf ./data
	@echo "✓ Database image built: $(TAG)"

# Create persistent volume for database data
.PHONY: volume
volume:
	@docker volume ls | grep -q $(DATA_VOLUME) || docker volume create $(DATA_VOLUME)
	@echo "✓ Volume $(DATA_VOLUME) is ready"

# Run the database container
.PHONY: run
run: network volume
	@echo "Starting PostgreSQL database container..."
	@docker run -d \
		--name $(CONTAINER_NAME) \
		--network $(NETWORK) \
		-p 5432:5432 \
		-v $(DATA_VOLUME):/var/lib/postgresql/data \
		--env-file .env \
		$(TAG)
	@echo "✓ Database container started: $(CONTAINER_NAME)"
	@echo "  - Database: loaded from .env"
	@echo "  - Port: 5432"
	@echo "  - User: loaded from .env"
	@echo "  - Password: loaded from .env"

# Stop and clean the container
.PHONY: clean
clean:
	@echo "Stopping and removing database container..."
	@docker stop $(CONTAINER_NAME) 2>/dev/null || true
	@docker rm $(CONTAINER_NAME) 2>/dev/null || true
	@echo "✓ Database container cleaned"

# Remove everything including the volume (WARNING: This deletes all data!)
.PHONY: clean-all
clean-all: clean
	@echo "⚠️  WARNING: This will delete all database data!"
	@read -p "Are you sure? (y/N): " confirm && [ "$$confirm" = "y" ] || exit 1
	@docker volume rm $(DATA_VOLUME) 2>/dev/null || true
	@echo "✓ Database volume removed"

# Push the image to registry
.PHONY: push
push:
	@echo "Pushing database image to registry..."
	@docker push $(TAG)
	@echo "✓ Image pushed: $(TAG)"

# Show container logs
.PHONY: logs
logs:
	@echo "Database container logs:"
	@docker logs $(CONTAINER_NAME) --tail 50 -f

# Get shell access to the running container
.PHONY: shell
shell:
	@echo "Connecting to database container shell..."
	@docker exec -it $(CONTAINER_NAME) /bin/sh

# Connect to PostgreSQL using psql
.PHONY: psql
psql:
	@echo "Connecting to PostgreSQL database..."
	@docker exec -it $(CONTAINER_NAME) psql -U postgres -d exago_viz

# Show database status and tables
.PHONY: status
status:
	@echo "Database container status:"
	@docker ps --filter name=$(CONTAINER_NAME) --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
	@echo ""
	@echo "Database tables:"
	@docker exec $(CONTAINER_NAME) psql -U postgres -d exago_viz -c "\dt" 2>/dev/null || echo "  Database not ready or no tables found"

# Initialize database with sample data (run after container is started)
.PHONY: setup
setup:
	@echo "Setting up database with sample data..."
	@sleep 5  # Wait for database to be ready
	@docker exec $(CONTAINER_NAME) psql -U postgres -d exago_viz -c "SELECT import_csv_data();" || echo "Data import may have failed - check logs"
	@echo "✓ Database setup completed"

# Show help
.PHONY: help
help:
	@echo "Available commands:"
	@echo "  build     - Build the PostgreSQL database image"
	@echo "  run       - Start the database container"
	@echo "  clean     - Stop and remove the container"
	@echo "  clean-all - Remove container and delete all data (WARNING!)"
	@echo "  push      - Push image to registry"
	@echo "  logs      - Show container logs"
	@echo "  shell     - Get shell access to container"
	@echo "  psql      - Connect to PostgreSQL database"
	@echo "  status    - Show container and database status"
	@echo "  setup     - Initialize database with sample data"
	@echo "  help      - Show this help message"
