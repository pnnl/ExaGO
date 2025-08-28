# Variables 
USERNAME := gradscaler
IMAGE_NAME := westmap-backend
VERSION := V1.0.3
PORT := 5000
CONTAINER_NAME := westmap-backend
NETWORK := ${USERNAME}-network

# Image tag
TAG := $(USERNAME)/$(IMAGE_NAME):$(VERSION)

# Build the Docker image
build: 
	docker build -t $(TAG) .

# Create network (if not exists)
network: 
	docker network inspect $(NETWORK) >/dev/null 2>&1 || docker network create $(NETWORK)

# Run the Docker container with volume mount for data
run: network
	docker run -d -p $(PORT):$(PORT) --network $(NETWORK) --env-file .env \
		-v "$(PWD)/../data:/app/data" --name $(CONTAINER_NAME) $(TAG)

# Run with database setup (alternative command that runs setup first)
run-with-setup: network
	docker run -d -p $(PORT):$(PORT) --network $(NETWORK) --env-file .env \
		-v "$(PWD)/../data:/app/data" --name $(CONTAINER_NAME) \
		--env SETUP_DB=true $(TAG) sh -c "python setup_database.py && python server.py"

# Stop and remove container
clean:
	docker stop $(CONTAINER_NAME) || true
	docker rm $(CONTAINER_NAME) || true

# Push the Docker image to registry
push:
	docker push $(TAG)

# Show logs
logs:
	docker logs -f $(CONTAINER_NAME)

# Shell into container
shell:
	docker exec -it $(CONTAINER_NAME) /bin/sh

# Run database setup only
setup-db: network
	docker run --rm --network $(NETWORK) --env-file .env \
		-v "$(PWD)/../data:/app/data" $(TAG) python setup_database.py
