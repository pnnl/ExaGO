# Variables 
USERNAME := gradscaler
IMAGE_NAME := westmap-frontend
VERSION := V1.0.4
PORT := 8080
CONTAINER_NAME := westmap-frontend
NETWORK := ${USERNAME}-network

# Image tag
TAG := $(USERNAME)/$(IMAGE_NAME):$(VERSION)

# Build the Docker image
build: 
	docker build -t $(TAG) .

# Create network (if not exists)
network: 
	docker network inspect $(NETWORK) >/dev/null 2>&1 || docker network create $(NETWORK)

# Run the Docker container
run: network
	docker run -d -p $(PORT):$(PORT) --network $(NETWORK) --env-file .env --name $(CONTAINER_NAME) $(TAG)

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
