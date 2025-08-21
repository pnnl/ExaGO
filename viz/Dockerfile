# Use Alpine-based Node.js image
FROM node:16-alpine

# Set working directory
WORKDIR /app

# Copy package files
COPY package*.json ./

# Install dependencies with legacy peer deps flag
RUN npm install --legacy-peer-deps

# Copy application code
COPY . .

# Create data directory and copy data files
RUN mkdir -p /app/data
COPY data/ /app/data/

# Expose port (default webpack-dev-server port)
EXPOSE 8080

# Start the application on all interfaces for Docker
CMD ["npx", "webpack-dev-server", "--host", "0.0.0.0", "--port", "8080", "--progress", "--hot"]
