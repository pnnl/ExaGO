# PostgreSQL with PostGIS for Westmap ChatGrid
FROM postgis/postgis:15-3.3-alpine

# Environment variables will be loaded from .env file during docker run
# No hardcoded secrets in Dockerfile

# Copy initialization script
COPY init-postgis.sql /docker-entrypoint-initdb.d/

# Copy data files for initialization
COPY data/ /var/lib/postgresql/csv-data/

# Expose PostgreSQL port
EXPOSE 5432

# Use the default postgres entrypoint
CMD ["postgres"]
