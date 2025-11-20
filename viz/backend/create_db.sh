#!/usr/bin/env bash
set -euo pipefail

# ---- Defaults (can be overridden by flags) ----
PGHOST="${PGHOST:-localhost}"
PGPORT="${PGPORT:-5432}"
PGUSER="${PGUSER:-$USER}"
PGPASSWORD="${PGPASSWORD:-}"           # export PGPASSWORD=... if needed
DB_NAME=""
SCHEMA_NAME="public"
TABLE_NAME=""
DELIM=","
QUOTE='"'
NULLSTR=""
ENCODING="UTF8"
SCHEMA_SQL=""
DROP_IF_EXISTS="false"
TRUNCATE_FIRST="false"

usage() {
  cat <<EOF
Usage: $(basename "$0") --db NAME --table NAME --csv /path/file.csv [options]

Required:
  --db NAME              Database name to create/use
  --table NAME           Target table name (without schema; uses --schema)

Optional:
  --schema NAME          Schema name (default: public)
  --schema-sql PATH      Execute this SQL file to define table (overrides auto-create)
  --drop                 Drop table if it already exists (dangerous)
  --truncate             TRUNCATE the table before loading (if it exists)
  --host HOST            PG host (default: ${PGHOST})
  --port PORT            PG port (default: ${PGPORT})
  --user USER            PG user (default: ${PGUSER})
  --delim CHAR           CSV delimiter (default: ",")
  --quote CHAR           CSV quote (default: '"')
  --null STR             CSV NULL string (default: empty -> '')
  --encoding ENC         Client encoding (default: UTF8)

Environment:
  PGPASSWORD can be exported for password auth.

Examples:
  PGPASSWORD=secret \\
  $(basename "$0") --db mydb --table sales --csv ./sales.csv --truncate

  $(basename "$0") --db mydb --table sales \\
    --schema analytics --csv ./sales.csv --delim ';'
EOF
  exit 1
}

# ---- Parse args ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db) DB_NAME="$2"; shift 2 ;;
    --csv) CSV_PATH="$2"; shift 2 ;;
    --schema) SCHEMA_NAME="$2"; shift 2 ;;
    --schema-sql) SCHEMA_SQL="$2"; shift 2 ;;
    --drop) DROP_IF_EXISTS="true"; shift 1 ;;
    --truncate) TRUNCATE_FIRST="true"; shift 1 ;;
    --host) PGHOST="$2"; shift 2 ;;
    --port) PGPORT="$2"; shift 2 ;;
    --user) PGUSER="$2"; shift 2 ;;
    --delim) DELIM="$2"; shift 2 ;;
    --quote) QUOTE="$2"; shift 2 ;;
    --null) NULLSTR="$2"; shift 2 ;;
    --encoding) ENCODING="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

[[ -z "${DB_NAME}" ]] && usage

export PGHOST PGPORT PGUSER PGPASSWORD

# Small helper to run psql with standard flags
psqlx() {
  psql --no-psqlrc --set ON_ERROR_STOP=1 "$@"
}

echo ">> Checking server connectivity..."
psqlx -d postgres -tAc 'SELECT 1;' >/dev/null

echo ">> Ensuring database '${DB_NAME}' exists..."
DB_EXISTS=$(psqlx -d postgres -tAc "SELECT 1 FROM pg_database WHERE datname='${DB_NAME}'" || true)
if [[ -z "${DB_EXISTS}" ]]; then
  psqlx -d postgres -c "CREATE DATABASE ${DB_NAME};"
  echo "   Created database ${DB_NAME}"
else
  echo "   Database already exists."
fi

echo ">> Ensuring schema '${SCHEMA_NAME}' exists..."
psqlx -d "${DB_NAME}" -c "CREATE SCHEMA IF NOT EXISTS ${SCHEMA_NAME};"




FQTN="${SCHEMA_NAME}.${TABLE_NAME}"

if [[ "${DROP_IF_EXISTS}" == "true" ]]; then
  echo ">> Dropping tables"
  psqlx -d "${DB_NAME}" -c "DROP TABLE IF EXISTS ${SCHEMA_NAME}.transmission_lines;"
  psqlx -d "${DB_NAME}" -c "DROP TABLE IF EXISTS ${SCHEMA_NAME}.buses;"
  psqlx -d "${DB_NAME}" -c "DROP TABLE IF EXISTS ${SCHEMA_NAME}.generators;"
  psqlx -d "${DB_NAME}" -c "DROP TABLE IF EXISTS ${SCHEMA_NAME}.states;"
  psqlx -d "${DB_NAME}" -c "DROP TABLE IF EXISTS ${SCHEMA_NAME}.counties;"
fi

# ---- Create table ----
echo ">> Creating table(s) via schema file: ${SCHEMA_SQL}"
[[ -f "${SCHEMA_SQL}" ]] || { echo "Schema SQL file not found: ${SCHEMA_SQL}"; exit 1; }
psqlx -d "${DB_NAME}" -f "${SCHEMA_SQL}"


# ---- Optionally TRUNCATE ----
if [[ "${TRUNCATE_FIRST}" == "true" ]]; then
  echo ">> Truncating tables"
  psqlx -d "${DB_NAME}" -c "TRUNCATE ${SCHEMA_NAME}.transmission_lines;"
  psqlx -d "${DB_NAME}" -c "TRUNCATE ${SCHEMA_NAME}.buses;"
  psqlx -d "${DB_NAME}" -c "TRUNCATE ${SCHEMA_NAME}.generators;"
  psqlx -d "${DB_NAME}" -c "TRUNCATE ${SCHEMA_NAME}.states;"
  psqlx -d "${DB_NAME}" -c "TRUNCATE ${SCHEMA_NAME}.counties;"
fi

# ---- Import CSV ----
echo ">> Importing CSV into ${FQTN}..."
# Ensure client encoding
psqlx -d "${DB_NAME}" -c "SET client_encoding TO '${ENCODING}';"

# Build COPY options
copy_opts="FORMAT csv, HEADER true, DELIMITER '${DELIM}', QUOTE '${QUOTE}'"
# NULL handling
if [[ -n "${NULLSTR}" ]]; then
  copy_opts="${copy_opts}, NULL '${NULLSTR}'"
fi

# Use \copy so the client reads the local file
psqlx -d "${DB_NAME}" -c "\COPY ${SCHEMA_NAME}.transmission_lines FROM 'transmission_line.csv' WITH (${copy_opts});"
psqlx -d "${DB_NAME}" -c "\COPY ${SCHEMA_NAME}.buses FROM 'bus.csv' WITH (${copy_opts});"
psqlx -d "${DB_NAME}" -c "\COPY ${SCHEMA_NAME}.generators FROM 'generation.csv' WITH (${copy_opts});"
psqlx -d "${DB_NAME}" -c "\COPY ${SCHEMA_NAME}.states FROM 'us_states.csv' WITH (${copy_opts});"
psqlx -d "${DB_NAME}" -c "\COPY ${SCHEMA_NAME}.counties FROM 'counties.csv' WITH (${copy_opts});"
