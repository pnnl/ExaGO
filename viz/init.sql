-- init.sql
DROP EXTENSION IF EXISTS postgis_tiger_geocoder;

CREATE TABLE IF NOT EXISTS bus (
    location geometry(Point,4326),
    bus_name character varying(255),
    kilovolt_levels double precision[],
    number_of_buses integer,
    vm double precision
);
CREATE TABLE IF NOT EXISTS counties (
    geom geometry(MultiPolygon,4326),
    geo_id character varying(50),
    fips integer,
    name character varying(255),
    censusarea double precision,
    state character varying(100)
);
CREATE TABLE IF NOT EXISTS generation (
    location geometry(Point,4326),
    power_generated double precision,
    power_capacity double precision,
    kv_levels double precision[],
    color character varying(50),
    generation_name character varying(255),
    number_of_buses integer,
    generation_type character varying(50)
);
CREATE TABLE IF NOT EXISTS transmission_line (
    geom geometry(LineString,4326),
    flow_capacity double precision,
    pf double precision,
    qf double precision,
    pt double precision,
    qt double precision,
    kilovolt double precision,
    line_name character varying,
    source character varying,
    target character varying,
    actual_flow double precision
);
CREATE TABLE IF NOT EXISTS us_states (
    geom geometry(MultiPolygon,4326),
    geoid character varying(50),
    state_code integer,
    state_name character varying(100),
    census_area double precision
);

COPY bus FROM '/mnt/bus.csv' DELIMITER ',' CSV HEADER;
COPY counties FROM '/mnt/counties.csv' DELIMITER ',' CSV HEADER;
COPY generation FROM '/mnt/generation.csv' DELIMITER ',' CSV HEADER;
COPY transmission_line FROM '/mnt/transmission_line.csv' DELIMITER ',' CSV HEADER;
COPY us_states FROM '/mnt/us states.csv' DELIMITER ',' CSV HEADER;
