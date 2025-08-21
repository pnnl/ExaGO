-- Load CSV data manually
\copy generation FROM '/var/lib/postgresql/csv-data/generation.csv' DELIMITER ',' CSV HEADER;
\copy bus FROM '/var/lib/postgresql/csv-data/bus.csv' DELIMITER ',' CSV HEADER;
\copy transmission_line FROM '/var/lib/postgresql/csv-data/transmission_line.csv' DELIMITER ',' CSV HEADER;
\copy counties FROM '/var/lib/postgresql/csv-data/counties.csv' DELIMITER ',' CSV HEADER;
\copy us_states FROM '/var/lib/postgresql/csv-data/us states.csv' DELIMITER ',' CSV HEADER;

-- Show imported data count
SELECT 'generation' as table_name, count(*) as rows FROM generation
UNION ALL
SELECT 'bus' as table_name, count(*) as rows FROM bus
UNION ALL
SELECT 'transmission_line' as table_name, count(*) as rows FROM transmission_line
UNION ALL
SELECT 'counties' as table_name, count(*) as rows FROM counties
UNION ALL
SELECT 'us_states' as table_name, count(*) as rows FROM us_states;
