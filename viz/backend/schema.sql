CREATE TABLE IF NOT EXISTS transmission_lines (
    coordinates TEXT NOT NULL,
    flow_capacity NUMERIC NOT NULL,
    pf NUMERIC NOT NULL,
    qf NUMERIC NOT NULL,
    pt NUMERIC NOT NULL,
    qt NUMERIC NOT NULL,
    kilovolt NUMERIC NOT NULL,
    line_name TEXT NOT NULL,
    source TEXT NOT NULL,
    target TEXT NOT NULL,
    actual_flow NUMERIC NOT NULL,
    loading_percentage NUMERIC NOT NULL
);

CREATE TABLE IF NOT EXISTS buses (
    coordinates TEXT NOT NULL,
    bus_name TEXT NOT NULL,
    kv_levels NUMERIC[] NOT NULL,
    number_of_buses INTEGER NOT NULL,
    vm NUMERIC NOT NULL,
    county_name TEXT,
    state_name TEXT
);



CREATE TABLE IF NOT EXISTS generators (
    coordinates TEXT NOT NULL,        
    power_generated NUMERIC NOT NULL,
    power_capacity NUMERIC NOT NULL,
    remaining_capacity NUMERIC NOT NULL,
    kv_levels NUMERIC[] NOT NULL,
    color TEXT NOT NULL,
    generation_name TEXT NOT NULL,
    number_of_buses INTEGER NOT NULL,
    generation_type TEXT NOT NULL,
    county_name TEXT,
    state_name TEXT
);

CREATE TABLE IF NOT EXISTS states (
    coordinates TEXT,
    geoid TEXT NOT NULL,
    state_code INTEGER,
    state_name TEXT,
    census_area NUMERIC NOT NULL
);

CREATE TABLE counties (
    coordinates TEXT NOT NULL,
    geoid TEXT NOT NULL,     
    state_code INTEGER NOT NULL,    
    county_name TEXT,       
    census_area NUMERIC NOT NULL, 
    state_name TEXT           
);
