-- Initialize PostGIS extension
CREATE EXTENSION IF NOT EXISTS postgis;

-- Import generation data
CREATE TABLE IF NOT EXISTS generation (
    NAME TEXT,
    nbus INTEGER,
    Pg FLOAT,
    Pcap FLOAT,
    gen_fuel TEXT,
    ngen INTEGER,
    KV TEXT,
    Pgcoal FLOAT,
    Pghydro FLOAT,
    Pgnuclear FLOAT,
    Pgng FLOAT,
    Pgsolar FLOAT,
    Pgwind FLOAT,
    Pgother FLOAT,
    Pgcoalcap FLOAT,
    Pghydrocap FLOAT,
    Pgnuclearcap FLOAT,
    Pgngcap FLOAT,
    Pgsolarcap FLOAT,
    Pgwindcap FLOAT,
    Pgothercap FLOAT,
    lat FLOAT,
    lon FLOAT,
    county TEXT,
    state TEXT
);

-- Import bus data
CREATE TABLE IF NOT EXISTS bus (
    BUS_I INTEGER,
    BUS_TYPE INTEGER,
    PD FLOAT,
    QD FLOAT,
    GS FLOAT,
    BS FLOAT,
    BUS_AREA INTEGER,
    VM FLOAT,
    VA FLOAT,
    BASE_KV FLOAT,
    ZONE INTEGER,
    VMAX FLOAT,
    VMIN FLOAT,
    LAM_P FLOAT,
    LAM_Q FLOAT,
    MU_VMAX FLOAT,
    MU_VMIN FLOAT,
    coordx FLOAT,
    coordy FLOAT,
    county TEXT,
    state TEXT
);

-- Import transmission line data
CREATE TABLE IF NOT EXISTS transmission_line (
    FBUS INTEGER,
    TBUS INTEGER,
    R FLOAT,
    X FLOAT,
    B FLOAT,
    RATEA FLOAT,
    RATEB FLOAT,
    RATEC FLOAT,
    TAP FLOAT,
    SHIFT FLOAT,
    BR_STATUS INTEGER,
    ANGMIN FLOAT,
    ANGMAX FLOAT,
    PF FLOAT,
    QF FLOAT,
    PT FLOAT,
    QT FLOAT,
    MU_SF FLOAT,
    MU_ST FLOAT,
    MU_ANGMIN FLOAT,
    MU_ANGMAX FLOAT,
    coordx_f FLOAT,
    coordy_f FLOAT,
    coordx_t FLOAT,
    coordy_t FLOAT
);

-- Import counties data
CREATE TABLE IF NOT EXISTS counties (
    NAME TEXT,
    STATEFP TEXT,
    COUNTYFP TEXT,
    COUNTYNS TEXT,
    GEOID TEXT,
    ALAND BIGINT,
    AWATER BIGINT,
    geometry TEXT
);

-- Import US states data
CREATE TABLE IF NOT EXISTS us_states (
    NAME TEXT,
    STATEFP TEXT,
    STATENS TEXT,
    AFFGEOID TEXT,
    GEOID TEXT,
    STUSPS TEXT,
    ALAND BIGINT,
    AWATER BIGINT,
    geometry TEXT
);
