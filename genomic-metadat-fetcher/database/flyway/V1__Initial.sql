CREATE SCHEMA IF NOT EXISTS ena;

CREATE TABLE ena.runs (
    id SERIAL PRIMARY KEY,
    run_accession VARCHAR(20) UNIQUE NOT NULL,
    instrument_platform VARCHAR(100) NOT NULL,
    instrument_model VARCHAR(100) NOT NULL,
    organism_name VARCHAR(100) NOT NULL,
    last_updated TIMESTAMP NOT NULL
);

CREATE TABLE ena.file_locations (
    id SERIAL PRIMARY KEY,
    run_id INTEGER REFERENCES ena.runs(id) ON DELETE CASCADE,
    file_path TEXT NOT NULL,
    is_local BOOLEAN NOT NULL DEFAULT FALSE,
    UNIQUE (run_id, file_path)
);

CREATE INDEX idx_runs_platform ON ena.runs(instrument_platform);
CREATE INDEX idx_runs_model ON ena.runs(instrument_model);
CREATE INDEX idx_runs_organism ON ena.runs(organism_name);
CREATE INDEX idx_runs_updated_date ON ena.runs(last_updated);
CREATE INDEX idx_file_locations_run_id ON ena.file_locations(run_id);
