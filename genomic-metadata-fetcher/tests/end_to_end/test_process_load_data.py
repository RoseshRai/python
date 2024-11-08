import pytest

from helpers.test_helper.setup import db_connection, ena_manager, test_db

EXPECTED_RUNS_DATA = {
    "run_accession": "ERR4796484",
    "instrument_platform": "ILLUMINA",
    "instrument_model": "NextSeq 500",
    "organism_name": "Mycobacterium tuberculosis",
    "last_updated": "2021-10-20 00:00:00",
}

EXPECTED_FILE_LOCATIONS_DATA = [
    {
        "run_id": 1,
        "file_path": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR479/004/ERR4796484/ERR4796484_1.fastq.gz",
        "is_local": True,
    },
    {
        "run_id": 1,
        "file_path": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR479/004/ERR4796484/ERR4796484_2.fastq.gz",
        "is_local": True,
    },
]

@pytest.mark.asyncio
@pytest.mark.parametrize("test_db", ["test_ena_runner_end_to_end"], indirect=True)
async def test_ena_runner_end_to_end(test_db, db_connection, ena_manager):
    async for result in ena_manager.process_and_load_data():
        assert "Processed run ID:" in result, f"Unexpected result: {result}"

    async with db_connection.get_connection() as conn:
        async with conn.cursor() as cur:
            await cur.execute(
                "SELECT run_accession, instrument_platform, instrument_model, organism_name, last_updated FROM ena.runs;",
                (),
            )
            runs_result = await cur.fetchone()
            assert runs_result is not None, "No data found in ena.runs table."
            assert runs_result[0] == EXPECTED_RUNS_DATA["run_accession"]
            assert runs_result[1] == EXPECTED_RUNS_DATA["instrument_platform"]
            assert runs_result[2] == EXPECTED_RUNS_DATA["instrument_model"]
            assert runs_result[3] == EXPECTED_RUNS_DATA["organism_name"]
            assert str(runs_result[4]) == EXPECTED_RUNS_DATA["last_updated"]

    async with db_connection.get_connection() as conn:
        async with conn.cursor() as cur:
            await cur.execute("SELECT run_id, file_path, is_local FROM ena.file_locations ORDER BY id;", ())
            file_locations_result = await cur.fetchall()
            assert len(file_locations_result) == len(
                EXPECTED_FILE_LOCATIONS_DATA
            ), "Unexpected number of entries in ena.file_locations table."

            for row, expected_data in zip(file_locations_result, EXPECTED_FILE_LOCATIONS_DATA):
                assert row[0] == expected_data["run_id"]
                assert row[1] == expected_data["file_path"]
                assert row[2] == expected_data["is_local"]
