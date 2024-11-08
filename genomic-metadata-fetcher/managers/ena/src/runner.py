import asyncio
from datetime import datetime
from typing import AsyncGenerator

from loguru import logger

from database.src.async_connection import Database

from .fetch_data import ENADataFetcher


class ENARunner:
    def __init__(self, db: Database, run_ids: list[str], base_url: str):
        self.db = db
        self.ena_fetcher = ENADataFetcher(run_ids, base_url)

    async def process_and_load_data(self) -> AsyncGenerator[str, None]:
        """Fetch metadata and load it into the database concurrently, yielding results as they complete."""

        async def process_record(record_data) -> str:
            await self._store_run_data(record_data)
            return f"Processed run ID: {record_data['run_accession']}"

        tasks = []

        async for record in self.ena_fetcher.fetch_metadata():
            task = asyncio.create_task(process_record(record))
            tasks.append(task)

        for task in asyncio.as_completed(tasks):
            result = await task
            logger.info(f"Task completed: {result}")
            yield result

    async def _store_run_data(self, record: dict[str, str]) -> None:
        """Store run metadata and file paths associated with the run asynchronously."""
        metadata_query = """
            INSERT INTO ena.runs (run_accession, instrument_platform, instrument_model, organism_name, last_updated)
            VALUES (%s, %s, %s, %s, %s)
            ON CONFLICT (run_accession) 
            DO UPDATE SET 
                instrument_platform = EXCLUDED.instrument_platform,
                instrument_model = EXCLUDED.instrument_model,
                organism_name = EXCLUDED.organism_name,
                last_updated = EXCLUDED.last_updated;
        """
        metadata_params: tuple[str, str, str, str, datetime] = (
            record["run_accession"],
            record["instrument_platform"],
            record["instrument_model"],
            record["scientific_name"],
            datetime.strptime(record["last_updated"], "%Y-%m-%d"),
        )

        run_accession: str = record["run_accession"]
        file_paths: list[str] = record.get("fastq_ftp", "").split(";") if record.get("fastq_ftp") else []
        file_path_queries: list[tuple[str, tuple[str, bool, str]]] = [
            (
                """
                INSERT INTO ena.file_locations (run_id, file_path, is_local)
                SELECT id, %s, %s FROM ena.runs WHERE run_accession = %s
                ON CONFLICT (run_id, file_path) 
                DO UPDATE SET 
                    is_local = EXCLUDED.is_local
                WHERE ena.file_locations.is_local != EXCLUDED.is_local;
                """,
                (path, not path.startswith("ftp://"), run_accession),
            )
            for path in file_paths
        ]

        try:
            async with self.db.get_cursor() as cur:
                await cur.execute(metadata_query, metadata_params)
                logger.debug(f"Inserted metadata for run ID: {run_accession}")

                for query, params in file_path_queries:
                    await cur.execute(query, params)
                    logger.debug(f"Inserted file path {params[0]} for run ID: {run_accession}")

                logger.info(f"Successfully stored data for run ID: {run_accession}")
        except Exception as e:
            logger.error(f"Failed to store data for run ID {run_accession}: {e}")
            raise
