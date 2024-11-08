import asyncio
from typing import AsyncGenerator, Coroutine

import aiohttp
from loguru import logger


class ENADataFetcher:
    def __init__(self, run_ids: list[str], base_url: str):
        self.run_ids = run_ids
        self.base_url = base_url
        self.semaphore = asyncio.Semaphore(10)

    async def fetch_metadata(self) -> AsyncGenerator[dict[str, str], None]:
        """Fetch metadata for each run ID with limited concurrency."""
        async with aiohttp.ClientSession() as session:

            async def fetch_single(run_id: str) -> dict[str, str] | None:
                """Fetch metadata for a single run ID."""
                async with self.semaphore:
                    params = {
                        "accession": run_id,
                        "result": "read_run",
                        "fields": "run_accession,instrument_platform,instrument_model,scientific_name,last_updated,fastq_ftp",
                        "format": "json",
                    }
                    try:
                        async with session.get(self.base_url, params=params, timeout=10) as response:
                            response.raise_for_status()
                            records = await response.json()
                            if records:
                                logger.debug(f"Fetched data for run ID {run_id}")
                                return records[0]
                    except aiohttp.ClientError as e:
                        logger.error(f"Failed to fetch metadata for run ID {run_id}: {e}")
                    except ValueError as e:
                        logger.error(f"Invalid JSON received for run ID {run_id}: {e}")
                    return None

            tasks: list[Coroutine[any, any, dict[str, str] | None]] = [
                fetch_single(run_id) for run_id in self.run_ids
            ]  # create all coroutine objects. Not yet executed
            for task in asyncio.as_completed(tasks):  # execute each task as a semaphore becomes available
                record = await task
                if record:
                    yield record
