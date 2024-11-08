import asyncio

from loguru import logger

from configuration.src.configuration import DatabaseConfig
from database.src.async_connection import Database
from managers.ena.src.runner import ENARunner

ENA_RUN_IDS = ["ERR4796484", "ERR5396044", "SRR2097047", "SRX5798521", "SRR9213986"]
BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"


async def run_pipeline():
    """Main function to run the data fetching and storing pipeline."""
    configuration = DatabaseConfig()
    postgres_url: str = configuration.get_postgres_url()

    db = Database(dsn=postgres_url)
    await db.check_connection()

    ena_runner = ENARunner(db, ENA_RUN_IDS, BASE_URL)

    try:
        # Consume the async generator to process all records
        async for result in ena_runner.process_and_load_data():
            logger.info(result)
    except Exception as e:
        logger.error(f"An error occurred during the pipeline execution: {e}")
    finally:
        await db.close()


if __name__ == "__main__":
    logger.info("Starting ENA Data Fetching Pipeline...")
    asyncio.run(run_pipeline())
