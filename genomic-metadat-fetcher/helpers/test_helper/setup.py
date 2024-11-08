from pathlib import Path

import pytest_asyncio

from database.src.async_connection import Database
from helpers.test_helper.async_database_helper import AsyncDatabaseHelper
from managers.ena.src.runner import ENARunner

BASE_DSN = "postgresql://postgres:postgres@localhost"
SCHEMA_DIR: Path = Path(__file__).parent.parent.parent / "database" / "flyway"

TEST_RUN_ID = ["ERR4796484"]
BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"


@pytest_asyncio.fixture
async def test_db(request):
    db_name = request.param
    helper = AsyncDatabaseHelper(BASE_DSN)
    dsn_with_db = await helper.create_test_db(db_name, SCHEMA_DIR)
    yield dsn_with_db
    await helper.drop_test_db(db_name)


@pytest_asyncio.fixture
async def db_connection(test_db):
    db = Database(dsn=test_db)
    await db.connect()
    yield db
    await db.close()


@pytest_asyncio.fixture
async def ena_manager(db_connection):
    return ENARunner(db=db_connection, run_ids=TEST_RUN_ID, base_url=BASE_URL)
