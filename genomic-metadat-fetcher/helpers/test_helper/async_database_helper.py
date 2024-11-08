from pathlib import Path

import psycopg
from loguru import logger
from psycopg.sql import SQL, Identifier


class AsyncDatabaseHelper:
    def __init__(self, base_dsn: str):
        self.base_dsn = base_dsn

    async def create_test_db(self, db_name: str, schema_dir: str | Path) -> str:
        await self.drop_test_db(db_name)
        conn = await psycopg.AsyncConnection.connect(self.base_dsn)
        await conn.set_autocommit(True)
        try:
            async with conn.cursor() as cur:
                await cur.execute(SQL("CREATE DATABASE {}").format(Identifier(db_name)))
        finally:
            await conn.close()

        dsn_with_db = f"{self.base_dsn}/{db_name}"
        logger.info(schema_dir)
        await self._apply_schema(dsn_with_db, schema_dir)
        return dsn_with_db

    async def drop_test_db(self, db_name: str) -> None:
        conn = await psycopg.AsyncConnection.connect(self.base_dsn)
        await conn.set_autocommit(True)
        try:
            async with conn.cursor() as cur:
                await cur.execute(SQL("DROP DATABASE IF EXISTS {}").format(Identifier(db_name)))
        finally:
            await conn.close()

    @staticmethod
    async def _apply_schema(dsn: str, schema_dir: str | Path) -> None:
        schema_files = sorted(Path(schema_dir).glob("*.sql"), key=lambda p: p.stem)
        conn = await psycopg.AsyncConnection.connect(dsn)
        await conn.set_autocommit(True)
        try:
            async with conn.cursor() as cur:
                for schema_file in schema_files:
                    with open(schema_file, "r") as file:
                        sql = file.read()
                        await cur.execute(SQL(sql))
                    logger.info(f"Successfully executed: {schema_file.name}")
        except Exception as e:
            logger.error(f"Error executing schema file {schema_file.name}: {e}")
        finally:
            await conn.close()
