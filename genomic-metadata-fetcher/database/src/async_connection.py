import functools
from contextlib import asynccontextmanager

from loguru import logger
from psycopg.errors import DatabaseError
from psycopg_pool import AsyncConnectionPool


def handle_exceptions(func: callable) -> callable:
    """Decorator to handle and log exceptions in async database methods."""
    @functools.wraps(func) # copies over the original function
    async def wrapper(*args, **kwargs):
        try:
            return await func(*args, **kwargs)
        except DatabaseError as e:
            logger.error(f"Database error in {func.__name__}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in {func.__name__}: {e}")
            return None

    return wrapper


class Database:
    def __init__(self, dsn: str) -> None:
        self.dsn = dsn
        if not self.dsn:
            raise ValueError("Database Data Service Name must be provided as an argument.")

        self.pool: AsyncConnectionPool | None = None
        logger.debug(f"Database initialized with Data Service Name: {self.dsn}")

    @handle_exceptions
    async def connect(self, min_size: int = 1, max_size: int = 5) -> None:
        """Initialize and open the connection pool."""
        self.pool = AsyncConnectionPool(conninfo=self.dsn, min_size=min_size, max_size=max_size, open=False)
        await self.pool.open()
        logger.info(f"Connection pool created successfully with min size {min_size} and max size {max_size}")

    @handle_exceptions
    async def close(self) -> None:
        """Close the connection pool gracefully."""
        if self.pool:
            await self.pool.close()
            logger.info("Connection pool closed successfully")

    @asynccontextmanager # clean way to ensure that the connection is always returned back to the pool. Guarantee cleanup
    async def get_connection(self):
        try:
            if not self.pool:
                logger.debug("Connection pool is not open; calling connect()")
                await self.connect()
            async with self.pool.connection() as conn:
                logger.debug("Acquired connection from pool")
                yield conn
        except DatabaseError as e:
            logger.error(f"Database error in get_connection: {e}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error in get_connection: {e}")
            raise

    @asynccontextmanager
    async def get_cursor(self):
        """Provide an async cursor from a connection in the pool."""
        try:
            async with self.get_connection() as conn:
                async with conn.cursor() as cur:
                    logger.debug("Acquired async cursor from connection")
                    yield cur
                    logger.debug("Async cursor released")
        except DatabaseError as e:
            logger.error(f"Database error in get_cursor: {e}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error in get_cursor: {e}")
            raise

    @handle_exceptions
    async def check_connection(self) -> None:
        """Run a simple query to verify the async connection to the database."""
        logger.debug("Checking async database connection health")
        async with self.get_connection() as cur:
            await cur.execute("SELECT 1")
            logger.info("Async database connection verified successfully")
