import os

from dotenv import load_dotenv
from loguru import logger


class DatabaseConfig:
    def __init__(self):
        load_dotenv()
        self.user: str = os.getenv("DATABASE_USER")
        self.password: str = os.getenv("DATABASE_PASSWORD")
        self.host: str = os.getenv("DATABASE_HOST")
        self.port: str = os.getenv("DATABASE_PORT")
        self.database: str = os.getenv("DATABASE_NAME")
        self._validate_config()

    def _validate_config(self):
        """Ensure all required database configuration values are present."""
        missing = [attr for attr in ["user", "password", "host", "port", "database"] if not getattr(self, attr)]
        if missing:
            logger.error(f"Missing database configuration values: {', '.join(missing)}")
            raise ValueError(f"Missing database configuration values: {', '.join(missing)}")

    def get_postgres_url(self) -> str:
        """Return the PostgresSQL connection string."""
        return f"postgresql://{self.user}:{self.password}@{self.host}:{self.port}/{self.database}"
