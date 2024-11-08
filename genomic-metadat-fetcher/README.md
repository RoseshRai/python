# Genomic Metadata fetcher

This project provides a data pipeline to download metadata and associated FastQ files from the European Nucleotide Archive (ENA) for specified sample runs. It saves metadata in a PostgreSQL database, with configuration settings and custom modules that organize the pipeline’s functionality.

## Project Structure

- **configuration**: Contains utilities for reading environment variables, particularly for database connection settings.
- **database**: Manages the database connection setup and utilities.
- **helpers**: Includes setup utilities, such as a test helper for initializing and managing a temporary test database.
- **main**: The entry point for the pipeline. The `main.py` file initializes the database, runs the pipeline, and processes data.
- **managers**: Contains core ENA data management classes, including classes to fetch and store ENA data.
- **tests**: Includes an integration/end-to-end test to verify pipeline functionality. Currently, only partial test coverage is implemented due to time constraints.


## Database Configuration

This project uses PostgreSQL for its database.

### Database Setup

Currently, the database schema is not fully managed by Flyway migrations. To run this project locally:
1. **Create the Database Manually**: Create a PostgreSQL database that matches the connection details specified in the `.env` file.
2. **Update the `.env` File**: Configure your local database connection in the `.env` file with values for `DATABASE_USER`, `DATABASE_PASSWORD`, `DATABASE_HOST`, `DATABASE_PORT`, and `DATABASE_NAME`.

## Running the Pipeline

To execute the data pipeline:
1. Ensure that the database is set up, and the `.env` file has the correct connection details.
2. Run the `main.py` file in the `main/src` directory.

### Running the Project in PyCharm

To run the project in PyCharm:
1. Simply press the **Play** button in the IDE.

### Running the Project in Terminal

If running from the terminal, ensure your Python path is set up correctly:
```bash
export PYTHONPATH="/home/{your_user}/{project_path}/genomic-metadata-fetcher:$PYTHONPATH"
python main/src/main.py
