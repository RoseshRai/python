#!/bin/sh

poetry run black .
if [ $? -ne 0 ]; then
    echo "Black failed, commit aborted."
    exit 1
fi

poetry run isort .
if [ $? -ne 0 ]; then
    echo "isort failed, commit aborted."
    exit 1
fi

poetry run mypy .
if [ $? -ne 0 ]; then
    echo "mypy failed, commit aborted."
    exit 1
fi