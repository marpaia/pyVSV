.PHONY: setup check lint format format-check typecheck test build

setup:
	uv sync --locked

check: lint format-check typecheck test

lint:
	uv run --locked ruff check .

format:
	uv run --locked ruff format .

format-check:
	uv run --locked ruff format --check .

typecheck:
	uv run --locked mypy

test:
	uv run --locked pytest

build:
	uv build
