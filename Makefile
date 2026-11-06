.PHONY: setup check lint format format-check typecheck test build check-dist

DIST_DIR ?= dist

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
	uv build --no-sources --out-dir "$(DIST_DIR)"

check-dist: build
	uv run --locked twine check --strict "$(DIST_DIR)"/*.whl "$(DIST_DIR)"/*.tar.gz
