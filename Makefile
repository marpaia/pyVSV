.PHONY: setup check lint format format-check typecheck test build api frontend dev

setup:
	uv sync --locked
	pnpm --dir frontend install --frozen-lockfile

check: lint format-check typecheck test

lint:
	uv run --locked ruff check .

format:
	uv run --locked ruff format .

format-check:
	uv run --locked ruff format --check .

typecheck:
	uv run --locked --group server mypy

test:
	uv run --locked --group server pytest

API_DEV = uv run --locked --group server uvicorn vsv.server:app --reload --port 8000
WEB_DEV = pnpm --dir frontend run dev

api:
	$(API_DEV)

frontend:
	$(WEB_DEV)

# Each server is its own process group. Ctrl-C kills both groups.
dev:
	@echo "API http://127.0.0.1:8000    App http://localhost:3000"
	@set -m; \
	cleanup() { \
		trap - INT TERM EXIT; \
		[ -n "$${api:-}" ] && kill -TERM -"$$api" 2>/dev/null || true; \
		[ -n "$${web:-}" ] && kill -TERM -"$$web" 2>/dev/null || true; \
		exit 130; \
	}; \
	trap cleanup INT TERM EXIT; \
	$(API_DEV) & api=$$!; \
	$(WEB_DEV) & web=$$!; \
	wait

build:
	uv build
