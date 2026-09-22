"""HTTP API. On Vercel the app is mounted at /api and also answers / when run alone."""

from __future__ import annotations

from fastapi import FastAPI

app = FastAPI(
    title="VSV",
    docs_url="/api/docs",
    redoc_url="/api/redoc",
    openapi_url="/api/openapi.json",
)


@app.get("/")
def read_root() -> dict[str, str]:
    """Identify the process when this app is serving its own origin."""
    return {"service": "vsv-api", "health": "/api/health"}


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "vsv-api"}
