"""Health routes for the HTTP API."""

from fastapi.testclient import TestClient

from vsv.server import app

client = TestClient(app)


def test_health() -> None:
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.json() == {"status": "ok", "service": "vsv-api"}


def test_root_points_at_health() -> None:
    response = client.get("/")
    assert response.status_code == 200
    assert response.json() == {"service": "vsv-api", "health": "/api/health"}
