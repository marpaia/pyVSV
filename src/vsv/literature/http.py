"""HTTP transport for literature services.

Live calls sleep between requests. Tests pass their own callables. Nothing here
runs at import.
"""

from __future__ import annotations

import json
import ssl
import time
import urllib.error
import urllib.request
from collections.abc import Callable, Mapping

try:
    import certifi
except ImportError:  # a system trust store is enough when certifi is not installed
    certifi = None  # type: ignore[assignment]

from vsv.errors import InputError

Fetch = Callable[[str, bytes | None], dict[str, object]]
FetchJson = Callable[[str], Mapping[str, object]]
FetchText = Callable[[str], str]

_USER_AGENT = "vsv-corpus"


def live_fetch(url: str, body: bytes | None = None) -> dict[str, object]:
    """GET or POST JSON, with one retry. Sleeps so a full corpus stays under NCBI's rate."""
    payload = _request(url, body)
    time.sleep(0.35)
    parsed = json.loads(payload)
    if not isinstance(parsed, dict):
        raise InputError("The literature service returned a response that is not an object.")
    return parsed


def live_fetch_text(url: str) -> str:
    """GET a document. Sleeps between calls."""
    text = _request(url, None).decode()
    time.sleep(0.35)
    return text


def _request(url: str, body: bytes | None) -> bytes:
    headers = {"User-Agent": _USER_AGENT}
    if body is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    request = urllib.request.Request(url, data=body, headers=headers)
    try:
        return _read(request)
    except urllib.error.URLError:
        pass
    try:
        return _read(request)
    except urllib.error.URLError as error:
        raise InputError(f"Literature request failed: {error}") from error


def _context() -> ssl.SSLContext:
    if certifi is not None:
        return ssl.create_default_context(cafile=certifi.where())
    return ssl.create_default_context()


def _read(request: urllib.request.Request) -> bytes:
    with urllib.request.urlopen(request, timeout=90, context=_context()) as response:
        payload = response.read()
    if not isinstance(payload, bytes):
        raise InputError("The literature service returned a response that is not bytes.")
    return payload
