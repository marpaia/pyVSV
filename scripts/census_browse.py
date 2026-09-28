"""Download publisher copies through the signed-in debugging Chrome.

Connects to Chrome at http://127.0.0.1:9222. Does not close that browser.
Writes data/census/acquired/{pmid}.html and .txt, sections the text, and
updates manifest.jsonl after each paper so a stopped run can resume.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import re
import sys
import time
from dataclasses import replace
from datetime import date
from pathlib import Path
from urllib.parse import urlparse

from playwright.sync_api import Page, sync_playwright
from playwright.sync_api import TimeoutError as PlaywrightTimeout

from vsv.ingest import redact
from vsv.literature.catalog import ManifestRow
from vsv.literature.index import CENSUS_ROOT, read_manifest, write_manifest
from vsv.literature.section import section_text, write_sections

CDP = "http://127.0.0.1:9222"
QUEUE = frozenset({"no-pmc", "front-matter-only", "error", "needs-login"})
DONE = frozenset({"xml-body", "browser-html", "browser-pdf", "unavailable"})
_NCBI = ("ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov")
_PII = re.compile(r"pii/([A-Za-z0-9]+)")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    root = CENSUS_ROOT
    acquired = root / "acquired"
    sections = root / "sections"
    acquired.mkdir(parents=True, exist_ok=True)
    rows = list(read_manifest(root))
    wanted = set(args.only)
    queue = [row for row in rows if _pending(row, acquired) and (not wanted or row.pmid in wanted)]
    if args.limit:
        queue = queue[: args.limit]
    print(f"queue {len(queue)}", flush=True)
    if not queue:
        return 0
    by_pmid = {row.pmid: index for index, row in enumerate(rows)}
    login_run = 0
    with sync_playwright() as playwright:
        browser = playwright.chromium.connect_over_cdp(CDP)
        page = browser.contexts[0].new_page()
        page.set_default_timeout(30_000)
        page.set_default_navigation_timeout(60_000)
        try:
            for row in queue:
                status, detail, source, text, html = _one(page, row.pmid, row.doi, acquired)
                digest = ""
                nbytes = 0
                html_path = acquired / f"{row.pmid}.html"
                if html:
                    html_path.write_text(html, encoding="utf-8")
                    raw = html_path.read_bytes()
                    digest = hashlib.sha256(raw).hexdigest()
                    nbytes = len(raw)
                if text:
                    (acquired / f"{row.pmid}.txt").write_text(redact(text), encoding="utf-8")
                    blocks = section_text(text)
                    write_sections(sections / f"{row.pmid}.json", row.pmid, blocks)
                updated = replace(
                    row,
                    source=source,
                    retrieved=date.today().isoformat(),
                    sha256=digest,
                    nbytes=nbytes,
                    body_characters=len(" ".join(text.split())),
                    status=status,
                    detail=detail,
                )
                rows[by_pmid[row.pmid]] = updated
                write_manifest(root, rows)
                print(f"{row.pmid} {status} {len(text)} {source[:120]}", flush=True)
                if status == "needs-login":
                    login_run += 1
                    if login_run >= 3:
                        print("stopped: three logins in a row", flush=True)
                        return 2
                else:
                    login_run = 0
                time.sleep(2)
        finally:
            page.close()
    return 0


def _one(page: Page, pmid: str, doi: str, acquired: Path) -> tuple[str, str, str, str, str]:
    try:
        page.goto(f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", wait_until="domcontentloaded")
        page.wait_for_timeout(1200)
        link = _publisher_link(page)
        if not link and doi:
            link = f"https://doi.org/{doi}"
        if not link:
            page.screenshot(path=str(acquired / f"{pmid}.png"))
            return "unavailable", "no publisher link", page.url, "", ""
        page.goto(link, wait_until="domcontentloaded")
        page.wait_for_timeout(2500)
        _accept_cookies(page)
        _deepen(page, link)
        _accept_cookies(page)
        response = page.url
        content_type = ""
        text = _article_text(page)
        html = page.content()
        status = _status(text)
        if status != "browser-html":
            page.screenshot(path=str(acquired / f"{pmid}.png"))
        return status, content_type, response, text, html
    except PlaywrightTimeout:
        with contextlib.suppress(Exception):
            page.screenshot(path=str(acquired / f"{pmid}.png"))
        return "error", "timeout", page.url, "", ""
    except Exception as error:  # a publisher page should not stop the queue
        return "error", str(error)[:200], "", "", ""


def _publisher_link(page: Page) -> str:
    found = page.eval_on_selector_all(
        ".full-text-links-list a",
        "els => els.map(a => a.href).filter(Boolean)",
    )
    for href in found:
        host = urlparse(href).netloc.removeprefix("www.")
        if host and not any(host.endswith(name) for name in _NCBI):
            return str(href)
    return ""


def _pending(row: ManifestRow, acquired: Path) -> bool:
    if row.status in {"xml-body", "browser-html", "browser-pdf"}:
        return False
    if row.status == "unavailable":
        text = acquired / f"{row.pmid}.txt"
        return not text.is_file() or text.stat().st_size < 4000
    return row.status in QUEUE


def _deepen(page: Page, origin: str) -> None:
    """Prefer the ScienceDirect HTML article over an abstract or linking hub."""
    url = page.url
    if "/abs/pii/" in url:
        page.goto(url.replace("/abs/pii/", "/pii/"), wait_until="domcontentloaded")
        page.wait_for_timeout(2000)
        return
    match = _PII.search(origin) or _PII.search(url)
    thin = len(_article_text(page)) < 4000
    if match and (thin or "cell.com" in url or "linkinghub.elsevier.com" in url):
        pii = match.group(1)
        page.goto(
            f"https://www.sciencedirect.com/science/article/pii/{pii}",
            wait_until="domcontentloaded",
        )
        page.wait_for_timeout(2500)


def _accept_cookies(page: Page) -> None:
    for label in ("Accept all cookies", "Accept cookies"):
        button = page.get_by_role("button", name=label)
        if button.count():
            try:
                button.first.click(timeout=1500)
            except Exception:
                return


def _article_text(page: Page) -> str:
    best = ""
    for selector in ("article", "#body", "main"):
        locator = page.locator(selector)
        if not locator.count():
            continue
        try:
            text = locator.first.inner_text(timeout=8000)
        except Exception:
            continue
        if len(text) > len(best):
            best = text
    if len(best) < 1500:
        with contextlib.suppress(Exception):
            best = page.locator("body").inner_text(timeout=8000)
    return best


def _status(text: str) -> str:
    folded = text.lower()
    head = folded[:2500]
    if len(text) >= 8000 and "view pdf" in folded:
        return "browser-html"
    if len(text) >= 12000:
        return "browser-html"
    wall = "access through your organization" in head or "purchase pdf" in head
    if wall:
        return "needs-login"
    if any(phrase in head for phrase in ("interlibrary loan", "request this article")):
        return "unavailable"
    if len(" ".join(text.split())) >= 4000:
        return "browser-html"
    return "unavailable"


if __name__ == "__main__":
    sys.exit(main())
