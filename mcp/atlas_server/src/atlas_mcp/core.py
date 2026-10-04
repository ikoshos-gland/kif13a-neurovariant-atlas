"""Shared plumbing: project paths, secrets, gene targets, HTTP with retries,
and provenance-tracked storage under data/raw.

Every file a tool writes gets a sidecar ``<file>.prov.json`` and a line in
``data/raw/PROVENANCE.jsonl`` (source URLs, parameters, UTC time, sha256), so a
figure can always be traced back to the exact upstream pull.
"""

from __future__ import annotations

import asyncio
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import logging
import os
import time
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlparse

import httpx
import yaml

PACKAGE_VERSION = "0.1.0"
USER_AGENT = f"kif13a-neurovariant-atlas/{PACKAGE_VERSION} (research; github.com/ikoshos-gland)"
RETRY_STATUS = {429, 500, 502, 503, 504}

log = logging.getLogger("atlas_mcp")


# --------------------------------------------------------------------------- paths


def project_root() -> Path:
    env = os.environ.get("ATLAS_PROJECT_ROOT")
    if env:
        return Path(env).resolve()
    # src/atlas_mcp/core.py -> src -> atlas_server -> mcp -> project root
    return Path(__file__).resolve().parents[4]


ROOT = project_root()
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def today() -> str:
    return dt.date.today().strftime("%Y%m%d")


# --------------------------------------------------------------------------- secrets


def load_dotenv() -> None:
    """Read KEY=VALUE lines from <project>/.env without overriding real env vars."""
    env_file = ROOT / ".env"
    if not env_file.exists():
        return
    for line in env_file.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key, value = key.strip(), value.strip().strip('"').strip("'")
        if value and key not in os.environ:
            os.environ[key] = value


def secret(name: str) -> str | None:
    value = os.environ.get(name, "").strip()
    return value or None


# --------------------------------------------------------------------------- targets


def load_targets() -> dict[str, Any]:
    with open(ROOT / "config" / "targets.yaml", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def target(gene: str) -> dict[str, Any]:
    """Return the verified identifiers for one of the project genes."""
    genes = load_targets()["genes"]
    for symbol, info in genes.items():
        if symbol.upper() == gene.strip().upper():
            return {"symbol": symbol, **info}
    raise ValueError(
        f"{gene!r} is not in config/targets.yaml (known: {', '.join(genes)}). "
        "Add it there with verified MANE identifiers first."
    )


def unversioned(identifier: str) -> str:
    return identifier.split(".")[0]


# --------------------------------------------------------------------------- HTTP

_host_locks: dict[str, asyncio.Lock] = {}
_host_last: dict[str, float] = {}


def _min_interval(host: str) -> float:
    if host.endswith("ncbi.nlm.nih.gov"):
        return 0.11 if secret("NCBI_API_KEY") else 0.35  # 10 vs 3 requests/s
    if host.endswith("ensembl.org"):
        return 0.07  # 15 requests/s
    return 0.0


async def request(
    method: str,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    json_body: Any = None,
    headers: dict[str, str] | None = None,
    timeout: float = 120.0,
    retries: int = 5,
) -> httpx.Response:
    """HTTP with per-host politeness spacing and exponential backoff."""
    host = urlparse(url).hostname or ""
    lock = _host_locks.setdefault(host, asyncio.Lock())
    hdrs = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    hdrs.update(headers or {})
    last_error: Exception | None = None

    for attempt in range(retries):
        async with lock:
            wait = _host_last.get(host, 0.0) + _min_interval(host) - time.monotonic()
            if wait > 0:
                await asyncio.sleep(wait)
            _host_last[host] = time.monotonic()
        try:
            async with httpx.AsyncClient(timeout=timeout, follow_redirects=True) as client:
                resp = await client.request(method, url, params=params, json=json_body, headers=hdrs)
        except (httpx.TransportError, httpx.TimeoutException) as exc:
            last_error = exc
            await asyncio.sleep(2**attempt)
            continue
        if resp.status_code in RETRY_STATUS and attempt < retries - 1:
            retry_after = resp.headers.get("Retry-After")
            delay = float(retry_after) if retry_after and retry_after.replace(".", "").isdigit() else 2**attempt
            log.warning("%s %s -> %s, retrying in %.1fs", method, url, resp.status_code, delay)
            await asyncio.sleep(min(delay, 60))
            continue
        resp.raise_for_status()
        return resp

    raise RuntimeError(f"{method} {url} failed after {retries} attempts: {last_error}")


async def get_json(url: str, **kw: Any) -> Any:
    return (await request("GET", url, **kw)).json()


async def post_json(url: str, body: Any, **kw: Any) -> Any:
    return (await request("POST", url, json_body=body, **kw)).json()


async def get_bytes(url: str, **kw: Any) -> bytes:
    kw.setdefault("headers", {"Accept": "*/*"})
    return (await request("GET", url, **kw)).content


# --------------------------------------------------------------------------- storage


def _write_provenance(path: Path, content: bytes, prov: dict[str, Any]) -> dict[str, Any]:
    record = {
        "file": rel(path),
        "sha256": hashlib.sha256(content).hexdigest(),
        "bytes": len(content),
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "server": f"kif13a-atlas-mcp {PACKAGE_VERSION}",
        **prov,
    }
    path.with_name(path.name + ".prov.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    RAW.mkdir(parents=True, exist_ok=True)
    with open(RAW / "PROVENANCE.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(record) + "\n")
    return record


def save_bytes(subdir: str, filename: str, content: bytes, *, base: Path = RAW, **prov: Any) -> dict[str, Any]:
    folder = base / subdir
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / filename
    path.write_bytes(content)
    record = _write_provenance(path, content, prov)
    return {"path": record["file"], "sha256": record["sha256"][:16], "bytes": record["bytes"]}


def save_json(subdir: str, filename: str, obj: Any, *, compress: bool = False, **prov: Any) -> dict[str, Any]:
    content = json.dumps(obj, indent=None if compress else 1, ensure_ascii=False).encode("utf-8")
    if compress:
        content = gzip.compress(content)
        filename = filename if filename.endswith(".gz") else filename + ".gz"
    return save_bytes(subdir, filename, content, **prov)


def save_tsv(
    subdir: str, filename: str, rows: list[dict[str, Any]], *, base: Path = RAW, **prov: Any
) -> dict[str, Any]:
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: _cell(v) for k, v in row.items()})
    prov.setdefault("n_records", len(rows))
    return save_bytes(subdir, filename, buf.getvalue().encode("utf-8"), base=base, **prov)


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(v) for v in value)
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False)
    return value


def inventory(limit: int = 200) -> list[dict[str, Any]]:
    ledger = RAW / "PROVENANCE.jsonl"
    if not ledger.exists():
        return []
    records = [json.loads(line) for line in ledger.read_text(encoding="utf-8").splitlines() if line.strip()]
    latest: dict[str, dict[str, Any]] = {}
    for rec in records:
        latest[rec["file"]] = rec
    keep = ("file", "retrieved_utc", "bytes", "n_records", "source", "tool")
    return [{k: r.get(k) for k in keep if k in r} for r in list(latest.values())[-limit:]]


# --------------------------------------------------------------------------- misc

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Pyl": "O",
}


def parse_positions(spec: str | None) -> set[int] | None:
    """'300-360,410' -> {300..360, 410}. None/empty -> None (no filter)."""
    if not spec:
        return None
    out: set[int] = set()
    for part in spec.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(part))
    return out


def chunks(items: list[Any], size: int) -> Iterable[list[Any]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]
