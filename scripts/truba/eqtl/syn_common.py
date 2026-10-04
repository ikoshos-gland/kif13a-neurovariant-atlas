"""Synapse REST helpers shared by the listing and download scripts.

The Personal Access Token is read programmatically from /arf/home/avural/.synapseConfig and is
never printed, logged or written anywhere. On HTTP 401 or 403 the caller stops that branch and
records only the status code and the endpoint (no retry loop), as agreed with the user.
"""
from __future__ import annotations

import configparser
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

REPO = "https://repo-prod.prod.sagebase.org/repo/v1"
FILE = "https://repo-prod.prod.sagebase.org/file/v1"
CFG = "/arf/home/avural/.synapseConfig"


class AuthError(RuntimeError):
    def __init__(self, status: int, endpoint: str):
        super().__init__(f"AUTH_FAIL status={status} endpoint={endpoint}")
        self.status = status
        self.endpoint = endpoint


def _token() -> str:
    cp = configparser.ConfigParser()
    with open(CFG, encoding="utf-8") as fh:
        cp.read_file(fh)
    for sec in ("authentication", "Authentication"):
        if cp.has_section(sec):
            for key in ("authtoken", "authToken", "auth_token"):
                if cp.has_option(sec, key):
                    tok = cp.get(sec, key).strip()
                    if tok:
                        return tok
    raise RuntimeError("authtoken not found in the [authentication] section of the config file")


_TOK: str | None = None


def call(method: str, url: str, body: dict | None = None, raw: bool = False, tries: int = 4):
    """One authenticated REST call. Retries only transient 5xx/timeouts, never 401/403."""
    global _TOK
    if _TOK is None:
        _TOK = _token()
    data = json.dumps(body).encode() if body is not None else None
    endpoint = url.split("?")[0].replace(REPO, "repo/v1").replace(FILE, "file/v1")
    for attempt in range(tries):
        req = urllib.request.Request(url, data=data, method=method)
        req.add_header("Authorization", "Bearer " + _TOK)
        req.add_header("Accept", "application/json")
        if data is not None:
            req.add_header("Content-Type", "application/json")
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                payload = r.read()
                return payload.decode() if raw else json.loads(payload)
        except urllib.error.HTTPError as e:
            if e.code in (401, 403):
                raise AuthError(e.code, method + " " + endpoint) from None
            if e.code in (404,):
                raise
            if attempt == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries - 1:
                raise
        time.sleep(5 * (attempt + 1))


def children(parent: str) -> list[dict]:
    out, tok = [], None
    while True:
        body = {"parentId": parent, "includeTypes": ["folder", "file", "link", "table"]}
        if tok:
            body["nextPageToken"] = tok
        res = call("POST", f"{REPO}/entity/children", body)
        out.extend(res.get("page", []))
        tok = res.get("nextPageToken")
        if not tok:
            return out


def file_handle(syn: str) -> dict:
    res = call("GET", f"{REPO}/entity/{syn}/filehandles")
    lst = res.get("list", [])
    # the entity's own data file handle is the non-preview one
    lst = [h for h in lst if not h.get("isPreview") and "Preview" not in h.get("concreteType", "")]
    return lst[0] if lst else {}


def presigned(syn: str) -> str:
    return call("GET", f"{REPO}/entity/{syn}/file?redirect=false", raw=True).strip().strip('"')


def download(syn: str, dest: str, expect_md5: str | None = None, resume_tries: int = 6) -> dict:
    """Download with HTTP Range resume of an existing .part file (presigned URLs are short-lived,
    so a fresh one is requested for every attempt). md5 is computed over the finished file."""
    tmp = dest + ".part"
    t0 = time.time()
    for attempt in range(resume_tries):
        have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        url = presigned(syn)
        req = urllib.request.Request(url)
        if have:
            req.add_header("Range", f"bytes={have}-")
        try:
            # presigned S3 URL: no Authorization header (the signature carries the grant)
            with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "ab" if have else "wb") as fh:
                if have and r.status != 206:
                    fh.seek(0)
                    fh.truncate()
                while True:
                    buf = r.read(1 << 22)
                    if not buf:
                        break
                    fh.write(buf)
            break
        except (OSError, urllib.error.URLError) as e:  # ssl eof, broken pipe, reset
            log(f"{syn}: transfer interrupted after {os.path.getsize(tmp)} bytes ({type(e).__name__}), resuming")
            time.sleep(10)
    md5 = hashlib.md5()
    n = 0
    with open(tmp, "rb") as fh:
        while True:
            buf = fh.read(1 << 24)
            if not buf:
                break
            md5.update(buf)
            n += len(buf)
    got = md5.hexdigest()
    ok = expect_md5 is None or got == expect_md5
    if ok:
        os.replace(tmp, dest)
    return {"syn": syn, "dest": dest, "bytes": n, "md5": got, "md5_expected": expect_md5,
            "md5_ok": ok, "seconds": round(time.time() - t0, 1)}


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)
    sys.stdout.flush()
