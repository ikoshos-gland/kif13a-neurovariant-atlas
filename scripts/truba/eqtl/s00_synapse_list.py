"""List the Zeng ... Roussos 2026 (Nat Genet) Synapse project tree with file sizes and md5.

Writes raw/synapse/listing.tsv and the wiki markdown of the root and of each folder that has one.
Never touches the token beyond syn_common.call(). Run inside a debug job (internet).
"""
from __future__ import annotations

import json
import os
import sys
import urllib.error

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import syn_common as sc  # noqa: E402

ROOT = "syn61929919"
OUT = "/arf/scratch/avural/kif13a_eqtl/raw/synapse"
MAXDEPTH = 4
MAXFILES_PER_FOLDER_SIZES = 400  # fetch sizes for at most this many files per folder


def wiki(syn: str) -> str | None:
    try:
        w = sc.call("GET", f"{sc.REPO}/entity/{syn}/wiki")
        return w.get("markdown")
    except urllib.error.HTTPError:
        return None


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    rows = []
    try:
        ent = sc.call("GET", f"{sc.REPO}/entity/{ROOT}")
        sc.log(f"root {ROOT}: type={ent.get('concreteType')} name={ent.get('name')!r}")
        md = wiki(ROOT)
        if md:
            open(f"{OUT}/wiki_{ROOT}.md", "w").write(md)
        stack = [(ROOT, ent.get("name", ROOT), 0)]
        while stack:
            pid, path, depth = stack.pop()
            kids = sc.children(pid)
            sc.log(f"{path} ({pid}): {len(kids)} children")
            nfile = 0
            for k in sorted(kids, key=lambda x: x["name"]):
                typ = k["type"].split(".")[-1]
                row = {"parent": pid, "path": f"{path}/{k['name']}", "id": k["id"], "type": typ,
                       "name": k["name"], "version": k.get("versionNumber"), "size": "", "md5": "",
                       "modified": k.get("modifiedOn", "")}
                if typ == "FileEntity":
                    nfile += 1
                    if nfile <= MAXFILES_PER_FOLDER_SIZES:
                        try:
                            h = sc.file_handle(k["id"])
                            row["size"] = h.get("contentSize", "")
                            row["md5"] = h.get("contentMd5", "")
                        except urllib.error.HTTPError as e:
                            row["size"] = f"HTTP{e.code}"
                elif typ == "Folder":
                    if depth + 1 < MAXDEPTH:
                        stack.append((k["id"], row["path"], depth + 1))
                    md = wiki(k["id"])
                    if md:
                        open(f"{OUT}/wiki_{k['id']}.md", "w").write(md)
                rows.append(row)
    except sc.AuthError as e:
        sc.log(str(e))
        open(f"{OUT}/AUTH_FAIL.txt", "w").write(f"status={e.status}\nendpoint={e.endpoint}\n")
        return 3
    cols = ["parent", "path", "id", "type", "name", "version", "size", "md5", "modified"]
    with open(f"{OUT}/listing.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    sc.log(f"wrote {len(rows)} rows to {OUT}/listing.tsv")
    json.dump({"root": ROOT, "n": len(rows)}, open(f"{OUT}/listing.done.json", "w"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
