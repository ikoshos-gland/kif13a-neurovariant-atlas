"""Download selected Zeng ... Roussos 2026 files from Synapse (ids from raw/synapse/listing.tsv).

usage: python s03_synapse_download.py small|class [n_parallel]
  small : result tables (coloc_class, specific_eQTL_class, mashr lfsr, specific_and_coloc ...)
  class : the 8 class-level cis-eQTL meta-analysis parquet files (~27 GB)
Each file is md5-checked against the Synapse file handle. Token handling: syn_common only.
"""
from __future__ import annotations

import concurrent.futures as cf
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import syn_common as sc  # noqa: E402

OUT = "/arf/scratch/avural/kif13a_eqtl/raw/synapse"
SMALL = ["coloc_class.tsv", "coloc_and_subclass_specific.tsv", "specific_and_coloc.tsv",
         "specific_eQTL_class.tsv.gz", "mashr_lfsr_class.tsv.gz", "merge_coloc_mashr.tsv"]
CLASSES = ["Oligo", "OPC", "Astro", "EN", "IN", "Immune", "Endo", "Mural"]


def main() -> int:
    what = sys.argv[1]
    npar = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    rows = list(csv.DictReader(open(f"{OUT}/listing.tsv"), delimiter="\t"))
    by_name = {r["name"]: r for r in rows if r["type"] == "FileEntity"}
    if what == "small":
        names = SMALL
    elif what == "class":
        names = [f"meta-{c}-class.parquet" for c in CLASSES]
    else:
        names = what.split(",")
    todo = []
    for n in names:
        r = by_name[n]
        dest = f"{OUT}/{n}"
        if os.path.exists(dest) and os.path.getsize(dest) == int(r["size"] or -1):
            sc.log(f"have {n}")
            continue
        todo.append((r["id"], dest, r["md5"] or None, n))
    results = []
    try:
        with cf.ThreadPoolExecutor(max_workers=npar) as ex:
            futs = {ex.submit(sc.download, i, d, m): n for i, d, m, n in todo}
            for f in cf.as_completed(futs):
                res = f.result()
                res["name"] = futs[f]
                sc.log(json.dumps(res))
                results.append(res)
    except sc.AuthError as e:
        sc.log(str(e))
        open(f"{OUT}/AUTH_FAIL.txt", "w").write(f"status={e.status}\nendpoint={e.endpoint}\n")
        return 3
    with open(f"{OUT}/download_{what if what in ('small', 'class') else 'custom'}.jsonl", "a") as fh:
        for r in results:
            fh.write(json.dumps(r) + "\n")
    bad = [r for r in results if not r["md5_ok"]]
    sc.log(f"done: {len(results)} downloaded, {len(bad)} md5 mismatches")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
