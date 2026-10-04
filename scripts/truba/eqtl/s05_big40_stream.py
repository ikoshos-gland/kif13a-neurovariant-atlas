"""Stream a BIG40 33k sumstats file and keep only the chr6 analysis window (debug job, internet).

BIG40 stats33k files are single-member gzip without BGZF sync points (docs/ak_madde_goruntuleme.md
trap 5), so the prefix up to chr6 ~19 Mb has to be streamed; the stream is closed as soon as the
window has been passed. Output has the same layout as s10_extract.extract_big40().
usage: python s05_big40_stream.py <pheno>
"""
from __future__ import annotations

import gzip
import sys
import time
import urllib.request
import zlib

W = "/arf/scratch/avural/kif13a_eqtl"
K_LO37, K_HI37 = 17_759_414, 17_987_854
WIN_LO, WIN_HI = K_LO37 - 1_150_000, K_HI37 + 1_150_000
BIG40 = "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k"


def main(p: str) -> int:
    url = f"{BIG40}/{p}.txt.gz"
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    tail = b""
    header, kept = None, []
    got, t0, done = 0, time.time(), False
    req = urllib.request.Request(url, headers={"User-Agent": "kif13a-atlas/1.0 (academic; TRUBA)"})
    with urllib.request.urlopen(req, timeout=600) as r:
        total = int(r.headers.get("Content-Length", 0))
        while not done:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            got += len(chunk)
            data = tail + dec.decompress(chunk)
            *lines, tail = data.split(b"\n")
            for raw in lines:
                f = raw.decode().split()
                if header is None:
                    header = f
                    continue
                if f[0] != "06":
                    if kept or (f[0].isdigit() and int(f[0]) > 6):
                        done = True
                        break
                    continue
                pos = int(f[2])
                if pos < WIN_LO:
                    continue
                if pos > WIN_HI:
                    done = True
                    break
                kept.append(f)
            if got % (50 << 20) < (1 << 20):
                print(f"{p}: {got/1048576:.0f}/{total/1048576:.0f} MB, {len(kept)} rows, {time.time()-t0:.0f}s", flush=True)
    with gzip.open(f"{W}/work/big40_{p}.tsv.gz", "wt") as fh:
        fh.write("\t".join(header) + "\n")
        for f in kept:
            fh.write("\t".join(f) + "\n")
    print(f"{p}: wrote {len(kept)} rows, streamed {got/1048576:.0f} of {total/1048576:.0f} MB in {time.time()-t0:.0f}s",
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
