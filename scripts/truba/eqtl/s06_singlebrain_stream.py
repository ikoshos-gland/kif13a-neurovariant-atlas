"""SingleBrain (Zenodo 14908182) narrow extraction with HTTP Range requests (debug job, internet).

The first attempt (streaming whole 4-5 GB files) failed: the files are multi-member gzip (BGZF
style) and the single zlib decoder stopped after the first member while unused_data grew until the
job was OOM-killed. This version never streams a whole file:
  1. probe: sample K evenly spaced byte offsets, sync to the next gzip member header, inflate a
     few members and read the first complete line -> (chr, pos, feature) per offset;
  2. if the samples are sorted by (chr, pos), bisect to the byte range that covers each target
     region and inflate only that range, keeping rows whose feature is a target gene.
usage: python s06_singlebrain_stream.py <celltype prefix, e.g. OD>
"""
from __future__ import annotations

import gzip
import os
import sys
import time
import urllib.request
import zlib

W = "/arf/scratch/avural/kif13a_eqtl"
ZEN = "https://zenodo.org/records/14908182/files"
UA = {"User-Agent": "kif13a-atlas/1.0 (academic; TRUBA)"}
# hg38 regions (gene +-1 Mb plus margin). KIF13A GRCh38 chr6:17,758,198-17,988,274 (MANE).
REGIONS = {"6": (16_600_000, 19_150_000)}
TARGETS = {"KIF13A", "NUP153", "NHLRC1", "TPMT", "ENSG00000137177", "ENSG00000124789",
           "ENSG00000187566", "ENSG00000137364"}
GZ_MAGIC = b"\x1f\x8b\x08"


def http_range(url: str, a: int, b: int, tries: int = 5) -> bytes:
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={**UA, "Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except OSError:
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"range {a}-{b} failed")


def size_of(url: str) -> int:
    req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        cr = r.headers.get("Content-Range", "")
        return int(cr.split("/")[-1])


def inflate_from(buf: bytes, start: int) -> tuple[bytes, int]:
    """Inflate consecutive gzip members from buf[start:]; return text and bytes consumed."""
    out, i = [], start
    while i < len(buf):
        if buf[i:i + 3] != GZ_MAGIC:
            break
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            out.append(d.decompress(buf[i:]))
        except zlib.error:
            break
        if not d.eof:
            break  # truncated member at the end of the buffer
        i = len(buf) - len(d.unused_data)
    return b"".join(out), i


def first_member(buf: bytes) -> int:
    j = 0
    while True:
        j = buf.find(GZ_MAGIC, j)
        if j < 0:
            return -1
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            d.decompress(buf[j:j + 70000])
            if d.eof or len(buf) - j < 70000:
                return j
        except zlib.error:
            pass
        j += 1


class Reader:
    def __init__(self, url: str, header: list[str]):
        self.url, self.h = url, header
        low = [x.lower() for x in header]
        self.fi = low.index("feature")
        self.ci = low.index("chr")
        self.pi = low.index("pos")

    def key_at(self, off: int):
        buf = http_range(self.url, off, off + 400_000)
        j = first_member(buf)
        if j < 0:
            return None
        txt, _ = inflate_from(buf, j)
        lines = txt.split(b"\n")
        for ln in lines[1:-1]:
            f = ln.decode().split("\t")
            if len(f) > self.pi:
                c = f[self.ci].replace("chr", "")
                try:
                    return (int(c) if c.isdigit() else 99, int(float(f[self.pi])), f[self.fi], off + j)
                except ValueError:
                    continue
        return None


def main(ct: str) -> int:
    url = f"{ZEN}/{ct}_eqtl_full_assoc.tsv.gz?download=1"
    out = f"{W}/work/singlebrain_{ct}.tsv.gz"
    size = size_of(url)
    head = http_range(url, 0, 300_000)
    txt, _ = inflate_from(head, 0)
    header = txt.split(b"\n", 1)[0].decode().split("\t")
    rd = Reader(url, header)
    print(f"{ct}: size {size/1e9:.2f} GB, header {header[:8]}", flush=True)
    K = 48
    samples = []
    for k in range(K):
        s = rd.key_at(int(size * k / K))
        if s:
            samples.append(s)
    keys = [(c, p) for c, p, _, _ in samples]
    sorted_pos = all(keys[i] <= keys[i + 1] for i in range(len(keys) - 1))
    feats = [f for _, _, f, _ in samples]
    sorted_feat = all(feats[i] <= feats[i + 1] for i in range(len(feats) - 1))
    print(f"{ct}: {len(samples)} samples, sorted by (chr,pos)={sorted_pos}, by feature={sorted_feat}", flush=True)
    for s in samples[::6]:
        print(f"   off {s[3]/1e9:.3f} GB  chr{s[0]}:{s[1]}  {s[2]}", flush=True)
    if not sorted_pos:
        print(f"{ct}: not position-sorted, narrow extraction impossible without an index; stopping", flush=True)
        return 2
    kept = 0
    with gzip.open(out + ".part", "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for chrom, (lo, hi) in REGIONS.items():
            ch = int(chrom)
            def bracket(target):
                a = max([s[3] for s in samples if (s[0], s[1]) < target] or [0])
                b = min([s[3] for s in samples if (s[0], s[1]) > target] or [size])
                for _ in range(14):
                    if b - a < 16_000_000:
                        break
                    m = (a + b) // 2
                    s = rd.key_at(m)
                    if s is None:
                        break
                    if (s[0], s[1]) < target:
                        a = s[3]
                    else:
                        b = m
                return a, b
            a = bracket((ch, lo))[0]
            b = bracket((ch, hi))[1]
            print(f"{ct}: chr{chrom}:{lo}-{hi} -> bytes {a}-{b} ({(b-a)/1e6:.0f} MB)", flush=True)
            pos = a
            carry = b""
            pending = b""
            while pos < b:
                end = min(b, pos + 32_000_000)
                buf = pending + http_range(url, pos, end - 1)
                j = first_member(buf) if not pending else 0
                txt, used = inflate_from(buf, j)
                pending = buf[used:]
                pos = end
                data = carry + txt
                *lines, carry = data.split(b"\n")
                for ln in lines:
                    f = ln.decode().split("\t")
                    if len(f) <= rd.pi:
                        continue
                    if f[rd.fi].split(".")[0] in TARGETS or f[rd.fi] in TARGETS:
                        fo.write(ln.decode() + "\n")
                        kept += 1
            print(f"{ct}: chr{chrom} kept {kept} rows", flush=True)
    os.replace(out + ".part", out)
    print(f"{ct}: done, {kept} rows", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
