"""Point liftOver (UCSC chain format) for one chromosome, hg38 -> hg19.

Only ungapped chain blocks are used (a position inside a gap does not map, as with UCSC liftOver).
Validated in s20_coloc.py build stage against Bryois snp_pos (which carries both hg38 and hg19
ids, lifted by the authors with the same UCSC chain).
"""
from __future__ import annotations

import bisect
import gzip

import numpy as np


class Lift:
    def __init__(self, chain_path: str, chrom: str):
        tchrom = chrom if chrom.startswith("chr") else "chr" + chrom
        starts, ends, offs, qchr, qstrand = [], [], [], [], []
        with gzip.open(chain_path, "rt") as fh:
            keep = False
            for line in fh:
                if line.startswith("chain"):
                    f = line.split()
                    tName, tStart = f[2], int(f[5])
                    qName, qSize, qStrand, qStart = f[7], int(f[8]), f[9], int(f[10])
                    keep = tName == tchrom and f[4] == "+"
                    t, q = tStart, qStart
                    continue
                if not keep:
                    continue
                f = line.split()
                if not f:
                    keep = False
                    continue
                size = int(f[0])
                starts.append(t)
                ends.append(t + size)
                if qStrand == "+":
                    offs.append(q - t)
                else:  # store negative-strand blocks but mark them; they are not used for SNP mapping
                    offs.append(None)
                qchr.append(qName)
                qstrand.append(qStrand)
                if len(f) == 3:
                    t += size + int(f[1])
                    q += size + int(f[2])
        order = np.argsort(starts)
        self.starts = [starts[i] for i in order]
        self.ends = [ends[i] for i in order]
        self.offs = [offs[i] for i in order]
        self.qchr = [qchr[i] for i in order]

    def map(self, pos1: int):
        """1-based position on the source build -> (chrom, 1-based position) or None."""
        p0 = pos1 - 1
        i = bisect.bisect_right(self.starts, p0) - 1
        hits = []
        # blocks can overlap between chains; check a few neighbours
        for j in range(max(0, i - 3), min(len(self.starts), i + 4)):
            if self.starts[j] <= p0 < self.ends[j] and self.offs[j] is not None:
                hits.append((self.qchr[j], p0 + self.offs[j] + 1))
        if len(set(hits)) == 1:
            return hits[0]
        return None
