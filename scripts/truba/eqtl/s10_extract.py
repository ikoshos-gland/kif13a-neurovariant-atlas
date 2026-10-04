"""Extract the analysis windows from the raw downloads (compute node only, never login).

Outputs (work/):
  kg_eur_chr6.npz            1000G phase 3 genotypes (all 2,504 samples + super-population), chr6 window (GRCh37)
  big40_<pheno>.tsv.gz       BIG40 33k sumstats in the window (GRCh37 as published)
  gwas_udipfa.tsv.gz         GCST90860937 window (build tested later, metadata says GRCh37)
  gwas_sulc.tsv.gz           GCST90095129 window
  gwas_ad_<locus>.tsv.gz     GCST90012877 (AD, Schwartzentruber 2021) positive-control loci
  bryois_<celltype>_<chr>.tsv.gz   Bryois 2022 nominal eQTL rows for the target genes
  bryois_snp_pos_<chr>.tsv.gz      Bryois SNP map (hg38 + hg19 ids, effect/other allele, MAF)
"""
from __future__ import annotations

import gzip
import io
import os
import sys
import time

import numpy as np
import pandas as pd

W = "/arf/scratch/avural/kif13a_eqtl"
RAW = f"{W}/raw"
WORK = f"{W}/work"
os.makedirs(WORK, exist_ok=True)

# KIF13A GRCh37 (Ensembl GRCh37 REST, data/raw/imaging_gwas/ensembl37_KIF13A.json): 17,759,414-17,987,854
K_LO37, K_HI37 = 17_759_414, 17_987_854
WIN_LO, WIN_HI = K_LO37 - 1_150_000, K_HI37 + 1_150_000  # +-1 Mb plus 150 kb margin for build shifts

TARGETS_CHR6 = {"KIF13A", "NUP153", "NHLRC1", "TPMT"}
TARGETS_CHR1 = {"CR1", "NDUFS2", "ADAMTS4"}
AD_LOCI = {"CR1": ("1", 207_302_548, 207_900_071), "ADAMTS4": ("1", 161_014_370, 161_228_800)}


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def extract_kg():
    out = f"{WORK}/kg_eur_chr6.npz"
    if os.path.exists(out):
        log("kg: have")
        return
    # All 2,504 phase 3 samples are kept (super-population label saved) so that LD can be computed
    # for EUR (Bryois, BIG40, UDIP-FA) and as an ancestry mixture (Zeng 2026 is multi-ancestry).
    panel = pd.read_csv(f"{RAW}/kg/panel.txt", sep="\t")
    sp_of = dict(zip(panel["sample"], panel["super_pop"]))
    eur_idx = None
    ids, pos, ref, alt, af, G = [], [], [], [], [], []
    t0 = time.time()
    with gzip.open(f"{RAW}/kg/chr6.vcf.gz", "rt") as fh:
        cols = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                cols = [i for i, s in enumerate(hdr) if s in sp_of]
                spop = np.array([sp_of[hdr[i]] for i in cols])
                eur_idx = np.where(spop == "EUR")[0]
                log(f"kg: {len(cols)} samples, {len(eur_idx)} EUR")
                continue
            # cheap position check before full split
            p = int(line[2:line.index("\t", 2)])
            if p < WIN_LO:
                continue
            if p > WIN_HI:
                break
            f = line.rstrip("\n").split("\t")
            if "," in f[4]:
                continue  # multi-allelic
            gt = f[9:]
            g = np.fromiter(((gt[i - 9][0] != "0") + (gt[i - 9][2] != "0") for i in cols), dtype=np.int8,
                            count=len(cols))
            a_all = g.mean() / 2
            a = g[eur_idx].mean() / 2
            if min(a_all, 1 - a_all) < 0.002:
                continue
            ids.append(f[2]); pos.append(p); ref.append(f[3]); alt.append(f[4]); af.append(a); G.append(g)
    G = np.vstack(G)
    np.savez_compressed(out, rsid=np.array(ids), pos=np.array(pos), ref=np.array(ref), alt=np.array(alt),
                        af=np.array(af), G=G, spop=spop)
    log(f"kg: {len(ids)} variants x {G.shape[1]} samples, {time.time()-t0:.0f}s")


def stream_filter(path, sep, keep_row, header_line=0, usecols=None):
    """Stream a (gzip) text table and keep rows where keep_row(fields)."""
    rows = []
    with gzip.open(path, "rt") as fh:
        hdr = None
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split(sep) if sep else line.split()
            if i == header_line:
                hdr = f
                continue
            r = keep_row(f, hdr)
            if r == "STOP":
                break
            if r:
                rows.append(f)
    return pd.DataFrame(rows, columns=hdr)


def extract_big40():
    for p in ("1769", "1925", "1475"):
        out = f"{WORK}/big40_{p}.tsv.gz"
        if os.path.exists(out):
            continue
        state = {"seen": False}

        def keep(f, hdr):
            if f[0] != "06":
                return "STOP" if state["seen"] else False
            state["seen"] = True
            pos = int(f[2])
            if pos > WIN_HI:
                return "STOP"
            return pos >= WIN_LO
        df = stream_filter(f"{RAW}/big40/{p}.txt.gz", None, keep)
        df.to_csv(out, sep="\t", index=False)
        log(f"big40 {p}: {len(df)} rows; header {list(df.columns)}")


def extract_gwas():
    jobs = [("udipfa", f"{RAW}/gwas/GCST90860937.tsv.gz", "chromosome", "base_pair_location", "6", WIN_LO, WIN_HI),
            ("sulc", f"{RAW}/gwas/GCST90095129_buildGRCh37.tsv.gz", "chromosome", "base_pair_location", "6", WIN_LO, WIN_HI)]
    for name, (c, lo, hi) in AD_LOCI.items():
        jobs.append((f"ad_{name}", f"{RAW}/gwas/GCST90012877_buildGRCh37.tsv.gz", "chromosome",
                     "base_pair_location", c, lo - 50_000, hi + 50_000))
    # Bellenguez 2022 AD (GRCh38) for the Zeng 2026 positive control: whole chr1 cis windows of
    # CR1 (GRCh38 ~207.5 Mb) and ADAMTS4 / NDUFS2 (~161.2 Mb), +-1.3 Mb
    for name, (lo, hi) in {"CR1": (206_200_000, 208_900_000), "ADAMTS4": (159_900_000, 162_500_000)}.items():
        jobs.append((f"bell_{name}", f"{RAW}/gwas/GCST90027158_buildGRCh38.tsv.gz", "chromosome",
                     "base_pair_location", "1", lo, hi))
    for tag, path, ccol, pcol, chrom, lo, hi in jobs:
        out = f"{WORK}/gwas_{tag}.tsv.gz"
        if os.path.exists(out):
            continue
        rows = []
        t0 = time.time()
        for chunk in pd.read_csv(path, sep="\t", chunksize=1_000_000, dtype=str):
            c = chunk[(chunk[ccol] == chrom)]
            if len(c):
                pp = c[pcol].astype(int)
                rows.append(c[(pp >= lo) & (pp <= hi)])
        df = pd.concat(rows) if rows else pd.DataFrame()
        df.to_csv(out, sep="\t", index=False)
        log(f"gwas {tag}: {len(df)} rows ({time.time()-t0:.0f}s)")


def extract_bryois():
    cts = ["Astrocytes", "Endothelial.cells", "Excitatory.neurons", "Inhibitory.neurons", "Microglia",
           "Oligodendrocytes", "OPCs...COPs", "Pericytes", "pb"]
    need_snps = {"6": set(), "1": set()}
    for chrom, targets, cel in (("6", TARGETS_CHR6, cts), ("1", TARGETS_CHR1, ["Oligodendrocytes", "OPCs...COPs"])):
        for ct in cel:
            out = f"{WORK}/bryois_{ct}_{chrom}.tsv.gz"
            if os.path.exists(out):
                d = pd.read_csv(out, sep="\t")
            else:
                rows = []
                with gzip.open(f"{RAW}/bryois/{ct}.{chrom}.gz", "rt") as fh:
                    for line in fh:
                        g = line.split(" ", 1)[0]
                        if g.split("_")[0] in targets:
                            rows.append(line.split())
                d = pd.DataFrame(rows, columns=["gene", "SNP", "dist_tss", "p", "beta"])
                d.to_csv(out, sep="\t", index=False)
                log(f"bryois {ct} chr{chrom}: {len(d)} rows, genes {sorted(d.gene.unique())}")
            need_snps[chrom] |= set(d.SNP)
    for chrom in ("6", "1"):
        out = f"{WORK}/bryois_snp_pos_{chrom}.tsv.gz"
        if os.path.exists(out):
            continue
        rows = []
        with gzip.open(f"{RAW}/bryois/snp_pos.txt.gz", "rt") as fh:
            hdr = fh.readline().split()
            for line in fh:
                s = line.split("\t", 1)[0]
                if s in need_snps[chrom]:
                    rows.append(line.split())
        pd.DataFrame(rows, columns=hdr).to_csv(out, sep="\t", index=False)
        log(f"bryois snp_pos chr{chrom}: {len(rows)} of {len(need_snps[chrom])} needed")


if __name__ == "__main__":
    what = sys.argv[1:] or ["kg", "gwas", "bryois"]  # Zeng: s12_zeng_extract.py; BIG40: s05
    for w in what:
        {"kg": extract_kg, "big40": extract_big40, "gwas": extract_gwas, "bryois": extract_bryois}[w]()
    log("extract done")
