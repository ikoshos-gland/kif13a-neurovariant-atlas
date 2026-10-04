"""Donor x cell type pseudobulk from a CELLxGENE h5ad, streamed from disk.

Runs on a TRUBA compute node (never on the login node). Reads the RAW count
matrix (raw/X when present, otherwise X), checks that the values are integers,
and writes, per group (donor x cell type x optional extra columns):

  pb_counts.npz      groups x genes, summed raw counts (float64)
  pb_detect.npz      groups x genes, number of cells with count > 0 (int32)
  groups.tsv         group keys, n_cells, summed UMI, mean genes per cell,
                     mean log(genes per cell), mean log(UMI per cell)
  genes.tsv          Ensembl id, symbol, per-gene total counts
  donors.tsv         every obs column that is constant within each donor
  focal_cells.npz    per cell: group index, UMI total, genes detected, and raw
                     counts of the focal genes (KIF13A, paralogs, controls)
  obs_summary.txt    obs columns, dtypes, number of levels (for the record)

log2(CPM + 1) is recomputed downstream from these raw sums. Nothing here uses
the normalised X of the h5ad, because CELLxGENE X can be in a different unit
(the mouse spinal cord analysis fell into exactly that trap).

Usage:
  python pseudobulk.py --h5ad raw/psychad_aging.h5ad --out out/psychad_aging \
      --group donor_id,cell_type [--chunk 20000] [--workers 6]
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

FOCAL = [
    # target and paralogs
    "KIF13A", "KIF13B", "KIF1A", "KIF1B", "KIF1C", "KIF16B",
    # OPC age controls, human white matter (Seeker 2023, PMID 37217978): lower in old OPCs
    "PDGFRA", "PLP1", "EBF1",
    # microglia age controls, same source: up CD74, NEAT1, HLA-C; down CX3CR1, P2RY12
    "CD74", "NEAT1", "HLA-C", "CX3CR1", "P2RY12",
    # mouse aging-related oligodendrocyte markers (Kaya 2022, PMID 36280798), human orthologs
    "C4A", "C4B", "SERPINA3", "B2M", "STAT1",
    # myelination program and lineage markers (reference curves, not controls)
    "MBP", "MOBP", "MOG", "MAG", "OPALIN", "CNP", "CSPG4", "VCAN", "BCAS1", "GPR37",
    # housekeeping
    "ACTB", "GAPDH",
]


def read_obs(f):
    obs = ad.io.read_elem(f["obs"]) if hasattr(ad, "io") else ad.experimental.read_elem(f["obs"])
    return obs


def matrix_group(f):
    if "raw" in f and "X" in f["raw"]:
        return "raw/X", "raw/var"
    return "X", "var"


def read_var(f, var_path):
    var = ad.io.read_elem(f[var_path]) if hasattr(ad, "io") else ad.experimental.read_elem(f[var_path])
    return var


def is_integer_sample(f, xpath, n=2_000_000):
    d = f[xpath]["data"]
    m = min(n, d.shape[0])
    x = d[:m]
    return bool(np.all(np.abs(x - np.round(x)) < 1e-6)), float(x.min()), float(x.max()), str(x.dtype)


def worker(args):
    (h5ad, xpath, start, stop, chunk, gidx_path, n_groups, n_genes, focal_idx, tmpdir, wid) = args
    gidx = np.load(gidx_path, mmap_mode="r")
    pb = np.zeros((n_groups, n_genes), dtype=np.float64)
    det = np.zeros((n_groups, n_genes), dtype=np.int32)
    tot_list, ng_list, foc_list = [], [], []
    with h5py.File(h5ad, "r") as f:
        g = f[xpath]
        indptr = g["indptr"]
        for a in range(start, stop, chunk):
            b = min(stop, a + chunk)
            ip = indptr[a:b + 1].astype(np.int64)
            lo, hi = ip[0], ip[-1]
            data = g["data"][lo:hi].astype(np.float64)
            ind = g["indices"][lo:hi]
            X = sp.csr_matrix((data, ind, ip - lo), shape=(b - a, n_genes))
            tot = np.asarray(X.sum(axis=1)).ravel()
            ng = np.diff(X.indptr)
            gi = np.asarray(gidx[a:b])
            keep = gi >= 0
            rows = gi[keep]
            cols = np.nonzero(keep)[0]
            G = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_groups, b - a))
            r = (G @ X).tocoo()
            r.sum_duplicates()
            pb[r.row, r.col] += r.data
            Xb = X.copy()
            Xb.data = (Xb.data > 0).astype(np.float64)
            r = (G @ Xb).tocoo()
            r.sum_duplicates()
            det[r.row, r.col] += np.rint(r.data).astype(np.int32)
            foc = X[:, focal_idx].toarray().astype(np.int32) if len(focal_idx) else np.zeros((b - a, 0), np.int32)
            tot_list.append(tot)
            ng_list.append(ng)
            foc_list.append(foc)
    np.save(os.path.join(tmpdir, f"pb_{wid}.npy"), pb)
    np.save(os.path.join(tmpdir, f"det_{wid}.npy"), det)
    np.savez(os.path.join(tmpdir, f"cells_{wid}.npz"), start=start, stop=stop,
             tot=np.concatenate(tot_list), ng=np.concatenate(ng_list), foc=np.vstack(foc_list))
    return wid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--group", default="donor_id,cell_type")
    ap.add_argument("--chunk", type=int, default=20000)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--inspect-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    with h5py.File(a.h5ad, "r") as f:
        obs = read_obs(f)
        xpath, varpath = matrix_group(f)
        var = read_var(f, varpath)
        enc = dict(f[xpath].attrs)
        shape = tuple(enc.get("shape", (obs.shape[0], var.shape[0])))
        isint, xmin, xmax, dt = is_integer_sample(f, xpath)
        # also look at the main X for the record
        xinfo = None
        if xpath != "X" and "X" in f:
            try:
                xinfo = is_integer_sample(f, "X", 200000)
            except Exception as e:  # dense X or other layout
                xinfo = str(e)
    with open(os.path.join(a.out, "obs_summary.txt"), "w") as fh:
        fh.write(f"file {a.h5ad}\nn_obs {obs.shape[0]} n_var {var.shape[0]}\n")
        fh.write(f"count matrix {xpath} encoding {enc.get('encoding-type')} shape {shape}\n")
        fh.write(f"count sample integer={isint} min={xmin} max={xmax} dtype={dt}\n")
        fh.write(f"main X sample (integer, min, max, dtype) = {xinfo}\n")
        fh.write(f"var columns {list(var.columns)}\n\n")
        for c in obs.columns:
            s = obs[c]
            nun = s.nunique(dropna=False)
            ex = list(pd.Series(s.unique()[:8]).astype(str))
            fh.write(f"{c}\t{s.dtype}\t{nun}\t{ex}\n")
    print(open(os.path.join(a.out, "obs_summary.txt")).read(), flush=True)
    if a.inspect_only:
        return
    if not isint:
        sys.exit("Count matrix is not integer: refusing to treat it as raw counts")
    if enc.get("encoding-type") != "csr_matrix":
        sys.exit(f"Expected csr_matrix, found {enc.get('encoding-type')}")

    gcols = a.group.replace(":", ",").split(",")  # ":" allowed because sbatch --export splits on commas
    key = obs[gcols].astype(str).agg("||".join, axis=1)
    bad = obs[gcols].isna().any(axis=1)
    codes, uniq = pd.factorize(key)
    codes = np.where(bad.values, -1, codes).astype(np.int64)
    n_groups = len(uniq)
    gidx_path = os.path.join(a.out, "_gidx.npy")
    np.save(gidx_path, codes)

    sym_col = "feature_name" if "feature_name" in var.columns else None
    symbols = var[sym_col].astype(str).values if sym_col else var.index.astype(str).values
    focal_idx, focal_names = [], []
    for s in FOCAL:
        hit = np.nonzero(symbols == s)[0]
        if len(hit) == 1:
            focal_idx.append(int(hit[0]))
            focal_names.append(s)
        else:
            print(f"focal gene {s}: {len(hit)} matches, skipped", flush=True)
    n_cells, n_genes = obs.shape[0], var.shape[0]

    # split cells into contiguous ranges, one per task
    n_tasks = max(a.workers * 3, 1)
    bounds = np.linspace(0, n_cells, n_tasks + 1).astype(int)
    tmpdir = os.path.join(a.out, "_tmp")
    os.makedirs(tmpdir, exist_ok=True)
    tasks = [(a.h5ad, xpath, int(bounds[i]), int(bounds[i + 1]), a.chunk, gidx_path, n_groups,
              n_genes, focal_idx, tmpdir, i) for i in range(n_tasks) if bounds[i + 1] > bounds[i]]
    print(f"{n_cells} cells, {n_genes} genes, {n_groups} groups, {len(tasks)} tasks, {a.workers} workers", flush=True)
    pb = np.zeros((n_groups, n_genes), dtype=np.float64)
    det = np.zeros((n_groups, n_genes), dtype=np.int64)
    tot = np.zeros(n_cells)
    ng = np.zeros(n_cells, dtype=np.int64)
    foc = np.zeros((n_cells, len(focal_idx)), dtype=np.int32)
    with Pool(a.workers) as pool:
        for wid in pool.imap_unordered(worker, tasks):
            pb += np.load(os.path.join(tmpdir, f"pb_{wid}.npy"))
            det += np.load(os.path.join(tmpdir, f"det_{wid}.npy"))
            z = np.load(os.path.join(tmpdir, f"cells_{wid}.npz"))
            s0, s1 = int(z["start"]), int(z["stop"])
            tot[s0:s1], ng[s0:s1], foc[s0:s1] = z["tot"], z["ng"], z["foc"]
            for fn in (f"pb_{wid}.npy", f"det_{wid}.npy", f"cells_{wid}.npz"):
                os.remove(os.path.join(tmpdir, fn))
            print(f"task {wid} merged, {time.time() - t0:.0f}s", flush=True)
    os.rmdir(tmpdir)

    # group table
    gt = pd.DataFrame([u.split("||") for u in uniq], columns=gcols)
    valid = codes >= 0
    cnt = np.bincount(codes[valid], minlength=n_groups)
    gt["n_cells"] = cnt
    gt["umi_sum"] = np.bincount(codes[valid], weights=tot[valid], minlength=n_groups)
    gt["mean_genes_per_cell"] = np.bincount(codes[valid], weights=ng[valid], minlength=n_groups) / np.maximum(cnt, 1)
    gt["mean_log_genes"] = np.bincount(codes[valid], weights=np.log(np.maximum(ng[valid], 1)), minlength=n_groups) / np.maximum(cnt, 1)
    gt["mean_log_umi"] = np.bincount(codes[valid], weights=np.log(np.maximum(tot[valid], 1)), minlength=n_groups) / np.maximum(cnt, 1)
    gt.index.name = "group"
    gt.to_csv(os.path.join(a.out, "groups.tsv"), sep="\t")

    genes = pd.DataFrame({"ensembl": var.index.astype(str), "symbol": symbols, "total_counts": pb.sum(axis=0)})
    genes.to_csv(os.path.join(a.out, "genes.tsv"), sep="\t", index=False)
    sp.save_npz(os.path.join(a.out, "pb_counts.npz"), sp.csr_matrix(pb))
    sp.save_npz(os.path.join(a.out, "pb_detect.npz"), sp.csr_matrix(det.astype(np.int32)))

    # donor-level metadata: columns constant within every donor
    if "donor_id" in obs.columns:
        dcols = []
        grp = obs.groupby("donor_id", observed=True)
        for c in obs.columns:
            try:
                if (grp[c].nunique(dropna=False) <= 1).all():
                    dcols.append(c)
            except TypeError:
                pass
        donors = grp[dcols].first() if dcols else pd.DataFrame(index=grp.size().index)
        donors["n_cells_total"] = grp.size()
        donors.to_csv(os.path.join(a.out, "donors.tsv"), sep="\t")

    np.savez_compressed(os.path.join(a.out, "focal_cells.npz"), group=codes, umi=tot, n_genes=ng,
                        counts=foc, genes=np.array(focal_names))
    os.remove(gidx_path)
    meta = {"h5ad": a.h5ad, "count_matrix": xpath, "group_cols": gcols, "n_cells": int(n_cells),
            "n_genes": int(n_genes), "n_groups": int(n_groups), "focal_genes": focal_names,
            "seconds": round(time.time() - t0, 1)}
    json.dump(meta, open(os.path.join(a.out, "meta.json"), "w"), indent=1)
    print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
