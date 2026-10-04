"""Donor-level KIF13A age analysis on the pseudobulk outputs of pseudobulk.py.

For each dataset in CONFIG and each cell type (oligodendrocyte, OPC, plus astrocyte and
microglia as context) this writes small tables to out/summary/<key>/:

  donor_values.tsv   one row per donor (x cell type): age, sex, covariates and, for every
                     focal gene, four measures:
                       pb      pseudobulk log2(CPM+1) from summed raw counts
                       cellmean mean over cells of log2(CPM+1), zeros included (the mouse metric)
                       det     fraction of cells with >= 1 count (detection rate)
                       lvl     mean log2(CPM+1) over detecting cells only (level where detected)
  age_tests.tsv      adult regression on age (per decade) with covariates, and the
                     young-vs-old contrast (Mann-Whitney, Welch, Cliff's delta), with power
                     to see the mouse effect (+0.25 log2)
  null_age.tsv       the same adult regression for every expressed gene (pseudobulk),
                     used as the empirical null for KIF13A
  depth_checks.tsv   depth confounding diagnostics
  development.tsv    age-bin summaries over the whole lifespan (if children are present)
  ranking.tsv        cell type ranking with detection / level split and depth-matched detection

Test unit is always the donor. When a dataset has several regions per donor, values are
region-centred and averaged within donor before testing (region acts as a fixed effect).

Run on a compute node: python analyze.py --base /arf/scratch/avural/kif13a_insan_yas --keys wm_all,...
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from stats_lib import (parse_age, ols, two_group, power_two_group, power_slope, lowess_curve)

MOUSE_DELTA = 0.25          # Allen aging mouse, mature oligodendrocyte, log2 units (aged - young)
MIN_CELLS = 20              # minimum nuclei of a cell type for a donor to enter a test
RNG = np.random.default_rng(13)

# Per-dataset settings. Cell type labels are the dataset's own annotation.
CONFIG = {
    "wm_all": dict(
        group=["donor_id", "broad_cell_type", "tissue"], ct="broad_cell_type", region="tissue",
        celltypes={"Oligo": "Oligo", "OPC": "OPC", "Astro": "Astrocyte", "Micro": "Microglia_Macrophages"},
        pmi="tissue_handling_interval", young=(30, 45), old=(60, 200), adult_min=20,
        extra_cov=[], subset=None),
    # SEA-AD: one file per cell type; assay (10x v3 / multiome) handled like a region (fixed effect).
    # Aged donors are 65+, so the adult test is a late-life slope adjusted for AD neuropathology.
    "seaad_mtg_oligo": dict(
        group=["donor_id", "Subclass", "assay"], ct="Subclass", region="assay",
        celltypes={"Oligo": "Oligodendrocyte"}, pmi="PMI", pmi_map={"3.2 to 5.9 hours": 4.55, "5.9 to 8.7 hours": 7.3, "8.7 to 11.4 hours": 10.05},
        young=(20, 45), old=(65, 200), adult_min=20, age90plus=90.0,
        extra_cov=["adnc"], cov_maps={"adnc": ("ADNC", {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3, "Reference": 0})},
        ranking=False),
    "seaad_mtg_opc": dict(
        group=["donor_id", "Subclass", "assay"], ct="Subclass", region="assay",
        celltypes={"OPC": "OPC"}, pmi="PMI", pmi_map={"3.2 to 5.9 hours": 4.55, "5.9 to 8.7 hours": 7.3, "8.7 to 11.4 hours": 10.05},
        young=(20, 45), old=(65, 200), adult_min=20, age90plus=90.0,
        extra_cov=["adnc"], cov_maps={"adnc": ("ADNC", {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3, "Reference": 0})},
        ranking=False),
    # Postnatal multiome (Clarence 2025): 10 donors, 4 regions; descriptive development only.
    "postnatal": dict(
        group=["donor_id", "cell_type", "tissue"], ct="cell_type", region="tissue",
        celltypes={"Oligo": "oligodendrocyte", "OPC": "oligodendrocyte precursor cell", "Astro": "astrocyte", "Micro": "microglial cell"},
        pmi=None, young=(20, 45), old=(60, 200), adult_min=20, extra_cov=[], infant_age=0.5),
}
# PsychAD lifespan cohort (the Nature 2026 lifespan atlas donors, CELLxGENE "Aging_Cohort").
# Age groups follow the atlas: young adulthood 20-39, late adulthood >= 60. No PMI in the h5ad;
# brain bank (Source H = HBCC, M = MSSM) is a covariate. Ages 90+ are coded 90 with a sensitivity run.
CONFIG["psychad_aging"] = dict(
    group=["donor_id", "subclass"], ct="subclass", region=None,
    celltypes={"Oligo": "Oligo", "OPC": "OPC", "Astro": "Astro", "Micro": "Micro"},
    pmi=None, young=(20, 39), old=(60, 200), adult_min=20, age90plus=90.0,
    extra_cov=["source_M"], cov_maps={"source_M": ("Source", {"M": 1, "H": 0})})
# PsychAD HBCC cohort: only donors that are NOT in the lifespan cohort and carry no diagnosis
# ("normal"), so that it could act as an independent replication within the same protocol.
# Checked 4 Oct 2026: all 146 "normal" HBCC donors are already in the lifespan cohort, so this
# leaves zero donors and HBCC is not used (kept here so the check can be rerun).
CONFIG["psychad_hbcc"] = dict(CONFIG["psychad_aging"], donor_filter=("disease", ["normal"]),
                              exclude_donors_from="psychad_aging")
# Velmeshev 2023 integrates four source studies with different chemistries. Each source study is
# analysed on its own (never pooled), development is descriptive.
for _src in ("Velmeshev", "Herring"):
    CONFIG[f"velmeshev_{_src.lower()}"] = dict(
        src="velmeshev", group=["donor_id", "cell_type", "dataset"], ct="cell_type", region=None,
        group_filter=("dataset", [_src]),
        celltypes={"Oligo": "oligodendrocyte", "OPC": "oligodendrocyte precursor cell", "Astro": "astrocyte", "Micro": "microglial cell"},
        pmi=None, young=(20, 39), old=(60, 200), adult_min=20, extra_cov=[], infant_age=0.5, ranking=False)
for _k in ("seaad_dlpfc_oligo", "seaad_dlpfc_opc"):
    CONFIG[_k] = dict(CONFIG[_k.replace("dlpfc", "mtg")])


def load(base, key):
    d = os.path.join(base, "out", key)
    groups = pd.read_csv(os.path.join(d, "groups.tsv"), sep="\t", index_col=0)
    donors = pd.read_csv(os.path.join(d, "donors.tsv"), sep="\t", index_col=0)
    genes = pd.read_csv(os.path.join(d, "genes.tsv"), sep="\t")
    pb = sp.load_npz(os.path.join(d, "pb_counts.npz")).toarray()
    det = sp.load_npz(os.path.join(d, "pb_detect.npz")).toarray()
    fc = np.load(os.path.join(d, "focal_cells.npz"), allow_pickle=True)
    return groups, donors, genes, pb, det, fc


def numeric_hours(x):
    try:
        return float(str(x).split()[0])
    except Exception:
        return np.nan


def per_group_focal(fc, n_groups):
    """Cell-level measures summed per group for each focal gene."""
    g = fc["group"]
    ok = g >= 0
    g = g[ok]
    umi = fc["umi"][ok]
    counts = fc["counts"][ok]
    names = list(fc["genes"])
    cpm = counts / umi[:, None] * 1e6
    l2 = np.log2(cpm + 1)
    detm = counts > 0
    out = {}
    n = np.bincount(g, minlength=n_groups).astype(float)
    for j, s in enumerate(names):
        sum_l2 = np.bincount(g, weights=l2[:, j], minlength=n_groups)
        n_det = np.bincount(g, weights=detm[:, j], minlength=n_groups)
        # Poisson expectation of detection given the group pseudobulk proportion and each cell's depth
        out[s] = dict(sum_l2=sum_l2, n_det=n_det, n=n)
    return out, names, g, umi, counts


def build_donor_values(cfg, groups, donors, genes, pb, fc, ct_label):
    gsel = groups[groups[cfg["ct"]] == ct_label].copy()
    gsel = gsel[gsel.n_cells >= MIN_CELLS]
    idx = gsel.index.values
    focal_stats, names, cell_g, cell_umi, cell_counts = per_group_focal(fc, len(groups))
    sym = genes.symbol.values
    rows = []
    for gi in idx:
        r = dict(group=gi, donor_id=gsel.at[gi, "donor_id"], n_cells=gsel.at[gi, "n_cells"],
                 umi_sum=gsel.at[gi, "umi_sum"], log_genes=gsel.at[gi, "mean_log_genes"],
                 log_umi=gsel.at[gi, "mean_log_umi"])
        if cfg.get("region"):
            r["region"] = gsel.at[gi, cfg["region"]]
        for s in names:
            fs = focal_stats[s]
            col = np.nonzero(sym == s)[0][0]
            cnt = pb[gi, col]
            r[f"{s}|pb"] = np.log2(cnt / gsel.at[gi, "umi_sum"] * 1e6 + 1)
            r[f"{s}|cellmean"] = fs["sum_l2"][gi] / fs["n"][gi]
            r[f"{s}|det"] = fs["n_det"][gi] / fs["n"][gi]
            r[f"{s}|lvl"] = fs["sum_l2"][gi] / fs["n_det"][gi] if fs["n_det"][gi] > 0 else np.nan
        rows.append(r)
    dv = pd.DataFrame(rows)
    # Poisson-expected detection for KIF13A: given each cell's UMI and the group proportion
    if "KIF13A" in names and len(dv):
        j = names.index("KIF13A")
        col = np.nonzero(sym == "KIF13A")[0][0]
        exp_det = []
        for gi in dv.group:
            m = cell_g == gi
            p = pb[gi, col] / groups.at[gi, "umi_sum"]
            exp_det.append(np.mean(1 - np.exp(-p * cell_umi[m])))
        dv["KIF13A|det_expected"] = exp_det
        dv["KIF13A|det_obs_over_exp"] = dv["KIF13A|det"] / dv["KIF13A|det_expected"]
    # donor covariates
    dv["age"] = dv.donor_id.map(lambda d: donor_age(donors, d, cfg))
    dv["age_censored"] = dv.donor_id.map(lambda d: str(donors.at[d, "development_stage"]) == "80 year-old and over stage")
    dv["sex_male"] = dv.donor_id.map(lambda d: 1.0 if str(donors.at[d, "sex"]) == "male" else 0.0)
    if cfg.get("pmi"):
        if cfg.get("pmi_map"):
            dv["pmi"] = dv.donor_id.map(lambda d: cfg["pmi_map"].get(str(donors.at[d, cfg["pmi"]]), np.nan))
        else:
            dv["pmi"] = dv.donor_id.map(lambda d: numeric_hours(donors.at[d, cfg["pmi"]]))
    for c in cfg.get("extra_cov", []):
        if c in cfg.get("cov_maps", {}):
            src, mp = cfg["cov_maps"][c]
            dv[c] = dv.donor_id.map(lambda d: mp.get(str(donors.at[d, src]), np.nan))
        else:
            dv[c] = dv.donor_id.map(lambda d: pd.to_numeric(donors.at[d, c], errors="coerce"))
    return dv


def donor_age(donors, d, cfg):
    if d not in donors.index:
        return np.nan
    if cfg.get("age_col"):
        v = pd.to_numeric(donors.at[d, cfg["age_col"]], errors="coerce")
        if v == v:
            return float(v)
    lab = str(donors.at[d, "development_stage"])
    if lab == "80 year-old and over stage" and cfg.get("age90plus"):
        return float(cfg["age90plus"])
    if lab == "infant stage" and cfg.get("infant_age") is not None:
        return float(cfg["infant_age"])
    return parse_age(lab)


def collapse_regions(dv, value_cols, cov_cols):
    """Region fixed effect: centre every value within region, then average per donor."""
    if "region" not in dv.columns:
        return dv
    x = dv.copy()
    for c in value_cols + ["log_genes", "log_umi"]:
        x[c] = x[c] - x.groupby("region")[c].transform("mean")
    agg = {c: "mean" for c in value_cols + ["log_genes", "log_umi"]}
    agg.update({c: "first" for c in cov_cols})
    agg["n_cells"] = "sum"
    out = x.groupby("donor_id").agg(agg).reset_index()
    out["n_regions"] = x.groupby("donor_id").size().values
    return out


def design(df, covs):
    cols = ["age_dec"] + covs
    X = np.column_stack([np.ones(len(df))] + [df[c].astype(float).values for c in cols])
    return X, ["intercept"] + cols


def adult_tests(dvc, cfg, measures, covs, label):
    out = []
    a = dvc[(dvc.age >= cfg["adult_min"])].copy()
    a = a.dropna(subset=["age"] + covs)
    a["age_dec"] = a.age / 10.0
    if len(a) < 8:
        return out, a
    X, names = design(a, covs)
    young = a[(a.age >= cfg["young"][0]) & (a.age <= cfg["young"][1])]
    old = a[(a.age >= cfg["old"][0]) & (a.age <= cfg["old"][1])]
    span = old.age.mean() - young.age.mean() if len(old) and len(young) else np.nan
    for m in measures:
        y = a[m].values.astype(float)
        ok = ~np.isnan(y)
        if ok.sum() < 8:
            continue
        B, se, t, p, df_, rsd = ols(y[ok], X[ok])
        # model without the depth covariate, for the depth check
        covs_nd = [c for c in covs if c not in ("log_genes", "log_umi")]
        Xn, _ = design(a[ok], covs_nd)
        Bn, sen, tn, pn, _, _ = ols(y[ok], Xn)
        slope_equiv = MOUSE_DELTA / span * 10 if span == span and span > 0 else np.nan
        pw_slope = power_slope(slope_equiv / 10, rsd[0], a.age.values[ok], n_cov=len(covs)) if slope_equiv == slope_equiv else np.nan
        tg = two_group(old[m], young[m])
        pw_tg = (power_two_group(MOUSE_DELTA, tg["sd_old"], tg["sd_young"], tg["n_old"], tg["n_young"])
                 if "sd_old" in tg else np.nan)
        se_tg = (np.sqrt(tg["sd_old"] ** 2 / tg["n_old"] + tg["sd_young"] ** 2 / tg["n_young"]) if "sd_old" in tg else np.nan)
        rho, prho = stats.spearmanr(a.age.values[ok], y[ok])
        out.append(dict(dataset=label[0], celltype=label[1], gene=m.split("|")[0], measure=m.split("|")[1],
                        n_adult=int(ok.sum()), age_min=a.age.min(), age_max=a.age.max(),
                        beta_per_decade=B[1, 0], se=se[1, 0], p=p[1, 0], t=t[1, 0],
                        beta_no_depth=Bn[1, 0], p_no_depth=pn[1, 0], resid_sd=rsd[0],
                        spearman_rho=rho, spearman_p=prho,
                        young_def=f"{cfg['young'][0]}-{cfg['young'][1]}", old_def=f">={cfg['old'][0]}",
                        mean_age_young=young.age.mean(), mean_age_old=old.age.mean(),
                        mouse_equiv_slope_per_decade=slope_equiv, power_slope_mouse_effect=pw_slope,
                        power_twogroup_mouse_effect=pw_tg,
                        mde80_twogroup=2.8 * se_tg if se_tg == se_tg else np.nan,
                        mde80_slope_per_decade=2.8 * se[1, 0],
                        covariates="+".join(covs), **{f"tg_{k}": v for k, v in tg.items()}))
    return out, a


def null_age(a, cfg, groups, genes, pb, covs, ct_label, label, kif_t):
    """Adult age regression for all expressed genes (pseudobulk log2CPM), same design."""
    if len(a) < 8:
        return None, None
    gi = a.group.values if "group" in a.columns else None
    if gi is None:
        return None, None
    lib = groups.loc[gi, "umi_sum"].values
    M = pb[gi, :] / lib[:, None] * 1e6
    expressed = (M.mean(axis=0) >= 1.0) & ((pb[gi, :] > 0).mean(axis=0) >= 0.8)
    Y = np.log2(M[:, expressed] + 1)
    X, names = design(a, covs)
    B, se, t, p, df_, rsd = ols(Y, X)
    res = pd.DataFrame({"symbol": genes.symbol.values[expressed], "mean_log2cpm": Y.mean(axis=0),
                        "beta_per_decade": B[1], "t": t[1], "p": p[1]})
    res.insert(0, "celltype", label[1])
    res.insert(0, "dataset", label[0])
    # empirical position of KIF13A
    k = res[res.symbol == "KIF13A"]
    summ = {}
    if len(k):
        kt = float(k.t.iloc[0])
        km = float(k.mean_log2cpm.iloc[0])
        allg = res[res.symbol != "KIF13A"]
        summ["n_genes_null"] = len(allg)
        summ["kif13a_t"] = kt
        summ["kif13a_beta"] = float(k.beta_per_decade.iloc[0])
        summ["emp_p_all_two_sided"] = float((np.abs(allg.t) >= abs(kt)).mean())
        summ["pct_rank_t_all"] = float((allg.t < kt).mean())
        # 1000 random genes matched on expression (within 0.5 log2 units)
        pool = allg[(allg.mean_log2cpm - km).abs() <= 0.5]
        if len(pool) >= 50:
            samp = pool.sample(n=min(1000, len(pool)), random_state=13)
            summ["n_matched"] = len(samp)
            summ["emp_p_matched_two_sided"] = float((np.abs(samp.t) >= abs(kt)).mean())
            summ["pct_rank_t_matched"] = float((samp.t < kt).mean())
            summ["matched_median_abs_t"] = float(np.median(np.abs(samp.t)))
        summ["all_median_abs_t"] = float(np.median(np.abs(allg.t)))
        summ["frac_genes_p_lt_0.05"] = float((allg.p < 0.05).mean())
    return res, summ


def lifespan_bins(dvc, measures, label):
    bins = [(-1, 0), (0, 1), (1, 6), (6, 12), (12, 20), (20, 40), (40, 60), (60, 200)]
    names = ["prenatal", "0-1", "1-5", "6-11", "12-19", "20-39", "40-59", "60+"]
    rows = []
    for (lo, hi), nm in zip(bins, names):
        s = dvc[(dvc.age >= lo) & (dvc.age < hi)] if lo >= 0 else dvc[dvc.age < 0]
        for m in measures:
            v = s[m].dropna()
            if len(v) == 0:
                continue
            rows.append(dict(dataset=label[0], celltype=label[1], gene=m.split("|")[0], measure=m.split("|")[1],
                             age_bin=nm, n_donors=len(v), median=v.median(), mean=v.mean(),
                             q25=v.quantile(0.25), q75=v.quantile(0.75)))
    return rows


def timing(dvc, measures, label, adult_ref=(20, 40)):
    """Age at which the LOWESS curve (x = log2(age+1)) covers half of its rise or fall from the
    earliest postnatal point to the young-adult level. Descriptive only."""
    rows = []
    s = dvc[dvc.age >= 0].dropna(subset=["age"])
    if len(s) < 15 or s.age.min() > 2:
        return rows
    for m in measures:
        x = np.log2(s.age.values + 1)
        y = s[m].values.astype(float)
        if np.isnan(y).all():
            continue
        grid = np.linspace(x.min(), np.log2(adult_ref[1] + 1), 80)
        gx, gy = lowess_curve(x, y, frac=0.4, grid=grid)
        start = gy[0]
        adult = np.nanmean(y[(s.age >= adult_ref[0]) & (s.age < adult_ref[1])])
        half = start + 0.5 * (adult - start)
        cross = np.nan
        for i in range(1, len(gy)):
            if (gy[i - 1] - half) * (gy[i] - half) <= 0:
                cross = 2 ** gx[i] - 1
                break
        kid = s[s.age < 20]
        rho, p = stats.spearmanr(kid.age, kid[m], nan_policy="omit") if len(kid) >= 6 else (np.nan, np.nan)
        rows.append(dict(dataset=label[0], celltype=label[1], gene=m.split("|")[0], measure=m.split("|")[1],
                         start_value=start, young_adult_value=adult, change=adult - start,
                         age_half_change=cross, n_under20=len(kid), spearman_rho_under20=rho, spearman_p_under20=p))
    return rows


def ranking(cfg, groups, donors, genes, pb, fc, label, adult_min=20):
    """Cell type ranking for KIF13A in adults: donor-averaged measures and depth-matched detection."""
    g = fc["group"]
    ok = g >= 0
    names = list(fc["genes"])
    j = names.index("KIF13A")
    umi = fc["umi"][ok]
    cnt = fc["counts"][ok, j]
    gg = g[ok]
    ct = groups[cfg["ct"]].values[gg]
    donor = groups["donor_id"].values[gg]
    age = pd.Series(donor).map(lambda d: donor_age(donors, d, cfg)).values
    adult = age >= adult_min
    if cfg.get("subset_donors") is not None:
        adult &= np.isin(donor, cfg["subset_donors"])
    umi, cnt, ct, donor = umi[adult], cnt[adult], ct[adult], donor[adult]
    l2 = np.log2(cnt / umi * 1e6 + 1)
    d = pd.DataFrame(dict(ct=ct, donor=donor, umi=umi, cnt=cnt, l2=l2, det=(cnt > 0)))
    # depth bins on the pooled UMI distribution
    edges = np.quantile(d.umi, np.linspace(0, 1, 11))
    edges[0], edges[-1] = -np.inf, np.inf
    d["ubin"] = np.digitize(d.umi, edges[1:-1])
    ref_w = d.ubin.value_counts(normalize=True).sort_index()
    rows = []
    for c, s in d.groupby("ct"):
        per_donor = s.groupby("donor").agg(n=("l2", "size"), cellmean=("l2", "mean"), det=("det", "mean"),
                                          umi=("umi", "median"), cnt=("cnt", "sum"), umis=("umi", "sum"))
        per_donor = per_donor[per_donor.n >= MIN_CELLS]
        if len(per_donor) < 3:
            continue
        lvl = s[s.det].groupby("donor").l2.mean().reindex(per_donor.index)
        pbcpm = per_donor.cnt / per_donor.umis * 1e6
        # depth-standardised detection: detection in each UMI decile, reweighted to the pooled distribution
        dr = s.groupby("ubin").det.mean()
        lv = s[s.det].groupby("ubin").l2.mean()
        w = ref_w.reindex(dr.index).fillna(0)
        det_std = float((dr * w).sum() / w.sum()) if w.sum() > 0 else np.nan
        w2 = ref_w.reindex(lv.index).fillna(0)
        lvl_std = float((lv * w2).sum() / w2.sum()) if w2.sum() > 0 else np.nan
        rows.append(dict(dataset=label, celltype=c, n_donors=len(per_donor), n_cells=int(per_donor.n.sum()),
                         median_umi=float(s.umi.median()),
                         cellmean_log2=per_donor.cellmean.mean(), det_frac=per_donor.det.mean(),
                         lvl_log2_detected=lvl.mean(), cpm_detected=float(np.mean(2 ** lvl - 1)),
                         pb_log2cpm=float(np.log2(pbcpm + 1).mean()),
                         det_frac_depth_std=det_std, lvl_log2_detected_depth_std=lvl_std))
    r = pd.DataFrame(rows).sort_values("cellmean_log2", ascending=False)
    for col in ["cellmean_log2", "det_frac", "lvl_log2_detected", "pb_log2cpm", "det_frac_depth_std"]:
        r[f"rank_{col}"] = r[col].rank(ascending=False).astype(int)
    return r


def run_dataset(base, key, cfg):
    outdir = os.path.join(base, "out", "summary", key)
    os.makedirs(outdir, exist_ok=True)
    groups, donors, genes, pb, det, fc = load(base, cfg.get("src", key))
    if cfg.get("group_filter"):
        col, vals = cfg["group_filter"]
        # groups outside the chosen source study get a cell type label that matches nothing
        groups[cfg["ct"]] = np.where(groups[col].astype(str).isin(vals), groups[cfg["ct"]], "__excluded__")
    keep = pd.Index(donors.index)
    if cfg.get("donor_filter"):
        col, vals = cfg["donor_filter"]
        keep = keep[donors.loc[keep, col].astype(str).isin(vals)]
    if cfg.get("exclude_donors_from"):
        other = pd.read_csv(os.path.join(base, "out", cfg["exclude_donors_from"], "donors.tsv"), sep="	", index_col=0).index
        keep = keep[~keep.isin(other)]
        print(f"{key}: {len(keep)} donors left after excluding donors of {cfg['exclude_donors_from']}", flush=True)
    if cfg.get("donor_filter") or cfg.get("exclude_donors_from"):
        cfg["subset_donors"] = np.array(keep)
    names = list(fc["genes"])
    all_tests, nulls, null_summ, devrows, timrows, depthrows, dvals = [], [], [], [], [], [], []
    for ctkey, ctlab in cfg["celltypes"].items():
        dv = build_donor_values(cfg, groups, donors, genes, pb, fc, ctlab)
        if cfg.get("subset_donors") is not None:
            dv = dv[dv.donor_id.isin(cfg["subset_donors"])]
        if len(dv) == 0:
            continue
        measures = [c for c in dv.columns if "|" in c]
        covs = ["sex_male", "log_genes"] + (["pmi"] if cfg.get("pmi") and dv["pmi"].notna().mean() > 0.9 else []) + cfg.get("extra_cov", [])
        cov_cols = ["age", "age_censored", "sex_male"] + (["pmi"] if "pmi" in dv.columns else []) + cfg.get("extra_cov", [])
        dvc = collapse_regions(dv, measures, cov_cols)
        if "group" not in dvc.columns and "region" not in dv.columns:
            pass
        dvc_out = dvc.copy()
        dvc_out.insert(0, "celltype", ctkey)
        dvals.append(dvc_out)
        tests, a = adult_tests(dvc, cfg, measures, covs, (key, ctkey))
        all_tests += tests
        if dvc.get("age_censored") is not None and dvc.age_censored.astype(bool).any():
            # sensitivity: drop donors whose age is only known as 90+
            t2, _ = adult_tests(dvc[~dvc.age_censored.astype(bool)], cfg, measures, covs, (key + "|no90plus", ctkey))
            all_tests += t2
        # null: only for single-region datasets (group index maps 1:1 to donor)
        if "region" not in dv.columns and len(a):
            ktest = [t for t in tests if t["gene"] == "KIF13A" and t["measure"] == "pb"]
            res, summ = null_age(a, cfg, groups, genes, pb, covs, ctlab, (key, ctkey), ktest[0]["t"] if ktest else np.nan)
            if res is not None:
                nulls.append(res)
                if summ:
                    summ.update(dataset=key, celltype=ctkey)
                    null_summ.append(summ)
        elif "region" in dv.columns:
            # region datasets: null on region-collapsed pseudobulk
            res, summ = null_age_regions(dv, cfg, groups, genes, pb, covs, (key, ctkey))
            if res is not None:
                nulls.append(res)
                summ.update(dataset=key, celltype=ctkey)
                null_summ.append(summ)
        # depth diagnostics
        if len(a) >= 8:
            rho, p = stats.spearmanr(a.age, a.log_genes)
            rho2, p2 = stats.spearmanr(a.age, a.log_umi)
            depthrows.append(dict(dataset=key, celltype=ctkey, n=len(a), rho_age_loggenes=rho, p_age_loggenes=p,
                                  rho_age_logumi=rho2, p_age_logumi=p2))
        devrows += lifespan_bins(dvc, measures, (key, ctkey))
        timrows += timing(dvc, [m for m in measures if m.split("|")[1] in ("pb", "cellmean")], (key, ctkey))
    pd.concat(dvals).to_csv(os.path.join(outdir, "donor_values.tsv"), sep="\t", index=False)
    pd.DataFrame(all_tests).to_csv(os.path.join(outdir, "age_tests.tsv"), sep="\t", index=False)
    if nulls:
        pd.concat(nulls).to_csv(os.path.join(outdir, "null_age.tsv.gz"), sep="\t", index=False)
    pd.DataFrame(null_summ).to_csv(os.path.join(outdir, "null_summary.tsv"), sep="\t", index=False)
    pd.DataFrame(depthrows).to_csv(os.path.join(outdir, "depth_checks.tsv"), sep="\t", index=False)
    pd.DataFrame(devrows).to_csv(os.path.join(outdir, "development.tsv"), sep="\t", index=False)
    pd.DataFrame(timrows).to_csv(os.path.join(outdir, "timing.tsv"), sep="\t", index=False)
    if cfg.get("ranking", True):
        r = ranking(cfg, groups, donors, genes, pb, fc, key)
        r.to_csv(os.path.join(outdir, "ranking.tsv"), sep="\t", index=False)
    composition(cfg, groups, donors, key, outdir)
    print(f"{key}: done", flush=True)


def null_age_regions(dv, cfg, groups, genes, pb, covs, label):
    """Null for multi-region datasets: region-centred pseudobulk log2CPM averaged per donor."""
    gi = dv.group.values
    lib = groups.loc[gi, "umi_sum"].values
    M = pb[gi, :] / lib[:, None] * 1e6
    expressed = (M.mean(axis=0) >= 1.0) & ((pb[gi, :] > 0).mean(axis=0) >= 0.8)
    Y = np.log2(M[:, expressed] + 1)
    Yd = pd.DataFrame(Y)
    Yd["region"] = dv.region.values
    Yd["donor"] = dv.donor_id.values
    vals = Yd.drop(columns=["donor"]).groupby("region").transform(lambda c: c - c.mean())
    vals["donor"] = Yd.donor
    Yc = vals.groupby("donor").mean()
    first_cols = ["age", "age_censored"] + [c for c in covs if c not in ("log_genes", "log_umi")]
    meta = dv.groupby("donor_id").agg(**{c: (c, "first") for c in first_cols}, log_genes=("log_genes", "mean"))
    meta = meta.loc[Yc.index]
    lg = dv.copy()
    lg["log_genes"] = lg.log_genes - lg.groupby("region").log_genes.transform("mean")
    meta["log_genes"] = lg.groupby("donor_id").log_genes.mean().loc[Yc.index]
    a = meta[meta.age >= cfg["adult_min"]].dropna(subset=["age"] + covs).copy()
    a["age_dec"] = a.age / 10
    if len(a) < 8:  # too few adult donors for a regression (postnatal set has 4)
        return None, None
    Ya = Yc.loc[a.index].values
    X, _ = design(a, covs)
    B, se, t, p, df_, rsd = ols(Ya, X)
    sym = genes.symbol.values[expressed]
    res = pd.DataFrame({"symbol": sym, "mean_log2cpm": Y.mean(axis=0), "beta_per_decade": B[1], "t": t[1], "p": p[1]})
    res.insert(0, "celltype", label[1])
    res.insert(0, "dataset", label[0])
    summ = {}
    k = res[res.symbol == "KIF13A"]
    if len(k):
        kt = float(k.t.iloc[0])
        km = float(k.mean_log2cpm.iloc[0])
        allg = res[res.symbol != "KIF13A"]
        summ.update(n_genes_null=len(allg), kif13a_t=kt, kif13a_beta=float(k.beta_per_decade.iloc[0]),
                    emp_p_all_two_sided=float((np.abs(allg.t) >= abs(kt)).mean()),
                    pct_rank_t_all=float((allg.t < kt).mean()),
                    all_median_abs_t=float(np.median(np.abs(allg.t))),
                    frac_genes_p_lt_0_05=float((allg.p < 0.05).mean()))
        pool = allg[(allg.mean_log2cpm - km).abs() <= 0.5]
        if len(pool) >= 50:
            samp = pool.sample(n=min(1000, len(pool)), random_state=13)
            summ.update(n_matched=len(samp), emp_p_matched_two_sided=float((np.abs(samp.t) >= abs(kt)).mean()),
                        pct_rank_t_matched=float((samp.t < kt).mean()))
    return res, summ


def composition(cfg, groups, donors, key, outdir):
    """Fraction of nuclei per cell type per donor (development positive control)."""
    g = groups[groups[cfg["ct"]] != "__excluded__"].copy()
    if cfg.get("subset_donors") is not None:
        g = g[g.donor_id.isin(cfg["subset_donors"])]
    keys = ["donor_id"] + ([cfg["region"]] if cfg.get("region") else [])
    tot = g.groupby(keys).n_cells.transform("sum")
    g["frac"] = g.n_cells / tot
    g = g[g[cfg["ct"]].isin(cfg["celltypes"].values())]
    comp = g.groupby(["donor_id", cfg["ct"]]).frac.mean().unstack(fill_value=0)
    comp["age"] = [donor_age(donors, d, cfg) for d in comp.index]
    comp.to_csv(os.path.join(outdir, "composition.tsv"), sep="\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/arf/scratch/avural/kif13a_insan_yas")
    ap.add_argument("--keys", default=",".join(CONFIG))
    a = ap.parse_args()
    for key in a.keys.replace(":", ",").split(","):
        run_dataset(a.base, key, dict(CONFIG[key]))


if __name__ == "__main__":
    main()
