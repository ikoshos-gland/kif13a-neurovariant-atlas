"""Combine per-dataset summaries: KIF13A table, controls, power against the mouse effect,
meta-analysis across independent datasets, ranking and development summaries.

Reads out/summary/<key>/*.tsv written by analyze.py, writes out/summary/combined/*.tsv.
Small and fast; still run as a SLURM job on TRUBA (never on the login node).
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

from stats_lib import meta_analysis, power_slope, power_two_group

MOUSE_DELTA = 0.25   # log2, aged minus young adult, mature oligodendrocyte (Allen aging mouse, donor level)
MOUSE = {"Oligo": dict(delta=0.250, fold=1.19, p=0.0002, cliffs=0.435, n_young=53, n_old=44),
         "OPC": dict(delta=0.131, fold=1.09, p=0.22, cliffs=0.146, n_young=53, n_old=44)}

PRIMARY = ["psychad_aging", "seaad_mtg_oligo", "seaad_mtg_opc", "seaad_dlpfc_oligo",
           "seaad_dlpfc_opc", "wm_all", "postnatal", "velmeshev_velmeshev", "velmeshev_herring"]
# independent adult cohorts for the meta-analysis (SEA-AD MTG and DLPFC share donors: MTG in, DLPFC as sensitivity)
META_SETS = {
    "main": {"Oligo": ["psychad_aging", "seaad_mtg_oligo", "wm_all"],
             "OPC": ["psychad_aging", "seaad_mtg_opc", "wm_all"]},
    "seaad_dlpfc_instead": {"Oligo": ["psychad_aging", "seaad_dlpfc_oligo", "wm_all"],
                            "OPC": ["psychad_aging", "seaad_dlpfc_opc", "wm_all"]},
}
CONTROLS = {
    "OPC": ["PDGFRA", "PLP1", "EBF1"],                      # Seeker 2023 (PMID 37217978): lower PDGFRA, PLP1 in old OPCs
    "Micro": ["CD74", "NEAT1", "HLA-C", "CX3CR1", "P2RY12"],  # Seeker 2023: up CD74, NEAT1, HLA-C; down CX3CR1, P2RY12
    "Oligo": ["C4A", "C4B", "SERPINA3", "B2M", "STAT1"],     # mouse aging oligodendrocytes (Kaya 2022, PMID 36280798)
}


def read(base, key, name):
    p = os.path.join(base, "out", "summary", key, name)
    if not os.path.exists(p):
        return pd.DataFrame()
    try:
        return pd.read_csv(p, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/arf/scratch/avural/kif13a_insan_yas")
    a = ap.parse_args()
    out = os.path.join(a.base, "out", "summary", "combined")
    os.makedirs(out, exist_ok=True)

    tests = pd.concat([read(a.base, k, "age_tests.tsv") for k in PRIMARY], ignore_index=True)
    nulls = pd.concat([read(a.base, k, "null_summary.tsv") for k in PRIMARY], ignore_index=True)
    tests.to_csv(os.path.join(out, "all_age_tests.tsv"), sep="\t", index=False)
    nulls.to_csv(os.path.join(out, "null_summary_all.tsv"), sep="\t", index=False)

    # 1. KIF13A / KIF13B table
    keep_cols = ["dataset", "celltype", "gene", "measure", "n_adult", "age_min", "age_max", "beta_per_decade", "se", "p",
                 "beta_no_depth", "p_no_depth", "covariates", "young_def", "old_def", "tg_n_young", "tg_n_old",
                 "tg_mean_young", "tg_mean_old", "tg_delta", "tg_fold", "tg_mannwhitney_p", "tg_welch_p",
                 "tg_cliffs_delta", "power_twogroup_mouse_effect", "mde80_twogroup", "mde80_slope_per_decade", "resid_sd"]
    k = tests[tests.gene.isin(["KIF13A", "KIF13B"])].reindex(columns=keep_cols)
    nulls = nulls.copy()
    if "frac_genes_p_lt_0.05" in nulls.columns:
        nulls["frac_genes_p_lt_0_05"] = nulls.get("frac_genes_p_lt_0_05", pd.Series(np.nan, index=nulls.index)).fillna(nulls["frac_genes_p_lt_0.05"])
        nulls = nulls.drop(columns=["frac_genes_p_lt_0.05"])
    k = k.merge(nulls.reindex(columns=["dataset", "celltype", "emp_p_all_two_sided", "emp_p_matched_two_sided",
                                       "pct_rank_t_all", "frac_genes_p_lt_0_05"])
                .assign(gene="KIF13A", measure="pb"), on=["dataset", "celltype", "gene", "measure"], how="left")
    k.to_csv(os.path.join(out, "kif13a_age_table.tsv"), sep="\t", index=False)

    # 2. power against the mouse effect, on a common footing
    pa = tests[(tests.dataset == "psychad_aging") & (tests.gene == "KIF13A") & (tests.measure == "pb") & (tests.celltype == "Oligo")]
    gap = float(pa.mean_age_old.iloc[0] - pa.mean_age_young.iloc[0]) if len(pa) else 45.0
    slope_mouse = MOUSE_DELTA / gap * 10  # log2 per decade if the +0.25 between young adult and old accrued linearly
    prow = []
    for key in PRIMARY:
        dv = read(a.base, key, "donor_values.tsv")
        if dv.empty:
            continue
        for ct in ["Oligo", "OPC"]:
            t = tests[(tests.dataset == key) & (tests.celltype == ct) & (tests.gene == "KIF13A") & (tests.measure == "pb")]
            if t.empty:
                continue
            t = t.iloc[0]
            ages = dv[(dv.celltype == ct) & (dv.age >= 20)]["age"].dropna().values / 10.0
            ncov = len(str(t.covariates).split("+")) if isinstance(t.covariates, str) else 0
            pw = power_slope(slope_mouse, t.resid_sd, ages, n_cov=ncov) if len(ages) > ncov + 3 else np.nan
            prow.append(dict(dataset=key, celltype=ct, n_adult=int(t.n_adult), age_range=f"{t.age_min:.0f}-{t.age_max:.0f}",
                             mouse_equiv_slope_per_decade=slope_mouse, observed_slope=t.beta_per_decade, se=t.se,
                             power_slope=pw, mde80_slope_per_decade=t.mde80_slope_per_decade,
                             observed_twogroup_delta=t.get("tg_delta", np.nan),
                             power_twogroup_delta_0_25=t.get("power_twogroup_mouse_effect", np.nan),
                             mde80_twogroup=t.get("mde80_twogroup", np.nan)))
    pw = pd.DataFrame(prow)
    pw.to_csv(os.path.join(out, "power_vs_mouse.tsv"), sep="\t", index=False)

    # 3. meta-analysis of the adult slope (pseudobulk log2CPM, covariate adjusted) and of the young-old contrast
    mrows = []
    for setname, sets in META_SETS.items():
        for ct, keys in sets.items():
            for gene in ["KIF13A", "KIF13B"] + (["PDGFRA"] if ct == "OPC" else []):
                s = tests[(tests.dataset.isin(keys)) & (tests.celltype == ct) & (tests.gene == gene) & (tests.measure == "pb")]
                s = s.dropna(subset=["se"])
                s = s[s.se > 0]
                if len(s) < 2:
                    continue
                m = meta_analysis(s.beta_per_decade.values, s.se.values)
                z = slope_mouse / m["fe_se"]
                power_meta = stats.norm.sf(1.96 - z) + stats.norm.cdf(-1.96 - z)
                mrows.append(dict(set=setname, effect="adult_slope_per_decade", celltype=ct, gene=gene,
                                  datasets=";".join(s.dataset), n_donors=int(s.n_adult.sum()),
                                  estimates=";".join(f"{x:.4f}" for x in s.beta_per_decade),
                                  ses=";".join(f"{x:.4f}" for x in s.se),
                                  power_fe_for_mouse_slope=power_meta, mouse_equiv_slope=slope_mouse, **m))
                # young vs old where both groups exist
                s2 = s.dropna(subset=["tg_delta", "tg_sd_young", "tg_sd_old"])
                if len(s2) >= 2:
                    se2 = np.sqrt(s2.tg_sd_old ** 2 / s2.tg_n_old + s2.tg_sd_young ** 2 / s2.tg_n_young)
                    m2 = meta_analysis(s2.tg_delta.values, se2.values)
                    z2 = MOUSE_DELTA / m2["fe_se"]
                    mrows.append(dict(set=setname, effect="old_minus_young_log2", celltype=ct, gene=gene,
                                      datasets=";".join(s2.dataset), n_donors=int((s2.tg_n_old + s2.tg_n_young).sum()),
                                      estimates=";".join(f"{x:.4f}" for x in s2.tg_delta),
                                      ses=";".join(f"{x:.4f}" for x in se2),
                                      power_fe_for_mouse_slope=stats.norm.sf(1.96 - z2) + stats.norm.cdf(-1.96 - z2),
                                      mouse_equiv_slope=MOUSE_DELTA, **m2))
    meta = pd.DataFrame(mrows)
    meta.to_csv(os.path.join(out, "meta_age.tsv"), sep="\t", index=False)

    # 4. controls
    crow = []
    for ct, genes in CONTROLS.items():
        s = tests[(tests.celltype == ct) & (tests.gene.isin(genes)) & (tests.measure == "pb") & ~tests.dataset.str.contains("no90")]
        crow.append(s.reindex(columns=keep_cols))
    pd.concat(crow).to_csv(os.path.join(out, "controls_positive.tsv"), sep="\t", index=False)
    neg = nulls.copy()
    neg.to_csv(os.path.join(out, "controls_negative_null.tsv"), sep="\t", index=False)
    depth = pd.concat([read(a.base, k, "depth_checks.tsv") for k in PRIMARY], ignore_index=True)
    depth.to_csv(os.path.join(out, "depth_checks_all.tsv"), sep="\t", index=False)

    # 5. ranking (adults)
    rk = pd.concat([read(a.base, k, "ranking.tsv") for k in PRIMARY], ignore_index=True)
    rk.to_csv(os.path.join(out, "ranking_all.tsv"), sep="\t", index=False)

    # 6. development
    dev = pd.concat([read(a.base, k, "development.tsv") for k in PRIMARY], ignore_index=True)
    dev.to_csv(os.path.join(out, "development_all.tsv"), sep="\t", index=False)
    tim = pd.concat([read(a.base, k, "timing.tsv") for k in PRIMARY], ignore_index=True)
    tim.to_csv(os.path.join(out, "timing_all.tsv"), sep="\t", index=False)
    comp = []
    for kk in PRIMARY:
        c = read(a.base, kk, "composition.tsv")
        if not c.empty:
            c.insert(0, "dataset", kk)
            comp.append(c)
    pd.concat(comp, ignore_index=True).to_csv(os.path.join(out, "composition_all.tsv"), sep="\t", index=False)

    pd.set_option("display.width", 250)
    print("mouse-equivalent slope per decade:", round(slope_mouse, 4), "gap years:", round(gap, 1))
    print(pw.round(4).to_string(index=False))
    print(meta[["set", "effect", "celltype", "gene", "datasets", "n_donors", "fe", "fe_se", "fe_p", "re", "re_p", "I2", "power_fe_for_mouse_slope"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
