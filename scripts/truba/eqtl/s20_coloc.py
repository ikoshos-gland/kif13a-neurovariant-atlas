"""Cell-type eQTL x white-matter GWAS colocalisation at the KIF13A locus (compute node only).

Stages (python s20_coloc.py <stage> ...):
  build    coordinate-build test of every dataset against 1000G phase 3 (GRCh37) by rsID,
           Zeng start/end convention against Bryois hg38 ids, and the python liftOver
  posctrl  positive controls: Bryois 2022 Table S5 (AD x Oligodendrocytes, CR1/NDUFS2/ADAMTS4)
           and Zeng 2026 coloc_class.tsv (AD x Oligo, CR1/ADAMTS4/NDUFS2)
  main     coloc.abf for every eQTL dataset x gene x GWAS trait (+ p12 sensitivity)
  lookup   eQTL statistics at the haplotype tag variants, implied effective N per dataset
  cond     stepwise conditional eQTL signals (1000G LD) and coloc on each conditional signal
  power    simulation: chance of PP.H4 >= 0.5 / 0.8 if the haplotype truly changed KIF13A
  export   R input files for the coloc / susieR reference run (s30_validate.R)
Outputs: /arf/scratch/avural/kif13a_eqtl/results/
Merging is by rsID plus unordered allele pair, the window is defined on GRCh37 GWAS positions.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coloc_py as cp  # noqa: E402

W = "/arf/scratch/avural/kif13a_eqtl"
WORK = f"{W}/work"
RES = f"{W}/results"
os.makedirs(RES, exist_ok=True)

K_LO37, K_HI37 = 17_759_414, 17_987_854  # KIF13A GRCh37 (Ensembl GRCh37 REST)
WLO, WHI = K_LO37 - 1_000_000, K_HI37 + 1_000_000  # analysis window: gene +-1 Mb (GRCh37)

TAGS = ["rs76091407", "rs144963738", "rs75428660", "rs115824472", "rs594115", "rs1924466", "rs143422697",
        "rs73373199", "rs112216920", "rs73724776", "rs16880094", "rs16880195"]
TAG_ROLE = {"rs76091407": "UDIP-FA GWS lead (GCST90860937)",
            "rs144963738": "BIG40 ak madde sinifi oncusu (L3 sol on korona radiata)",
            "rs75428660": "sulkus derinligi GWS (GCST90095129)",
            "rs115824472": "MetaBrain korteks KIF13A eQTL sinyal 3 oncusu",
            "rs594115": "MetaBrain korteks KIF13A eQTL sinyal 1 oncusu",
            "rs1924466": "MetaBrain korteks KIF13A eQTL sinyal 2 oncusu",
            "rs143422697": "BIG40 tek GWS hiti (gri-ak kontrasti), ayri haplotip",
            "rs73373199": "Bryois setindeki sulkus derinligi oncusu (haplotip vekili?)",
            "rs112216920": "Bryois setindeki BIG40 L3 oncusu (haplotip test edilmemis)",
            "rs73724776": "Bryois setindeki UDIP-FA oncusu (haplotip test edilmemis)",
            "rs16880094": "Zeng EN kosullu sinyal 3 (LD paneline duyarli)",
            "rs16880195": "SingleBrain EN kosullu sinyal 2"}
GENES = ["KIF13A", "NUP153", "NHLRC1", "TPMT"]
BRYOIS_CT = {"Oligodendrocytes": "Oligo", "OPCs...COPs": "OPC", "Astrocytes": "Astro",
             "Excitatory.neurons": "EN", "Inhibitory.neurons": "IN", "Microglia": "Immune",
             "Endothelial.cells": "Endo", "Pericytes": "Mural", "pb": "pseudobulk"}
ZENG_CT = ["Oligo", "OPC", "Astro", "EN", "IN", "Immune", "Endo", "Mural"]
SBRAIN_CT = {"OD": "Oligo", "OPC": "OPC", "Ast": "Astro", "Ext": "EN", "IN": "IN", "MG": "Immune", "End": "Endo"}
ENSG = {"ENSG00000137177": "KIF13A", "ENSG00000124789": "NUP153", "ENSG00000187566": "NHLRC1",
        "ENSG00000137364": "TPMT"}


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def ukey(rsid, a, b):
    a, b = str(a).upper(), str(b).upper()
    return f"{rsid}:{min(a, b)}:{max(a, b)}"


# ----------------------------------------------------------------------------------------------
# loaders: GWAS -> rsid, pos (GRCh37), a_eff, a_oth, beta, se, p, n, eaf ; eQTL -> rsid, a_eff,
# a_oth, beta, se, p, (maf), symbol
# ----------------------------------------------------------------------------------------------

_KG = None


def pkey(pos, a, b):
    a, b = str(a).upper(), str(b).upper()
    return f"{int(pos)}:{min(a, b)}:{max(a, b)}"


def load_kg():
    """1000G phase 3 chr6 window. The v5b release VCF carries no rsIDs (ID = '.'), so rsIDs are
    attached through position + allele pair (GRCh37) from the BIG40 and UDIP-FA tables."""
    global _KG
    if _KG is None:
        z = np.load(f"{WORK}/kg_eur_chr6.npz", allow_pickle=False)
        kg = pd.DataFrame({"pos37": z["pos"], "ref": z["ref"], "alt": z["alt"], "af": z["af"]})
        kg["pkey"] = [pkey(p, a, b) for p, a, b in zip(kg.pos37, kg.ref, kg.alt)]
        rs = {}
        for t in (load_big40("1769"), load_udipfa(), load_sulc()):
            for r, p, a, b in zip(t.rsid, t.pos, t.a_eff, t.a_oth):
                if isinstance(r, str) and r.startswith("rs"):
                    rs.setdefault(pkey(p, a, b), r)
        kg["rsid"] = kg.pkey.map(rs)
        kg["key"] = [ukey(r, a, b) if isinstance(r, str) else None for r, a, b in zip(kg.rsid, kg.ref, kg.alt)]
        _KG = (kg, z["G"], z["spop"])
    return _KG


def load_big40(pheno):
    d = pd.read_csv(f"{WORK}/big40_{pheno}.tsv.gz", sep="\t")
    d.columns = ["chr", "rsid", "pos", "a1", "a2", "beta", "se", "nlog10p"]
    out = pd.DataFrame({"rsid": d.rsid, "pos": d.pos.astype(int), "a_eff": d.a2, "a_oth": d.a1,
                        "beta": d.beta.astype(float), "se": d.se.astype(float),
                        "p": 10 ** (-d.nlog10p.astype(float))})
    out["n"] = 31356 if pheno != "1925" else 31355
    out["eaf"] = np.nan
    return out


def load_udipfa():
    d = pd.read_csv(f"{WORK}/gwas_udipfa.tsv.gz", sep="\t", na_values=["#NA", "NA"])
    return pd.DataFrame({"rsid": d.rsid, "pos": d.base_pair_location.astype(int), "a_eff": d.effect_allele,
                         "a_oth": d.other_allele, "beta": np.nan, "se": np.nan, "p": d.p_value.astype(float),
                         "eaf": d.effect_allele_frequency.astype(float), "n": d.n.astype(float)})


def load_sulc():
    d = pd.read_csv(f"{WORK}/gwas_sulc.tsv.gz", sep="\t")
    return pd.DataFrame({"rsid": d.variant_id, "pos": d.base_pair_location.astype(int), "a_eff": d.effect_allele,
                         "a_oth": d.other_allele, "beta": np.nan, "se": np.nan, "p": d.p_value.astype(float),
                         "z": d.Z.astype(float), "n": d.N.astype(float), "eaf": np.nan})


def load_bryois(ct, chrom="6"):
    e = pd.read_csv(f"{WORK}/bryois_{ct}_{chrom}.tsv.gz", sep="\t")
    sp = pd.read_csv(f"{WORK}/bryois_snp_pos_{chrom}.tsv.gz", sep="\t")
    e = e.merge(sp, on="SNP", how="inner")
    e["pos37"] = e.SNP_id_hg19.str.split(":").str[1].astype(float)
    e["pos38"] = e.SNP_id_hg38.str.split(":").str[1].astype(float)
    e["symbol"] = e.gene.str.split("_").str[0]
    e["p"] = e.p.astype(float)
    e["beta"] = e.beta.astype(float)
    with np.errstate(divide="ignore"):
        e["se"] = np.abs(e.beta / stats.norm.ppf(e.p / 2))  # exactly as Bryois 04_coloc.Rmd
    e = e.rename(columns={"SNP": "rsid", "effect_allele": "a_eff", "other_allele": "a_oth", "MAF": "maf"})
    return e


def load_zeng(cls):
    d = pd.read_csv(f"{WORK}/zeng_{cls}.tsv.gz", sep="\t")
    d = d.rename(columns={"variant": "rsid", "assessed_allele": "a_eff", "nonassessed_allele": "a_oth",
                          "gene": "symbol"})
    d["p"] = 2 * stats.norm.sf(np.abs(d.z.astype(float)))
    d["maf"] = np.nan
    d["pos38"] = d.end.astype(float)
    return d


def load_sbrain(ct):
    """SingleBrain (Zenodo 14908182): fixed-effect meta beta/sd across the 4 cohorts; rows whose
    fixed_sd is 0 are cohort-specific allele codings without a meta estimate and are dropped."""
    d = pd.read_csv(f"{WORK}/singlebrain_{ct}.tsv.gz", sep="	")
    d = d[d.fixed_sd > 0].copy()
    d["symbol"] = d.feature.str.split(".").str[0].map(ENSG)
    d["a_eff"] = d.Allele
    d["a_oth"] = np.where(d.Allele == d.alt, d.ref, d.alt)
    d["beta"] = d.fixed_beta
    d["se"] = d.fixed_sd
    d["p"] = d.Fixed_P.astype(float)
    d["rsid"] = np.where(d.variant_id.str.startswith("rs"), d.variant_id, None)
    d["maf"] = np.nan
    d["pos38"] = d.pos
    return d[d.rsid.notna()]


def gwas_sets():
    return {"BIG40 L3 sol on korona radiata (UKB 25319)": (load_big40("1769"), "beta"),
            "BIG40 ICVF sol on korona radiata (UKB 25367)": (load_big40("1925"), "beta"),
            "BIG40 FA sol on korona radiata (UKB 25079)": (load_big40("1475"), "beta"),
            "UDIP-FA cok degiskenli (GCST90860937)": (load_udipfa(), "p"),
            "Sulkus derinligi (GCST90095129)": (load_sulc(), "p")}


def eqtl_sets():
    out = {}
    for cls in ZENG_CT:
        if os.path.exists(f"{WORK}/zeng_{cls}.tsv.gz"):
            out[("Zeng 2026", cls)] = load_zeng(cls)
    for ct, lab in BRYOIS_CT.items():
        if os.path.exists(f"{WORK}/bryois_{ct}_6.tsv.gz"):
            out[("Bryois 2022", lab)] = load_bryois(ct)
    for ct, lab in SBRAIN_CT.items():
        if os.path.exists(f"{WORK}/singlebrain_{ct}.tsv.gz"):
            out[("SingleBrain 2025", lab)] = load_sbrain(ct)
    return out


# ----------------------------------------------------------------------------------------------
# build test
# ----------------------------------------------------------------------------------------------

def build_test():
    from liftover_py import Lift
    kg, _, _ = load_kg()
    kgkeys = set(kg.pkey)
    rows = []

    def test(name, pos, a, b, declared):
        """Share of the dataset's (position, allele pair) keys present in 1000G phase 3 GRCh37."""
        d = pd.DataFrame({"pos": pos, "a": a, "b": b}).dropna()
        d = d[(d.pos >= 16_609_414) & (d.pos <= 19_137_854)]
        k = [pkey(p, x, y) for p, x, y in zip(d.pos, d.a, d.b)]
        hit = np.mean([x in kgkeys for x in k]) if k else np.nan
        k1 = [pkey(p + 1, x, y) for p, x, y in zip(d.pos, d.a, d.b)]
        hit1 = np.mean([x in kgkeys for x in k1]) if k1 else np.nan
        rows.append({"veri": name, "beyan_edilen_yapi": declared, "karsilastirma": "1000G faz 3 GRCh37 (konum+allel)",
                     "varyant": len(k), "GRCh37_eslesme_orani": hit, "konum+1_eslesme_orani": hit1})

    def test_rs(name, rsid, pos, declared, against, against_name):
        m = pd.DataFrame({"rsid": rsid, "pos": pos}).dropna()
        m = m[m.rsid.isin(against.index)].drop_duplicates("rsid")
        d = (m.pos.astype(int).values - against.loc[m.rsid].values.astype(int))
        vals, cnt = np.unique(d, return_counts=True)
        top = sorted(zip(cnt, vals), reverse=True)[:3]
        rows.append({"veri": name, "beyan_edilen_yapi": declared, "karsilastirma": against_name,
                     "varyant": len(m), "GRCh37_eslesme_orani": float((d == 0).mean()) if len(d) else np.nan,
                     "konum+1_eslesme_orani": "en sik farklar " + "; ".join(f"{c}:{v}" for c, v in top)})

    for p, lab in (("1769", "BIG40 L3 ACR L"), ("1925", "BIG40 ICVF ACR L"), ("1475", "BIG40 FA ACR L")):
        b = load_big40(p)
        test(lab, b.pos, b.a_eff, b.a_oth, "GRCh37 (BIG40 sayfasi)")
    u = load_udipfa()
    test("GCST90860937 UDIP-FA", u.pos, u.a_eff, u.a_oth, "GRCh37 (GWAS Catalog meta.yaml)")
    s_ = load_sulc()
    test("GCST90095129 sulkus derinligi", s_.pos, s_.a_eff, s_.a_oth, "GRCh37 (dosya adi)")
    e = load_bryois("Oligodendrocytes").drop_duplicates("rsid")
    test("Bryois snp_pos SNP_id_hg19", e.pos37, e.a_eff, e.a_oth, "hg19 (yazarlarin liftOver'i)")
    test("Bryois snp_pos SNP_id_hg38", e.pos38, e.a_eff, e.a_oth, "hg38 (genotipler)")
    lf = Lift(f"{W}/raw/ucsc/hg38ToHg19.over.chain.gz", "6")
    z = load_zeng("Oligo")
    z = z[z.chr == 6].drop_duplicates("rsid")
    test("Zeng 2026 Oligo 'end' (ham)", z.end, z.a_eff, z.a_oth, "? (hg38 bekleniyor)")
    z["pos37_lift"] = [(lambda h: h[1] if h and h[0] == "chr6" else np.nan)(lf.map(int(p))) for p in z.end]
    test("Zeng 'end' -> python liftOver", z.pos37_lift, z.a_eff, z.a_oth, "hg19 (bizim liftOver)")
    b38 = e.set_index("rsid").pos38
    test_rs("Zeng 2026 Oligo 'start'", z.rsid, z.start, "?", b38, "Bryois hg38 (rsID)")
    test_rs("Zeng 2026 Oligo 'end'", z.rsid, z.end, "?", b38, "Bryois hg38 (rsID)")
    sb = load_sbrain("OD").drop_duplicates("rsid")
    test("SingleBrain OD 'pos' (ham)", sb.pos, sb.a_eff, sb.a_oth, "hg38 (Zenodo aciklamasi)")
    test_rs("SingleBrain OD 'pos'", sb.rsid, sb.pos, "hg38", b38, "Bryois hg38 (rsID)")
    e["pos37_lift"] = [(lambda h: h[1] if h and h[0] == "chr6" else np.nan)(lf.map(int(p))) for p in e.pos38]
    test_rs("Bryois hg38 -> python liftOver", e.rsid, e.pos37_lift, "hg19 (bizim liftOver)",
            e.set_index("rsid").pos37, "Bryois hg19 (yazarlarin liftOver'i)")
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/build_test.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 60):
        print(out.to_string())
    tags = kg[kg.rsid.isin(TAGS)][["rsid", "pos37", "ref", "alt", "af"]].copy()
    tags["rol"] = tags.rsid.map(TAG_ROLE)
    tags.to_csv(f"{RES}/tag_variants_panel.tsv", sep="\t", index=False)
    print(tags.to_string())


# ----------------------------------------------------------------------------------------------
# harmonised merge + coloc
# ----------------------------------------------------------------------------------------------

def gwas_labf(m, kind):
    if kind == "beta":
        return cp.labf_beta(m.beta.values, m.se.values ** 2, sdY=1.0, kind="quant")
    return cp.labf_p(m.p.values, m.maf_ref.values, m.n.values, sdY=1.0, kind="quant")


def merge(gw, eq):
    kg, _, _ = load_kg()
    gw = gw[(gw.pos >= WLO) & (gw.pos <= WHI)].copy()
    gw["key"] = [ukey(r, a, b) for r, a, b in zip(gw.rsid, gw.a_eff, gw.a_oth)]
    eq = eq.copy()
    eq["key"] = [ukey(r, a, b) for r, a, b in zip(eq.rsid, eq.a_eff, eq.a_oth)]
    m = gw.merge(eq[["key", "rsid", "beta", "se", "p", "a_eff", "maf"]].rename(
        columns={"rsid": "rsid_e", "beta": "beta_e", "se": "se_e", "p": "p_e", "a_eff": "a_eff_e", "maf": "maf_e"}),
        on="key", how="inner")
    m = m.merge(kg[["key", "af"]].drop_duplicates("key").rename(columns={"af": "af_kg"}), on="key", how="left")
    maf = np.minimum(m.eaf.values, 1 - m.eaf.values)
    fill = np.isnan(maf)
    maf[fill] = np.minimum(m.af_kg.values, 1 - m.af_kg.values)[fill]
    fill = np.isnan(maf)
    maf[fill] = m.maf_e.values[fill]
    m["maf_ref"] = maf
    m = m[np.isfinite(m.se_e) & (m.se_e > 0) & np.isfinite(m.p)]
    m = m[~(m.maf_ref.isna() & m.beta.isna())]
    m = m.drop_duplicates("key")
    return m.sort_values("pos").reset_index(drop=True)


def r2_between(rs1, rs2, pops=("EUR",)):
    kg, G, spop = load_kg()
    idx = {r: i for i, r in enumerate(kg.rsid)}
    if rs1 not in idx or rs2 not in idx:
        return np.nan
    cols = np.isin(spop, pops)
    a, b = G[idx[rs1], cols].astype(float), G[idx[rs2], cols].astype(float)
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1] ** 2)


def run_main():
    gsets = gwas_sets()
    rows = []
    for (src, ct), eq in eqtl_sets().items():
        for gene in GENES:
            e = eq[eq.symbol == gene]
            if len(e) < 10:
                rows.append({"eqtl_kaynak": src, "hucre_tipi": ct, "gen": gene, "gwas": "-", "nsnps": len(e),
                             "not": "gen bu hucre tipinde test edilmemis (ifade esigi alti)"})
                continue
            for gname, (gw, kind) in gsets.items():
                m = merge(gw, e)
                if len(m) < 10:
                    continue
                l1 = gwas_labf(m, kind)
                l2 = cp.labf_beta(m.beta_e.values, m.se_e.values ** 2, sdY=1.0, kind="quant")
                r = cp.coloc_abf(l1, l2)
                r5 = cp.coloc_abf(l1, l2, p12=5e-6)
                r1 = cp.coloc_abf(l1, l2, p12=1e-6)
                ie, ig = int(np.argmin(m.p_e.values)), int(np.argmin(m.p.values))
                itop = int(np.argmax(r["snp_pp_h4"]))
                rows.append({
                    "eqtl_kaynak": src, "hucre_tipi": ct, "gen": gene, "gwas": gname, "nsnps": r["nsnps"],
                    "PP.H0": r["PP.H0"], "PP.H1": r["PP.H1"], "PP.H2": r["PP.H2"], "PP.H3": r["PP.H3"],
                    "PP.H4": r["PP.H4"], "PP.H4_p12_5e-6": r5["PP.H4"], "PP.H4_p12_1e-6": r1["PP.H4"],
                    "eqtl_min_p": m.p_e.min(), "eqtl_oncu": m.rsid_e.iloc[ie], "eqtl_oncu_pos37": int(m.pos.iloc[ie]),
                    "gwas_min_p": m.p.min(), "gwas_oncu": m.rsid.iloc[ig], "gwas_oncu_pos37": int(m.pos.iloc[ig]),
                    "r2_oncu_EUR": r2_between(m.rsid_e.iloc[ie], m.rsid.iloc[ig]),
                    "eqtl_p_gwas_oncusunde": m.p_e.iloc[ig],
                    "H4_en_yuksek_snp": m.rsid_e.iloc[itop], "H4_snp_pp": float(r["snp_pp_h4"][itop]),
                    "H4>=0.5_icin_gereken_p12": cp.p12_for_h4(l1, l2, 0.5)})
            log(f"main {src} {ct} {gene} done")
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/coloc_main.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        k = out[(out.gen == "KIF13A") & (out.gwas != "-")]
        print(k[["eqtl_kaynak", "hucre_tipi", "gwas", "nsnps", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
                 "eqtl_min_p", "eqtl_oncu", "gwas_min_p", "gwas_oncu", "r2_oncu_EUR", "eqtl_p_gwas_oncusunde"]].to_string())


# ----------------------------------------------------------------------------------------------
# positive controls
# ----------------------------------------------------------------------------------------------

POSCTRL_BRYOIS = [  # Bryois 2022 Table S5 - AD (41593_2022_1128_MOESM3_ESM.xlsx), read 2026-10-03
    {"gene": "CR1", "locus": ("1", 207_302_548, 207_900_071), "nsnps": 955, "H3": 0.006599, "H4": 0.993401},
    {"gene": "NDUFS2", "locus": ("1", 161_014_370, 161_228_800), "nsnps": 285, "H3": 0.011350, "H4": 0.926255},
    {"gene": "ADAMTS4", "locus": ("1", 161_014_370, 161_228_800), "nsnps": 285, "H3": 0.019159, "H4": 0.767585},
]


def run_posctrl():
    rows = []
    for pc in POSCTRL_BRYOIS:
        c, lo, hi = pc["locus"]
        tag = "CR1" if pc["gene"] == "CR1" else "ADAMTS4"
        gw = pd.read_csv(f"{WORK}/gwas_ad_{tag}.tsv.gz", sep="\t")
        gw = gw[(gw.base_pair_location >= lo) & (gw.base_pair_location <= hi)].dropna(subset=["p_value"])
        gw["SNP_id_hg19"] = "chr" + gw.chromosome.astype(str) + ":" + gw.base_pair_location.astype(str)
        e = load_bryois("Oligodendrocytes", chrom="1")
        e = e[e.symbol == pc["gene"]]
        # Bryois 04_coloc.Rmd: eQTL rows whose hg19 id is in the GWAS locus, chr:pos match, no allele check
        m = gw.merge(e[["SNP_id_hg19", "beta", "se", "p"]].rename(columns={"beta": "beta_e", "se": "se_e", "p": "p_e"}),
                     on="SNP_id_hg19", how="inner")
        m = m[np.isfinite(m.se_e) & (m.se_e > 0)]
        l1 = cp.labf_beta(m.beta.values, m.standard_error.values ** 2, kind="cc")
        l2 = cp.labf_beta(m.beta_e.values, m.se_e.values ** 2, sdY=1.0, kind="quant")
        r = cp.coloc_abf(l1, l2)
        rows.append({"kaynak": "Bryois 2022 Table S5", "gwas": "AD Schwartzentruber 2021 (GCST90012877)",
                     "gen": pc["gene"], "hucre_tipi": "Oligodendrocytes", "p12": 1e-5,
                     "nsnps_bizim": r["nsnps"], "nsnps_yayin": pc["nsnps"],
                     "PP.H3_bizim": r["PP.H3"], "PP.H3_yayin": pc["H3"],
                     "PP.H4_bizim": r["PP.H4"], "PP.H4_yayin": pc["H4"]})
    # Zeng 2026: coloc_class.tsv, AD x Oligo. Their run_coloc.R: beta = se * z, varbeta = se^2,
    # GWAS type cc, p1 = p2 = 1e-4, p12 = 1e-6 (the preprint text says "default"); full cis window.
    zc = pd.read_csv(f"{W}/raw/synapse/coloc_class.tsv", sep="\t")
    if os.path.exists(f"{WORK}/zeng_Oligo.tsv.gz"):
        z = load_zeng("Oligo")
        for gene, tag in (("CR1", "CR1"), ("ADAMTS4", "ADAMTS4"), ("NDUFS2", "ADAMTS4")):
            pub = zc[(zc.Trait == "AD") & (zc.CellType == "Oligo") & (zc.Gene == gene)]
            if not os.path.exists(f"{WORK}/gwas_bell_{tag}.tsv.gz"):
                continue
            gw = pd.read_csv(f"{WORK}/gwas_bell_{tag}.tsv.gz", sep="\t")
            gw = gw.dropna(subset=["beta", "standard_error"])
            e = z[z.symbol == gene].copy()
            e["pos38"] = e.end.astype(int)
            for match in ("pos+allel", "rsid"):
                if match == "rsid":
                    m = gw.merge(e[["rsid", "a_eff", "a_oth", "beta", "se"]].rename(
                        columns={"rsid": "variant_id", "beta": "beta_e", "se": "se_e"}), on="variant_id")
                else:
                    gw["k"] = [f"{p}:{min(a,b)}:{max(a,b)}" for p, a, b in zip(gw.base_pair_location, gw.effect_allele, gw.other_allele)]
                    e["k"] = [f"{p}:{min(a,b)}:{max(a,b)}" for p, a, b in zip(e.pos38, e.a_eff, e.a_oth)]
                    m = gw.merge(e[["k", "beta", "se"]].rename(columns={"beta": "beta_e", "se": "se_e"}), on="k")
                m = m.drop_duplicates("variant_id")
                l1 = cp.labf_beta(m.beta.values, m.standard_error.values ** 2, kind="cc")
                l2 = cp.labf_beta(m.beta_e.values, m.se_e.values ** 2, sdY=1.0, kind="quant")
                for p12 in (1e-6, 1e-5):
                    r = cp.coloc_abf(l1, l2, p12=p12)
                    rows.append({"kaynak": f"Zeng 2026 coloc_class.tsv (eslesme: {match})",
                                 "gwas": "AD Bellenguez 2022 (GCST90027158)", "gen": gene, "hucre_tipi": "Oligo",
                                 "p12": p12, "nsnps_bizim": r["nsnps"], "nsnps_yayin": np.nan,
                                 "PP.H3_bizim": r["PP.H3"], "PP.H3_yayin": float(pub.ppH3.iloc[0]) if len(pub) else np.nan,
                                 "PP.H4_bizim": r["PP.H4"], "PP.H4_yayin": float(pub.ppH4.iloc[0]) if len(pub) else np.nan})
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/posctrl.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250):
        print(out.to_string())


# ----------------------------------------------------------------------------------------------
# targeted lookup at the haplotype tags + implied N
# ----------------------------------------------------------------------------------------------

def run_lookup():
    kg, _, _ = load_kg()
    af = kg.drop_duplicates("rsid").set_index("rsid").af
    rows, nrows = [], []
    for (src, ct), eq in eqtl_sets().items():
        for gene in GENES:
            e = eq[eq.symbol == gene]
            if not len(e):
                continue
            ie = int(np.argmin(e.p.values))
            for tag in TAGS:
                t = e[e.rsid == tag]
                if len(t):
                    r = t.iloc[0]
                    rows.append({"eqtl_kaynak": src, "hucre_tipi": ct, "gen": gene, "varyant": tag,
                                 "rol": TAG_ROLE[tag], "test_edildi": True, "etki_alleli": r.a_eff,
                                 "beta": r.beta, "se": r.se, "p": r.p,
                                 "gen_oncusu": e.rsid.iloc[ie], "gen_oncu_p": e.p.iloc[ie]})
                else:
                    rows.append({"eqtl_kaynak": src, "hucre_tipi": ct, "gen": gene, "varyant": tag,
                                 "rol": TAG_ROLE[tag], "test_edildi": False,
                                 "gen_oncusu": e.rsid.iloc[ie], "gen_oncu_p": e.p.iloc[ie]})
        e = eq.copy()
        f = e.maf.values if e.maf.notna().any() else e.rsid.map(af).values
        f = np.minimum(f, 1 - f)
        ok = (f > 0.05) & (e.se.values > 0) & np.isfinite(e.se.values)
        neff = 1.0 / (2 * f[ok] * (1 - f[ok]) * e.se.values[ok] ** 2)
        nrows.append({"eqtl_kaynak": src, "hucre_tipi": ct, "genler": ",".join(sorted(eq.symbol.unique())),
                      "N_etkin_ortanca": float(np.median(neff)) if ok.any() else np.nan,
                      "N_etkin_IQR": f"{np.percentile(neff,25):.0f}-{np.percentile(neff,75):.0f}" if ok.any() else "",
                      "maf_kaynagi": "veri" if eq.maf.notna().any() else "1000G EUR (yaklasik)"})
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/tag_lookup.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 300):
        print(out[out.gen == "KIF13A"].to_string())
    pd.DataFrame(nrows).to_csv(f"{RES}/n_eff.tsv", sep="\t", index=False)
    r2rows = []
    for a in TAGS:
        for b in TAGS:
            if a < b:
                r2rows.append({"v1": a, "v2": b, "r2_EUR": r2_between(a, b),
                               "r2_1000G_tumu": r2_between(a, b, pops=("EUR", "AFR", "AMR", "EAS", "SAS"))})
    pd.DataFrame(r2rows).to_csv(f"{RES}/tag_r2.tsv", sep="\t", index=False)
    print(pd.DataFrame(nrows).to_string())


# ----------------------------------------------------------------------------------------------
# conditional analysis
# ----------------------------------------------------------------------------------------------

def ld_matrix(G, idx, cols=None):
    X = G[idx][:, cols].astype(float) if cols is not None else G[idx].astype(float)
    X = X - X.mean(1, keepdims=True)
    sd = X.std(1, keepdims=True)
    sd[sd == 0] = 1
    X = X / sd
    return (X @ X.T) / X.shape[1]


def conditional_z(z, R, cond):
    if not cond:
        return z.copy()
    C = np.array(cond)
    Rcc = R[np.ix_(C, C)] + 1e-3 * np.eye(len(C))
    Rjc = R[:, C]
    A = Rjc @ np.linalg.inv(Rcc)
    num = z - A @ z[C]
    den = 1 - np.einsum("ij,ij->i", A, Rjc)
    den = np.where(den < 0.05, np.nan, den)
    return num / np.sqrt(den)


def oriented(m, kg_idx, kg):
    """eQTL and GWAS z oriented to the 1000G alt allele."""
    alt = kg.alt.values[kg_idx].astype(str)
    s_e = np.where(m.a_eff_e.astype(str).str.upper().values == np.char.upper(alt), 1.0, -1.0)
    z_e = s_e * m.beta_e.values / m.se_e.values
    if m.beta.notna().all():
        s_g = np.where(m.a_eff.astype(str).str.upper().values == np.char.upper(alt), 1.0, -1.0)
        z_g = s_g * m.beta.values / m.se.values
    else:
        z_g = None
    return z_e, z_g


def ld_sign_check(z, R, rsids, thr=0.9, zmin=1.5):
    """LD-consistency diagnostic: among variant pairs with panel r^2 >= thr and both |z| >= zmin,
    the share whose z signs disagree with the sign of r. Near 0 means the panel fits the data."""
    iu = np.triu_indices(len(z), 1)
    r = R[iu]
    keep = (r ** 2 >= thr) & (np.abs(z[iu[0]]) >= zmin) & (np.abs(z[iu[1]]) >= zmin)
    if not keep.any():
        return 0, np.nan
    agree = np.sign(z[iu[0]][keep]) * np.sign(z[iu[1]][keep]) == np.sign(r[keep])
    return int(keep.sum()), float(1 - agree.mean())


def run_cond():
    """Stepwise conditional eQTL signals with reference LD and coloc on each conditional signal.
    The eQTL ABF uses the dataset's own varbeta (se^2 / (1 - R_jC R_CC^-1 R_Cj)), i.e. the same
    scale as the main coloc.abf. A step is flagged as an LD-mismatch suspect when the selected
    lead's conditional |z| exceeds its marginal |z| by more than 3 (conditioning created it)."""
    kg, G, spop = load_kg()
    kidx = {k: i for i, k in enumerate(kg.key) if k}
    gsets = gwas_sets()
    ld_sets = {"EUR": np.where(spop == "EUR")[0], "1000G tumu": np.arange(len(spop))}
    rows, diag = [], []
    for (src, ct), eq in eqtl_sets().items():
        if ct == "pseudobulk":
            continue
        e = eq[eq.symbol == "KIF13A"]
        if len(e) < 10:
            continue
        for gname in ("BIG40 L3 sol on korona radiata (UKB 25319)", "UDIP-FA cok degiskenli (GCST90860937)"):
            gw, kind = gsets[gname]
            m = merge(gw, e)
            m = m[m.key.isin(kidx)].reset_index(drop=True)
            if len(m) < 50:
                continue
            idx = [kidx[k] for k in m.key]
            z_e, _ = oriented(m, idx, kg)
            se2 = m.se_e.values ** 2
            l1 = gwas_labf(m, kind)
            ig = int(np.argmin(m.p.values))
            for ldname, cols in ld_sets.items():
                if src != "Zeng 2026" and ldname != "EUR":
                    continue
                R = ld_matrix(G, idx, cols)
                npair, disc = ld_sign_check(z_e, R, m.rsid_e.values)
                if gname.startswith("BIG40"):
                    diag.append({"eqtl_kaynak": src, "hucre_tipi": ct, "LD_paneli": ldname,
                                 "r2>=0.9_cift": npair, "isaret_uyumsuz_oran": disc})
                conds, step = [], 0
                while step < 4:
                    zc = conditional_z(z_e, R, conds)
                    Rc = R[:, conds] if conds else None
                    if conds:
                        C = np.array(conds)
                        A = R[:, C] @ np.linalg.inv(R[np.ix_(C, C)] + 1e-3 * np.eye(len(C)))
                        den = np.clip(1 - np.einsum("ij,ij->i", A, R[:, C]), 0.05, 1)
                    else:
                        den = np.ones(len(z_e))
                    pcnd = 2 * stats.norm.sf(np.abs(np.nan_to_num(zc, nan=0.0)))
                    j = int(np.argmin(np.where(np.isnan(zc), np.inf, pcnd)))
                    V = se2 / den
                    rr = 0.0225 / (0.0225 + V)
                    l2 = 0.5 * np.log(1 - rr) + 0.5 * rr * np.nan_to_num(zc, nan=0.0) ** 2
                    r = cp.coloc_abf(l1, l2)
                    suspect = abs(np.nan_to_num(zc[j])) > abs(z_e[j]) + 3
                    rows.append({"eqtl_kaynak": src, "hucre_tipi": ct, "gwas": gname, "LD_paneli": ldname,
                                 "adim": step, "kosullanan": ",".join(m.rsid_e.iloc[conds]) if conds else "-",
                                 "sinyal_oncusu": m.rsid_e.iloc[j], "sinyal_p_marjinal": m.p_e.iloc[j],
                                 "sinyal_p_kosullu": pcnd[j], "LD_uyumsuzlugu_supheli": bool(suspect),
                                 "r2_gwas_oncusuyle": float(R[j, ig] ** 2), "gwas_oncu": m.rsid.iloc[ig],
                                 "eqtl_p_marjinal_gwas_oncusunde": m.p_e.iloc[ig],
                                 "eqtl_p_kosullu_gwas_oncusunde": pcnd[ig],
                                 "nsnps": len(m), "PP.H3": r["PP.H3"], "PP.H4": r["PP.H4"]})
                    if pcnd[j] > 1e-4 or suspect:
                        break
                    conds.append(j)
                    step += 1
            log(f"cond {src} {ct} {gname[:10]} done")
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/coloc_conditional.tsv", sep="	", index=False)
    pd.DataFrame(diag).to_csv(f"{RES}/ld_uyum_tanisi.tsv", sep="	", index=False)
    with pd.option_context("display.width", 300):
        print(out.to_string())
        print(pd.DataFrame(diag).to_string())


# ----------------------------------------------------------------------------------------------
# power simulation
# ----------------------------------------------------------------------------------------------

def run_power(nsim=1000, seed=7):
    """Power of the tests used here if the white-matter haplotype truly changed KIF13A expression in
    a cell type by b SD per allele (MetaBrain bulk cortex signal 3: b = 0.388).

    Simulated on the BIG40 L3 SNP set with 1000G EUR LD; the causal variant is the GWAS lead
    (rs144963738). Scenario 'yalin': the haplotype is the only eQTL. Scenario 'birincil eQTL var':
    the observed Zeng Oligo primary eQTL (rs113815188, z as observed) is added, which is the real
    situation; there coloc.abf is run both marginally and after conditioning on the primary lead.
    Also the analytic power of the single-variant lookup at the tag (two-sided alpha 0.05)."""
    rng = np.random.default_rng(seed)
    kg, G, spop = load_kg()
    kidx = {k: i for i, k in enumerate(kg.key) if k}
    gsets = gwas_sets()
    cols = np.where(spop == "EUR")[0]
    gname = "BIG40 L3 sol on korona radiata (UKB 25319)"
    gw, kind = gsets[gname]
    e = load_zeng("Oligo")
    e = e[e.symbol == "KIF13A"]
    m = merge(gw, e)
    m = m[m.key.isin(kidx)].reset_index(drop=True)
    idx = [kidx[k] for k in m.key]
    R = ld_matrix(G, idx, cols)
    w, U = np.linalg.eigh(0.999 * R + 0.001 * np.eye(len(R)))
    Lc = U * np.sqrt(np.clip(w, 0, None))
    z_obs, _ = oriented(m, idx, kg)
    l1 = gwas_labf(m, kind)
    c = int(np.argmin(m.p.values))
    j0 = int(np.argmin(m.p_e.values))
    fc = float(m.maf_ref.values[c])
    rows = []
    for N in (700, 1100):
        V = 1.0 / N  # standardized scale, as in Zeng 2026 (se ~ 1/sqrt(N))
        rr = 0.0225 / (0.0225 + V)
        for b in (0.1, 0.2, 0.3, 0.388):
            ncp = b * np.sqrt(2 * fc * (1 - fc) * N)
            for scen in ("yalin", "birincil eQTL var"):
                mu = R[:, c] * ncp
                if scen != "yalin":
                    mu = mu + R[:, j0] * z_obs[j0]
                h4m, h4c, hit = np.empty(nsim), np.empty(nsim), np.empty(nsim)
                for s in range(nsim):
                    z = mu + Lc @ rng.standard_normal(len(mu))
                    l2 = 0.5 * np.log(1 - rr) + 0.5 * rr * z ** 2
                    h4m[s] = cp.coloc_abf(l1, l2)["PP.H4"]
                    if scen != "yalin":
                        zc = np.nan_to_num(conditional_z(z, R, [j0]), nan=0.0)
                        l2c = 0.5 * np.log(1 - rr) + 0.5 * rr * zc ** 2
                        h4c[s] = cp.coloc_abf(l1, l2c)["PP.H4"]
                    else:
                        h4c[s] = np.nan
                    hit[s] = abs(z[c]) > 1.96
                rows.append({"senaryo": scen, "N_etkin": N, "beta_SD_allel_basina": b, "maf_nedensel": fc,
                             "beklenen_z_etiket": ncp, "guc_etiket_p<0.05": float(hit.mean()),
                             "guc_PP.H4>=0.5_marjinal": float((h4m >= 0.5).mean()),
                             "guc_PP.H4>=0.8_marjinal": float((h4m >= 0.8).mean()),
                             "guc_PP.H4>=0.5_kosullu": float(np.nanmean(h4c >= 0.5)) if scen != "yalin" else np.nan,
                             "guc_PP.H4>=0.8_kosullu": float(np.nanmean(h4c >= 0.8)) if scen != "yalin" else np.nan})
                log(f"power N={N} b={b} {scen}: {rows[-1]}")
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/power_sim.tsv", sep="	", index=False)
    print(out.to_string())
    # allele frequency of the tag by super-population (for the multi-ancestry Zeng cohort)
    fr = []
    for tag in ("rs144963738", "rs75428660", "rs115824472", "rs76091407"):
        ii = np.where(kg.rsid.values == tag)[0]
        if len(ii):
            g = G[ii[0]]
            fr.append({"varyant": tag, **{p: float(g[spop == p].mean() / 2) for p in ("EUR", "AFR", "AMR", "EAS", "SAS")}})
    pd.DataFrame(fr).to_csv(f"{RES}/tag_frekans_1000G.tsv", sep="	", index=False)
    print(pd.DataFrame(fr).to_string())


def run_sdy():
    """Sensitivity of the main coloc.abf to the phenotype SD. Zeng 2026 betas are standardised
    (se ~ 1/sqrt(N) whatever the MAF), so sdY = 1 holds there by construction. SingleBrain and Bryois
    report per-allele betas on residualised expression whose SD is not 1: estimate it with coloc's
    sdY.est (regression of 2 N f (1-f) on 1/varbeta) at the nominal donor count, then rerun."""
    kg, _, _ = load_kg()
    af = kg.dropna(subset=["rsid"]).drop_duplicates("rsid").set_index("rsid").af
    gsets = gwas_sets()
    rows = []
    nominal = {"SingleBrain 2025": 983, "Bryois 2022": 192}
    for (src, ct), eq in eqtl_sets().items():
        if src not in nominal:
            continue
        e = eq[eq.symbol == "KIF13A"].copy()
        if len(e) < 10:
            continue
        f = e.maf.values if e.maf.notna().any() else e.rsid.map(af).values
        f = np.minimum(f, 1 - f)
        ok = np.isfinite(f) & (f > 0.01) & (e.se.values > 0)
        nvx = 2 * nominal[src] * f[ok] * (1 - f[ok])
        oneover = 1 / e.se.values[ok] ** 2
        sdy = float(np.sqrt(np.sum(nvx * oneover) / np.sum(oneover ** 2)))
        for gname in ("BIG40 L3 sol on korona radiata (UKB 25319)", "UDIP-FA cok degiskenli (GCST90860937)"):
            gw, kind = gsets[gname]
            m = merge(gw, e)
            l1 = gwas_labf(m, kind)
            r1 = cp.coloc_abf(l1, cp.labf_beta(m.beta_e.values, m.se_e.values ** 2, sdY=1.0))
            rs = cp.coloc_abf(l1, cp.labf_beta(m.beta_e.values, m.se_e.values ** 2, sdY=sdy))
            rows.append({"eqtl_kaynak": src, "hucre_tipi": ct, "gwas": gname, "N_nominal": nominal[src],
                         "sdY_tahmini": sdy, "PP.H3_sdY1": r1["PP.H3"], "PP.H4_sdY1": r1["PP.H4"],
                         "PP.H3_sdYtahmin": rs["PP.H3"], "PP.H4_sdYtahmin": rs["PP.H4"]})
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/coloc_sdY_duyarlilik.tsv", sep="	", index=False)
    with pd.option_context("display.width", 250):
        print(out.to_string())


# ----------------------------------------------------------------------------------------------
# export for the R validation
# ----------------------------------------------------------------------------------------------

def run_export():
    kg, G, spop = load_kg()
    kidx = {k: i for i, k in enumerate(kg.key)}
    gsets = gwas_sets()
    os.makedirs(f"{RES}/r_input", exist_ok=True)
    cols = np.where(spop == "EUR")[0]
    es = eqtl_sets()
    picks = [("Zeng 2026", "Oligo"), ("Zeng 2026", "OPC"), ("Zeng 2026", "EN"), ("Zeng 2026", "Astro"),
             ("Bryois 2022", "Oligo"), ("Bryois 2022", "OPC"), ("SingleBrain 2025", "Oligo"),
             ("SingleBrain 2025", "OPC")]
    for src, ct in picks:
        if (src, ct) not in es:
            continue
        e = es[(src, ct)]
        e = e[e.symbol == "KIF13A"]
        for gname, tag in (("UDIP-FA cok degiskenli (GCST90860937)", "udipfa"),
                           ("BIG40 L3 sol on korona radiata (UKB 25319)", "big40L3")):
            gw, kind = gsets[gname]
            m = merge(gw, e)
            m = m[m.key.isin(kidx)].reset_index(drop=True)
            idx = [kidx[k] for k in m.key]
            z_e, z_g = oriented(m, idx, kg)
            m["z_e_alt"] = z_e
            m["z_g_alt"] = z_g if z_g is not None else np.nan
            m["snp"] = m.key
            fn = f"{tag}__{src.split()[0]}_{ct}"
            m[["snp", "rsid", "pos", "beta", "se", "p", "n", "maf_ref", "beta_e", "se_e", "p_e", "z_e_alt",
               "z_g_alt"]].to_csv(f"{RES}/r_input/{fn}.tsv", sep="\t", index=False)
            np.savetxt(f"{RES}/r_input/{fn}.ld", ld_matrix(G, idx, cols), fmt="%.5f")
    log("export done")


if __name__ == "__main__":
    for st in sys.argv[1:]:
        {"build": build_test, "posctrl": run_posctrl, "main": run_main, "lookup": run_lookup,
         "cond": run_cond, "power": run_power, "export": run_export, "sdy": run_sdy}[st]()
