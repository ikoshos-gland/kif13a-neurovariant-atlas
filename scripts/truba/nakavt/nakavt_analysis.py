"""C ekseni (nakavt): KIF13A icin cekinik (iki kopya kaybi) sorunun hesaplari.

TRUBA compute node'da SLURM isi olarak calisir (login node'da asla). Girdiler fetch_sources.py ve
inspect_sources.py'nin urettikleri; ciktilar interim/nakavt/ altina TSV + JSON.

Bolumler
  1. gnomAD v4.1: populasyon bazli LoF frekanslari, NMD kurali, varyant siniflari
  2. Genes & Health (Kim 2026, ST16): kohort ici kalibrasyon ve KIF13A beklenen/gozlenen
  3. RGC-ME (Sun 2024, ST2 + ST6): ayni hesap, 983.578 kisi
  4. gnomAD v4.1 homozigotlari: sinonim varyantlarla notr kalibrasyon, sonra gercek LoF
  5. Chundru 2024: gelisimsel bozukluk triolarinda gen duzeyi cekinik sonuc, ARDD zenginlesme carpani
  6. popEVE: KIF13A ve KIF1A skorlari, katman A ve KIF1A patojenik seti
  7. Ozet tablo
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

BASE = Path("/arf/scratch/avural/kif13a_nakavt")
RAW = BASE / "raw"
SH = BASE / "interim" / "sheets"
PRJ = BASE / "interim" / "project"
OUT = BASE / "interim" / "nakavt"
OUT.mkdir(parents=True, exist_ok=True)
TAG = "20261003"
S: dict = {}  # makine okunur ozet

GENES = {"KIF13A": "ENST00000259711", "KIF1A": "ENST00000498729", "KIF13B": "ENST00000524189"}
LOF_CSQ = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"}
SNV_CSQ = {"stop_gained", "splice_donor_variant", "splice_acceptor_variant"}
N_GH = 44028
N_RGC = 983578
# RGC-ME Ek Tablo 1A, iliskili analiz seti (>%50 kita kokeni)
RGC_ANC = {"AFR": 55281, "EAS": 5849, "EUR": 722477, "IAM": 85125, "MEA": 2229, "SAS": 30752, "UNK": 81865}
ANC2POP = {"AFR": "afr", "EAS": "eas", "EUR": "nfe", "IAM": "amr", "MEA": "mid", "SAS": "sas", "UNK": None}
# Varsayilan F (otozigotluk) degerleri. Kaynak: Chundru 2024 Tablo 1 (DD probandlarinin F_ROH'u,
# genel populasyondan biraz yuksek olmasi beklenir) ve G&H ana metin (genomun %1,8'i ROH'ta).
F_POP = {"afr": 0.0003, "amr": 0.0020, "asj": 0.0010, "eas": 0.0008, "fin": 0.0010, "mid": 0.0350,
         "nfe": 0.0005, "sas": 0.0125, "ami": 0.0200, "remaining": 0.0010}
POPS = list(F_POP)
F_GH_PAPER = 0.018
POPEVE_SEVERE = -5.056


def jdump(name: str, obj) -> None:
    (OUT / name).write_text(json.dumps(obj, indent=2, default=float), encoding="utf-8")


def read_rows(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def table_from(path: Path, first_col: str) -> pd.DataFrame:
    """Ek tablo sayfalarinda baslik satiri bastan birkac satir asagida olabilir."""
    rows = read_rows(path)
    h = next(i for i, r in enumerate(rows) if r and r[0].strip() == first_col)
    head = rows[h]
    body = [r + [""] * (len(head) - len(r)) for r in rows[h + 1:] if any(c.strip() for c in r)]
    return pd.DataFrame([r[:len(head)] for r in body], columns=head)


def num(s):
    return pd.to_numeric(s, errors="coerce")


def ko_model(q, F, n):
    """Beklenen iki kopya kayipli birey: N[(1-F)q^2 + Fq]."""
    q = np.asarray(q, float)
    return n * ((1 - F) * q ** 2 + F * q)


def fit_F(q, ko, n) -> dict:
    q = np.asarray(q, float)
    ko = np.asarray(ko, float)

    def nll(F):
        e = ko_model(q, F, n)
        return -(ko * np.log(e) - e).sum()

    r = optimize.minimize_scalar(nll, bounds=(1e-7, 0.3), method="bounded")
    F = float(r.x)
    e = ko_model(q, F, n)
    return {"F_eff": F, "n_genes": int(len(q)), "obs_total": float(ko.sum()), "exp_total": float(e.sum()),
            "obs_exp_F0": float(ko.sum() / ko_model(q, 0.0, n).sum())}


def pois_le(k, mu):
    return float(stats.poisson.cdf(k, mu))


def pois_ge(k, mu):
    return float(stats.poisson.sf(k - 1, mu)) if k > 0 else 1.0


# ----------------------------------------------------------------------------------- 0. referanslar


def load_constraint() -> pd.DataFrame:
    c = pd.read_csv(RAW / "gnomad_v41_constraint__gnomad.v4.1.constraint_metrics.tsv", sep="\t",
                    usecols=["gene", "gene_id", "transcript", "canonical", "mane_select", "lof.oe",
                             "lof.oe_ci.upper", "lof.pLI", "lof.obs", "lof.exp", "cds_length"],
                    low_memory=False)
    c = c[c["mane_select"].astype(str).str.lower() == "true"].copy()
    c = c.drop_duplicates("gene")
    return c.rename(columns={"lof.oe_ci.upper": "loeuf"})


def load_ddg2p() -> pd.DataFrame:
    with gzip.open(RAW / "ddg2p__DDG2P_2026-09-28.csv.gz", "rt", encoding="utf-8") as fh:
        d = pd.read_csv(fh)
    d["gene"] = d["gene symbol"]
    d["biallelic"] = d["allelic requirement"].fillna("").str.contains("biallelic_autosomal")
    d["monoallelic"] = d["allelic requirement"].fillna("").str.contains("monoallelic_autosomal")
    d["confident"] = d["confidence"].fillna("").isin(["definitive", "strong"])
    d["lof_mech"] = (d["variant consequence"].fillna("").str.contains("absent gene product")
                     | (d["molecular mechanism"].fillna("") == "loss of function"))
    d["dd_panel"] = d["panel"].fillna("").str.contains("DD")
    d["ardd_lof"] = d["biallelic"] & d["confident"] & d["lof_mech"] & d["dd_panel"]
    return d


def ddg2p_flags(d: pd.DataFrame) -> pd.DataFrame:
    g = d.groupby("gene").agg(ddg2p_any=("gene", "size"), ardd_lof=("ardd_lof", "any"),
                              mono_conf=("monoallelic", "any"))
    g["ddg2p_any"] = True
    return g.reset_index()


# ----------------------------------------------------------------------------------- 1. gnomAD


CDS_POS = re.compile(r"^c\.(?:\*|-)?(\d+)")


def structure(enst: str) -> dict:
    raw = json.loads((RAW / "ensembl" / f"{enst}_structure.json").read_text())
    tr = raw["Translation"]
    lo_c, hi_c = min(tr["start"], tr["end"]), max(tr["start"], tr["end"])
    exons, cum = [], 0
    for i, ex in enumerate(raw["Exon"], start=1):
        lo, hi = min(ex["start"], ex["end"]), max(ex["start"], ex["end"])
        coding = max(0, min(hi, hi_c) - max(lo, lo_c) + 1)
        exons.append({"n": i, "coding": coding, "cds_start": cum + 1 if coding else None,
                      "cds_end": cum + coding if coding else None})
        cum += coding
    last_junction = sum(e["coding"] for e in exons[:-1])
    return {"n_exons": len(exons), "cds_length": cum, "protein_length": tr["length"], "exons": exons,
            "last_junction_cds": last_junction, "escape_from_cds": last_junction - 50 + 1}


def exon_of(cds: int | None, st: dict):
    if cds is None:
        return None
    for e in st["exons"]:
        if e["cds_start"] and e["cds_start"] <= cds <= e["cds_end"]:
            return e["n"]
    return None


def gnomad_table(sym: str, kind: str = "pop") -> pd.DataFrame:
    enst = GENES[sym]
    fn = f"{sym}_{enst}_gnomad_r4_{'pop' if kind == 'pop' else 'allvariants_pop'}.json"
    payload = json.loads((RAW / "gnomad" / fn).read_text())
    st = structure(enst)
    rows = []
    for v in payload["variants"]:
        j = v.get("joint") or {}
        filters = j.get("filters") or []
        pops = {p["id"]: p for p in (j.get("populations") or []) if p.get("id") in F_POP}
        hgvsc = (v.get("hgvsc") or "").split(":")[-1]
        m = CDS_POS.match(hgvsc)
        cds = int(m.group(1)) if m else None
        mp = re.search(r"p\.[A-Za-z]{3}(\d+)", v.get("hgvsp") or "")
        r = {"gene": sym, "variant_id": v["variant_id"], "consequence": v.get("consequence"),
             "hgvsc": hgvsc, "hgvsp": (v.get("hgvsp") or "").split(":")[-1],
             "ref": v.get("ref") or v["variant_id"].split("-")[2],
             "alt": v.get("alt") or v["variant_id"].split("-")[3],
             "loftee": v.get("lof") or "", "lof_filter": v.get("lof_filter") or "",
             "cds_pos": cds, "exon": exon_of(cds, st), "n_exons": st["n_exons"],
             "protein_pos": int(mp.group(1)) if mp else None,
             "joint_ac": j.get("ac") or 0, "joint_an": j.get("an") or 0, "joint_hom": j.get("homozygote_count") or 0,
             "joint_filters": ";".join(filters)}
        r["joint_af"] = r["joint_ac"] / r["joint_an"] if r["joint_an"] else 0.0
        for p in POPS:
            pp = pops.get(p) or {}
            ac, an = pp.get("ac") or 0, pp.get("an") or 0
            r[f"ac_{p}"], r[f"an_{p}"], r[f"hom_{p}"] = ac, an, pp.get("homozygote_count") or 0
            r[f"af_{p}"] = ac / an if an else 0.0
        rows.append(r)
    df = pd.DataFrame(rows)
    df["confident"] = (df["joint_ac"] >= 1) & ~df["joint_filters"].str.contains("AC0")
    df["csq_lof"] = df["consequence"].isin(LOF_CSQ)
    df["nmd_escape"] = df["cds_pos"].apply(lambda c: bool(c is not None and not pd.isna(c) and c >= st["escape_from_cds"]))
    df["is_snv"] = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    df["cls_all_plof"] = df["confident"] & df["csq_lof"]
    df["cls_hc_any"] = df["cls_all_plof"] & (df["loftee"] == "HC")
    df["cls_hc_snv"] = df["cls_hc_any"] & df["consequence"].isin(SNV_CSQ) & df["is_snv"]
    df["cls_true_lof"] = df["cls_hc_any"] & (df["lof_filter"] == "") & ~df["nmd_escape"]
    df["cls_true_lof_snv"] = df["cls_true_lof"] & df["consequence"].isin(SNV_CSQ) & df["is_snv"]
    rgc_af = 0.0
    for anc, n in RGC_ANC.items():
        p = ANC2POP[anc]
        rgc_af = rgc_af + n * (df[f"af_{p}"] if p else df["joint_af"])
    df["af_rgc_weighted"] = rgc_af / sum(RGC_ANC.values())
    df.attrs["structure"] = st
    return df


CLASSES = ["cls_all_plof", "cls_hc_any", "cls_hc_snv", "cls_true_lof", "cls_true_lof_snv"]


def q_by_pop(df: pd.DataFrame) -> dict:
    out = {}
    for c in CLASSES:
        sub = df[df[c]]
        d = {"n_variants": int(len(sub)), "allele_count": int(sub["joint_ac"].sum()),
             "q_joint": float(sub["joint_af"].sum()), "q_rgc_weighted": float(sub["af_rgc_weighted"].sum())}
        for p in POPS:
            d[f"q_{p}"] = float(sub[f"af_{p}"].sum())
            d[f"hom_{p}"] = int(sub[f"hom_{p}"].sum())
        d["hom_joint"] = int(sub["joint_hom"].sum())
        out[c] = d
    return out


# ----------------------------------------------------------------------------------- main


def main() -> None:
    cons = load_constraint()
    dd = load_ddg2p()
    ddf = ddg2p_flags(dd)
    ann = cons[["gene", "loeuf", "lof.oe", "lof.pLI", "cds_length"]].merge(ddf, on="gene", how="left")
    ann[["ddg2p_any", "ardd_lof", "mono_conf"]] = ann[["ddg2p_any", "ardd_lof", "mono_conf"]].astype("boolean").fillna(False).astype(bool)

    # DDG2P'de hedef genler
    kd = dd[dd["gene"].isin(GENES)][["g2p id", "gene", "disease name", "allelic requirement", "confidence",
                                      "variant consequence", "molecular mechanism", "panel"]]
    kd.to_csv(OUT / f"ddg2p_target_genes_{TAG}.tsv", sep="\t", index=False)
    S["ddg2p_target_rows"] = kd.to_dict("records")
    S["ddg2p_n_ardd_lof_genes"] = int(dd.loc[dd["ardd_lof"], "gene"].nunique())
    print("DDG2P hedef satirlari:\n", kd.to_string(), flush=True)

    # ---------------------------------------------------------------- 1. gnomAD
    gq, gtabs = {}, {}
    for sym in GENES:
        t = gnomad_table(sym, "pop")
        gtabs[sym] = t
        gq[sym] = q_by_pop(t)
        st = t.attrs["structure"]
        gq[sym]["structure"] = {k: st[k] for k in ("n_exons", "cds_length", "protein_length",
                                                   "last_junction_cds", "escape_from_cds")}
        print(f"{sym}: rows {len(t)}, confident {int(t['confident'].sum())}, "
              f"all_plof {int(t['cls_all_plof'].sum())}, hc_any {int(t['cls_hc_any'].sum())}, "
              f"true_lof {int(t['cls_true_lof'].sum())}, escape_from_cds {st['escape_from_cds']}", flush=True)
    S["gnomad_q"] = gq
    k = gtabs["KIF13A"]
    keep = ["variant_id", "consequence", "hgvsc", "hgvsp", "exon", "cds_pos", "nmd_escape", "loftee", "lof_filter",
            "joint_ac", "joint_af", "joint_hom", "af_rgc_weighted"] + [f"af_{p}" for p in POPS] + \
           [f"hom_{p}" for p in POPS] + CLASSES
    kv = k[k["confident"] & k["csq_lof"]][keep].sort_values("joint_ac", ascending=False)
    kv.to_csv(OUT / f"KIF13A_gnomad_lof_by_pop_{TAG}.tsv", sep="\t", index=False)
    for sym in ("KIF1A", "KIF13B"):
        t = gtabs[sym]
        t[t["confident"] & t["csq_lof"]][keep].sort_values("joint_ac", ascending=False).to_csv(
            OUT / f"{sym}_gnomad_lof_by_pop_{TAG}.tsv", sep="\t", index=False)

    # Ekzon 38 (son intronun verici bolgesi c.4581+1) cerceve hesabi
    st = gtabs["KIF13A"].attrs["structure"]
    ex = {e["n"]: e for e in st["exons"]}
    S["KIF13A_exon38"] = {"exon38_coding_nt": ex[38]["coding"], "exon38_cds": [ex[38]["cds_start"], ex[38]["cds_end"]],
                         "exon39_cds": [ex[39]["cds_start"], ex[39]["cds_end"]],
                         "exon38_skip_in_frame": ex[38]["coding"] % 3 == 0}

    # ---------------------------------------------------------------- 2. Genes & Health
    gh = table_from(SH / "gh_supp_tables__41588_2026_2553_MOESM4_ESM" / "ST16.tsv", "Symbol")
    for c in ["Number of heterozygous carriers", "Number of homozygous carriers", "Number of compound heterozygous carriers"]:
        gh[c] = num(gh[c]).fillna(0)
    gh = gh.rename(columns={"Symbol": "gene", "Number of heterozygous carriers": "het",
                            "Number of homozygous carriers": "hom", "Number of compound heterozygous carriers": "chet"})
    gh["ko"] = gh["hom"] + gh["chet"]
    gh["q"] = (gh["het"] + 2 * gh["hom"] + 2 * gh["chet"]) / (2 * N_GH)
    gh = gh.merge(ann, on="gene", how="left")
    for c in ["ddg2p_any", "ardd_lof", "mono_conf"]:
        gh[c] = gh[c].astype("boolean").fillna(False).astype(bool)
    S["gh_n_genes"] = int(len(gh))
    S["gh_n_ko_genes"] = int((gh["ko"] > 0).sum())
    tol = gh[(gh["loeuf"] > 0.8) & ~gh["ddg2p_any"] & (gh["q"] > 0)]
    band = tol[(tol["q"] >= 2e-5) & (tol["q"] <= 2e-3)]
    fit_all = fit_F(band["q"], band["ko"], N_GH)
    near = tol[(tol["het"] >= 10) & (tol["het"] <= 50)]
    fit_near = fit_F(near["q"], near["ko"], N_GH)
    S["gh_fit_tolerant_band"] = fit_all
    S["gh_fit_tolerant_het10_50"] = fit_near
    print("G&H F fit (toleransli, q 2e-5..2e-3):", fit_all, flush=True)
    print("G&H F fit (toleransli, het 10..50):", fit_near, flush=True)

    def gh_group(df, lo, hi):
        g = df[(df["het"] >= lo) & (df["het"] <= hi)]
        return {"n_genes": int(len(g)), "frac_ko_ge1": float((g["ko"] > 0).mean()) if len(g) else None,
                "mean_ko": float(g["ko"].mean()) if len(g) else None,
                "exp_mean_ko_Feff": float(ko_model(g["q"], fit_near["F_eff"], N_GH).mean()) if len(g) else None}

    ardd = gh[gh["ardd_lof"] & (gh["q"] > 0)]
    S["gh_empirical_het15_30"] = {"tolerant": gh_group(tol, 15, 30), "ardd_lof": gh_group(ardd, 15, 30),
                                  "all_genes": gh_group(gh[gh["q"] > 0], 15, 30)}
    S["gh_empirical_het3_10"] = {"tolerant": gh_group(tol, 3, 10), "ardd_lof": gh_group(ardd, 3, 10)}
    a_band = ardd[(ardd["q"] >= 2e-5) & (ardd["q"] <= 2e-3)]
    e_a = ko_model(a_band["q"], fit_all["F_eff"], N_GH).sum()
    S["gh_ardd_depletion"] = {"n_genes": int(len(a_band)), "obs": float(a_band["ko"].sum()), "exp": float(e_a),
                              "obs_over_exp": float(a_band["ko"].sum() / e_a) if e_a else None,
                              "p_le_obs": pois_le(a_band["ko"].sum(), e_a)}
    gh_rows = {}
    for sym in GENES:
        r = gh[gh["gene"] == sym].iloc[0]
        gh_rows[sym] = {"het": float(r["het"]), "hom": float(r["hom"]), "chet": float(r["chet"]), "q": float(r["q"]),
                        "loeuf": float(r["loeuf"]),
                        "exp_F_eff_band": float(ko_model(r["q"], fit_all["F_eff"], N_GH)),
                        "exp_F_eff_het10_50": float(ko_model(r["q"], fit_near["F_eff"], N_GH)),
                        "exp_F_paper_0.018": float(ko_model(r["q"], F_GH_PAPER, N_GH)),
                        "exp_HWE_F0": float(ko_model(r["q"], 0.0, N_GH)),
                        "in_RGC_ME": r.get("In RGC-ME"), "in_gnomAD_v4": r.get("In gnomAD v4")}
    # gercek LoF payi: gnomAD SAS'ta NMD'ye duyarli HC / (G&H tanimina denk sinif)
    for sym in GENES:
        qs = gq[sym]
        for denom in ("cls_hc_any", "cls_all_plof"):
            pi = qs["cls_true_lof"]["q_sas"] / qs[denom]["q_sas"] if qs[denom]["q_sas"] else float("nan")
            qt = gh_rows[sym]["q"] * pi
            gh_rows[sym][f"true_lof_share_sas_vs_{denom}"] = pi
            gh_rows[sym][f"exp_true_lof_Feff_het10_50_vs_{denom}"] = float(ko_model(qt, fit_near["F_eff"], N_GH))
    S["gh_target_genes"] = gh_rows
    gh[gh["gene"].isin(GENES)].to_csv(OUT / f"gh_ST16_target_genes_{TAG}.tsv", sep="\t", index=False)
    print("G&H hedef genler:", json.dumps(gh_rows, indent=1, default=float), flush=True)

    # ---------------------------------------------------------------- 3. RGC-ME
    st2 = table_from(SH / "rgc_st2__41586_2024_7556_MOESM4_ESM" / "table_s2.tsv", "GeneName")
    st6 = table_from(SH / "rgc_st6__41586_2024_7556_MOESM8_ESM" / "suppl_table_6_HKOs.tsv", "CHR")
    for c in ["n", "N_total", "mean", "MAF"]:
        st2[c] = num(st2[c])
    for c in ["allhetRA", "sum_homAA", "homAlt_aaf", "comp_het", "compHet_aaf", "sum_KO"]:
        st6[c] = num(st6[c])
    S["rgc_st6_n_genes"] = int(len(st6))
    rg = st2.rename(columns={"GeneName": "gene"}).merge(
        st6.rename(columns={"GeneName": "gene", "TranscriptId": "TranscriptId_st6"})[
            ["gene", "TranscriptId_st6", "allhetRA", "sum_homAA", "homAlt_aaf", "comp_het", "sum_KO"]],
        on="gene", how="left")
    rg["ko"] = rg["sum_KO"].fillna(0)
    rg["q"] = rg["n"] / rg["N_total"]
    rg = rg.merge(ann.drop(columns=["loeuf"]).merge(cons[["gene", "loeuf"]], on="gene", how="left"),
                  on="gene", how="left", suffixes=("", "_v41"))
    for c in ["ddg2p_any", "ardd_lof", "mono_conf"]:
        rg[c] = rg[c].astype("boolean").fillna(False).astype(bool)
    tol_r = rg[(rg["loeuf_v41"] > 0.8) & ~rg["ddg2p_any"] & (rg["q"] > 0)]
    band_r = tol_r[(tol_r["q"] >= 2e-5) & (tol_r["q"] < 1e-3)]
    fit_r = fit_F(band_r["q"], band_r["ko"], N_RGC)
    near_r = tol_r[(tol_r["n"] >= 230) & (tol_r["n"] <= 920)]
    fit_rn = fit_F(near_r["q"], near_r["ko"], N_RGC)
    S["rgc_fit_tolerant_band"] = fit_r
    S["rgc_fit_tolerant_n230_920"] = fit_rn
    print("RGC F fit:", fit_r, fit_rn, flush=True)

    def rg_group(df, lo, hi):
        g = df[(df["n"] >= lo) & (df["n"] <= hi)]
        return {"n_genes": int(len(g)), "frac_ko_ge1": float((g["ko"] > 0).mean()) if len(g) else None,
                "mean_ko": float(g["ko"].mean()) if len(g) else None,
                "exp_mean_ko_Feff": float(ko_model(g["q"], fit_rn["F_eff"], N_RGC).mean()) if len(g) else None}

    ardd_r = rg[rg["ardd_lof"] & (rg["q"] > 0)]
    S["rgc_empirical_n300_650"] = {"tolerant": rg_group(tol_r, 300, 650), "ardd_lof": rg_group(ardd_r, 300, 650),
                                   "all_genes": rg_group(rg[rg["q"] > 0], 300, 650)}
    S["rgc_empirical_n10_40"] = {"tolerant": rg_group(tol_r, 10, 40), "ardd_lof": rg_group(ardd_r, 10, 40)}
    a_r = ardd_r[(ardd_r["q"] >= 2e-5) & (ardd_r["q"] < 1e-3)]
    e_ar = ko_model(a_r["q"], fit_r["F_eff"], N_RGC).sum()
    S["rgc_ardd_depletion"] = {"n_genes": int(len(a_r)), "obs": float(a_r["ko"].sum()), "exp": float(e_ar),
                               "obs_over_exp": float(a_r["ko"].sum() / e_ar) if e_ar else None,
                               "p_le_obs": pois_le(a_r["ko"].sum(), e_ar)}
    rgc_rows = {}
    for sym in GENES:
        r = rg[rg["gene"] == sym].iloc[0]
        e1 = float(ko_model(r["q"], fit_r["F_eff"], N_RGC))
        e2 = float(ko_model(r["q"], fit_rn["F_eff"], N_RGC))
        # gercek LoF: ayni q tanimi (HC SNV), gnomAD'den RGC kokenine agirliklandirilmis, RGC olcegine tasinmis
        qs = gq[sym]
        kscale = r["q"] / qs["cls_hc_snv"]["q_rgc_weighted"] if qs["cls_hc_snv"]["q_rgc_weighted"] else float("nan")
        q_true_snv = qs["cls_true_lof_snv"]["q_rgc_weighted"] * kscale
        rgc_rows[sym] = {"n_hc_snv_alleles": float(r["n"]), "N_total": float(r["N_total"]), "q": float(r["q"]),
                         "shet_mean": float(r["mean"]), "ko_obs": float(r["ko"]),
                         "sum_homAA": None if pd.isna(r["sum_homAA"]) else float(r["sum_homAA"]),
                         "comp_het": None if pd.isna(r["comp_het"]) else float(r["comp_het"]),
                         "homAlt_aaf": None if pd.isna(r["homAlt_aaf"]) else float(r["homAlt_aaf"]),
                         "allhetRA": None if pd.isna(r["allhetRA"]) else float(r["allhetRA"]),
                         "transcript_st2": r["TranscriptId"], "transcript_st6": r["TranscriptId_st6"],
                         "exp_Feff_band": e1, "exp_Feff_n230_920": e2,
                         "p_ge_obs_Feff_n230_920": pois_ge(int(r["ko"]), e2),
                         "rgc_over_gnomad_hcsnv_scale": float(kscale),
                         "q_true_lof_snv_rgc_scale": float(q_true_snv),
                         "exp_true_lof_Feff_n230_920": float(ko_model(q_true_snv, fit_rn["F_eff"], N_RGC))}
    S["rgc_target_genes"] = rgc_rows
    rg[rg["gene"].isin(GENES)].to_csv(OUT / f"rgc_ST2_ST6_target_genes_{TAG}.tsv", sep="\t", index=False)
    print("RGC hedef genler:", json.dumps(rgc_rows, indent=1, default=float), flush=True)

    # RGC'deki 2 homozigotun hangi varyant olabilecegi: kumulatif homozigot AAF ile gnomAD frekanslari
    target = rgc_rows["KIF13A"]["homAlt_aaf"]
    kv2 = kv.copy()
    kv2["ratio_to_rgc_homAlt_aaf"] = kv2["af_rgc_weighted"] / target
    # Kokene gore beklenen homozigot (RGC iliskili seti), varsayilan F'lerle
    def exp_hom_rgc(row, fmult=1.0):
        e = 0.0
        for anc, n in RGC_ANC.items():
            p = ANC2POP[anc]
            af = row[f"af_{p}"] if p else row["joint_af"]
            F = (F_POP[p] if p else 0.001) * fmult
            if anc == "SAS":
                F = 0.030 * fmult  # 19.509 kisilik akraba evliligi yuksek Pakistan kohortu SAS'in cogunu olusturur
            e += n * ((1 - F) * af ** 2 + F * af)
        return e
    kv2["exp_hom_rgc_Fdefault"] = kv2.apply(exp_hom_rgc, axis=1)
    cand = kv2[(kv2["ratio_to_rgc_homAlt_aaf"] > 0.33) & (kv2["ratio_to_rgc_homAlt_aaf"] < 3)]
    cand.to_csv(OUT / f"KIF13A_rgc_hom_candidates_{TAG}.tsv", sep="\t", index=False)
    S["rgc_KIF13A_hom_candidates"] = cand[["variant_id", "hgvsc", "hgvsp", "consequence", "loftee", "lof_filter",
                                           "nmd_escape", "joint_ac", "af_rgc_weighted", "ratio_to_rgc_homAlt_aaf",
                                           "exp_hom_rgc_Fdefault", "af_sas", "af_amr", "af_nfe", "af_mid",
                                           "cls_true_lof"]].to_dict("records")
    print("RGC homozigot adaylari:\n", cand[["variant_id", "hgvsc", "hgvsp", "loftee", "lof_filter", "nmd_escape",
                                             "joint_ac", "af_rgc_weighted", "ratio_to_rgc_homAlt_aaf",
                                             "exp_hom_rgc_Fdefault", "af_sas", "af_amr", "af_nfe"]].to_string(), flush=True)

    # ---------------------------------------------------------------- 4. gnomAD homozigotlari
    hom_res = {}
    syn_obs = syn_exp = 0.0
    alltabs = {}
    for sym in GENES:
        f = RAW / "gnomad" / f"{sym}_{GENES[sym]}_gnomad_r4_allvariants_pop.json"
        if not f.exists():
            print("  gnomAD all-variant dosyasi yok:", f.name, flush=True)
            continue
        t = gnomad_table(sym, "all")
        alltabs[sym] = t

    def exp_hom_gnomad(t: pd.DataFrame, fmult: float) -> np.ndarray:
        e = np.zeros(len(t))
        for p in POPS:
            n = t[f"an_{p}"].values / 2.0
            af = t[f"af_{p}"].values
            F = F_POP[p] * fmult
            e += n * ((1 - F) * af ** 2 + F * af)
        return e

    if alltabs:
        syn = pd.concat([t[(t["consequence"] == "synonymous_variant") & t["confident"] & (t["joint_af"] < 0.005)]
                         for t in alltabs.values()])
        res = {}
        for fm in (0.0, 0.5, 1.0, 2.0):
            res[str(fm)] = {"obs_hom": float(syn["joint_hom"].sum()), "exp_hom": float(exp_hom_gnomad(syn, fm).sum())}
        lam = optimize.brentq(lambda fm: exp_hom_gnomad(syn, fm).sum() - syn["joint_hom"].sum(), 0.0, 50.0) \
            if syn["joint_hom"].sum() > exp_hom_gnomad(syn, 0.0).sum() else 0.0
        S["gnomad_syn_calibration"] = {"n_variants": int(len(syn)), "by_F_multiplier": res, "F_multiplier_fit": lam}
        print("gnomAD sinonim kalibrasyonu:", S["gnomad_syn_calibration"], flush=True)
        for sym in GENES:
            t = gtabs[sym]
            for cls in ("cls_true_lof", "cls_hc_any", "cls_all_plof"):
                sub = t[t[cls]]
                hom_res[f"{sym}:{cls}"] = {
                    "obs_hom": float(sub["joint_hom"].sum()),
                    "exp_hom_F1": float(exp_hom_gnomad(sub, 1.0).sum()),
                    "exp_hom_Ffit": float(exp_hom_gnomad(sub, max(lam, 1e-9)).sum()),
                    "exp_hom_F0": float(exp_hom_gnomad(sub, 0.0).sum())}
        w = gtabs["KIF13A"][gtabs["KIF13A"]["hgvsp"] == "p.Trp1800Ter"]
        hom_res["KIF13A:p.Trp1800Ter"] = {"obs_hom": float(w["joint_hom"].sum()),
                                          "exp_hom_F1": float(exp_hom_gnomad(w, 1.0).sum()),
                                          "exp_hom_Ffit": float(exp_hom_gnomad(w, max(lam, 1e-9)).sum()),
                                          "hom_by_pop": {p: int(w[f"hom_{p}"].sum()) for p in POPS}}
    S["gnomad_hom"] = hom_res
    print("gnomAD homozigot:", json.dumps(hom_res, indent=1), flush=True)

    # ---------------------------------------------------------------- 5. Chundru 2024
    rows = read_rows(SH / "RecessiveDD_SuppDataFile" / "Sheet1.tsv")
    r1, r2, r3, r4 = rows[0], rows[1], rows[2], rows[3]
    width = max(len(r) for r in rows[:6])
    cols, a, c_, s_ = [], "", "", ""
    for i in range(width):
        x1 = r1[i].strip() if i < len(r1) else ""
        x2 = r2[i].strip() if i < len(r2) else ""
        x3 = r3[i].strip() if i < len(r3) else ""
        x4 = r4[i].strip() if i < len(r4) else ""
        if x1:
            a, c_, s_ = x1, "", ""
        if x2:
            c_, s_ = x2, ""
        if x3:
            s_ = x3
        cols.append((a, c_, s_, x4))
    data = {}
    for r in rows[4:]:
        if len(r) > 1 and r[1].strip():
            data[r[1].strip()] = r
    main_a = next(x[0] for x in cols if x[0])
    print("Chundru ana analiz:", main_a, flush=True)
    print("Chundru sinonim sutunlari:", [(i, c) for i, c in enumerate(cols) if c[0] == main_a and c[1].startswith("Synonymous")][:6], flush=True)

    idx_cache: dict = {}

    def col_index(cls_prefix, sub, stat, analysis):
        key = (cls_prefix, sub, stat, analysis)
        if key not in idx_cache:
            idx_cache[key] = next((i for i, (aa, cc, ss, tt) in enumerate(cols)
                                   if aa == analysis and cc.startswith(cls_prefix) and ss == sub and tt == stat), None)
        return idx_cache[key]

    def chv(gene, cls_prefix, sub, stat, analysis=main_a):
        r = data.get(gene)
        i = col_index(cls_prefix, sub, stat, analysis)
        if i is None and stat == "Expected":  # kaynak dosyada sinonim blogunun basligi "Exepected" yazilmis
            i = col_index(cls_prefix, sub, "Exepected", analysis)
        if r is None or i is None or i >= len(r):
            return None
        try:
            return float(r[i])
        except ValueError:
            return None

    classes = ["All", "LoF/LoF", "LoF/Functional", "Functional/Functional", "Synonymous/Synonymous"]
    ch = {}
    for g in GENES:
        ch[g] = {"known_recessive_status": data[g][2] if g in data else None}
        for cl in classes:
            ch[g][cl] = {"obs": chv(g, cl, "Total", "Observed"), "exp": chv(g, cl, "Total", "Expected"),
                         "p": chv(g, cl, "Total", "p-value")}
        ch[g]["LoF/LoF_SAS"] = {s: {"obs": chv(g, "LoF/LoF", s, "Observed"), "exp": chv(g, "LoF/LoF", s, "Expected")}
                                for s in ("SAS3", "SAS4", "SAS5")}
    # ARDD zenginlesme carpani
    st5 = read_rows(SH / "chundru_supp_tables__41588_2024_1910_MOESM4_ESM" / "Supp_Table_5.tsv")
    consensus = {r[2].strip() for r in st5 if len(r) > 5 and r[5].strip() == "Consensus"}
    discordant = {r[2].strip() for r in st5 if len(r) > 5 and r[5].strip() == "Discordant"}
    ardd_set = set(dd.loc[dd["ardd_lof"], "gene"])
    allg = []
    for g, r in data.items():
        allg.append({"gene": g, "status": r[2],
                     "lof_obs": chv(g, "LoF/LoF", "Total", "Observed"), "lof_exp": chv(g, "LoF/LoF", "Total", "Expected"),
                     "all_obs": chv(g, "All", "Total", "Observed"), "all_exp": chv(g, "All", "Total", "Expected"),
                     "syn_obs": chv(g, "Synonymous", "Total", "Observed"), "syn_exp": chv(g, "Synonymous", "Total", "Expected")})
    cg = pd.DataFrame(allg)
    cg["consensus"] = cg["gene"].isin(consensus)
    cg["ardd_lof_ddg2p"] = cg["gene"].isin(ardd_set)
    cg.to_csv(OUT / f"chundru_gene_level_{TAG}.tsv", sep="\t", index=False)

    def ratio(df, o, e):
        df = df.dropna(subset=[o, e])
        return {"n_genes": int(len(df)), "obs": float(df[o].sum()), "exp": float(df[e].sum()),
                "ratio": float(df[o].sum() / df[e].sum()) if df[e].sum() else None}

    kexp = ch["KIF13A"]["LoF/LoF"]["exp"] or 0.0
    lofcons = cg[cg["consensus"] & cg["ardd_lof_ddg2p"]]
    simil = lofcons[(lofcons["lof_exp"] >= kexp / 3) & (lofcons["lof_exp"] <= kexp * 3)]
    S["chundru"] = {
        "target_genes": ch,
        "all_genes_syn": ratio(cg, "syn_obs", "syn_exp"),
        "all_genes_lof": ratio(cg, "lof_obs", "lof_exp"),
        "consensus_lof": ratio(cg[cg["consensus"]], "lof_obs", "lof_exp"),
        "consensus_and_ddg2p_ardd_lof_lof": ratio(lofcons, "lof_obs", "lof_exp"),
        "consensus_all": ratio(cg[cg["consensus"]], "all_obs", "all_exp"),
        "similar_exp_consensus_ardd": {"n_genes": int(len(simil)),
                                       "frac_lof_obs_ge1": float((simil["lof_obs"] > 0).mean()) if len(simil) else None,
                                       "exp_range": [kexp / 3, kexp * 3]},
        "n_consensus": len(consensus), "n_discordant": len(discordant)}
    R = S["chundru"]["consensus_and_ddg2p_ardd_lof_lof"]["ratio"]
    if R and kexp:
        S["chundru"]["KIF13A_under_H_C"] = {"R_ardd_lof": R, "exp_null": kexp, "exp_H_C": R * kexp,
                                           "p0_under_H_C": math.exp(-R * kexp)}
        med = lofcons.dropna(subset=["lof_exp"])
        med = med[med["lof_exp"] > 0]
        Rm = float(np.median(med["lof_obs"] / med["lof_exp"]))
        S["chundru"]["KIF13A_under_H_C"]["R_median_gene"] = Rm
    print("Chundru:", json.dumps(S["chundru"], indent=1, default=float), flush=True)

    # ---------------------------------------------------------------- 6. popEVE
    pk = pd.read_csv(RAW / "popeve_KIF13A_NP_071396.4__NP_071396.4.csv")
    p1 = pd.read_csv(RAW / "popeve_KIF1A_NP_004312.2__NP_004312.2.csv")
    for d_ in (pk, p1):
        d_["pos"] = d_["mutant"].str.extract(r"^[A-Z](\d+)[A-Z]$")[0].astype(float)
        d_["ref"] = d_["mutant"].str[0]
    st2p = table_from(SH / "popeve_supp_tables__41588_2025_2400_MOESM3_ESM" / "SuppTable02.tsv", "chrom")
    st2p = st2p[st2p["gene"].isin(GENES)]
    pe = {"supp_table2": st2p.to_dict("records"),
          "KIF13A_n_mutants": int(len(pk)), "KIF13A_max_pos": float(pk["pos"].max()),
          "KIF1A_n_mutants": int(len(p1)), "KIF1A_max_pos": float(p1["pos"].max())}
    for name, d_ in (("KIF13A", pk), ("KIF1A", p1)):
        s = d_["popEVE"].dropna()  # bosluk orani yuksek konumlarda filtreli popEVE bos birakilmis
        pe[f"{name}_n_scored_filtered"] = int(len(s))
        pe[f"{name}_n_positions_scored"] = int(d_.loc[d_["popEVE"].notna(), "pos"].nunique())
        pe[f"{name}_scored_pos_range"] = [float(d_.loc[d_["popEVE"].notna(), "pos"].min()),
                                         float(d_.loc[d_["popEVE"].notna(), "pos"].max())]
        pe[f"{name}_frac_severe_le"] = float((s <= POPEVE_SEVERE).mean())
        pe[f"{name}_frac_severe_lt"] = float((s < POPEVE_SEVERE).mean())
        row = st2p[st2p["gene"] == name].iloc[0]
        sev, mod = float(row["%_severe"]), float(row["%_moderate"])
        pe[f"{name}_st2_pct_severe"] = sev
        pe[f"{name}_moderate_threshold_implied"] = float(np.quantile(s, sev + mod))
        pe[f"{name}_n_severe"] = int((s <= POPEVE_SEVERE).sum())
        pe[f"{name}_min_unfiltered"] = float(d_["popEVE (unfiltered)"].min())
    # konum bazinda kapsam (motor 5-352, kuyruk 520-1805)
    sc_pos = pk.loc[pk["popEVE"].notna(), "pos"]
    pe["KIF13A_scored_positions_motor_5_352"] = int(sc_pos[(sc_pos >= 5) & (sc_pos <= 352)].nunique())
    pe["KIF13A_scored_positions_tail_520_1805"] = int(sc_pos[(sc_pos >= 520)].nunique())
    t_mod = float(np.mean([pe["KIF13A_moderate_threshold_implied"], pe["KIF1A_moderate_threshold_implied"]]))
    pe["moderate_threshold_used"] = t_mod
    # katman A (gnomAD'de gorulen) ve doygunluk tablosu katman A
    cand = pd.read_csv(PRJ / "KIF13A_candidates_20260915.tsv", sep="\t")
    sat = pd.read_csv(PRJ / "KIF13A_saturation_lookup_20260916.tsv", sep="\t", low_memory=False)
    pkm = pk.set_index("mutant")["popEVE"]

    def tier_stats(df, label):
        sc = df["protein_change"].map(pkm)
        return {"set": label, "n": int(len(df)), "n_scored": int(sc.notna().sum()),
                "n_severe": int((sc <= POPEVE_SEVERE).sum()),
                "n_moderate": int(((sc > POPEVE_SEVERE) & (sc <= t_mod)).sum()),
                "median": float(sc.median()), "min": float(sc.min())}

    pe["KIF13A_tiers"] = [tier_stats(cand[cand["tier"] == t], f"candidates tier {t}") for t in "ABCD"] + \
                         [tier_stats(sat[sat["tier"] == "A"], "saturation tier A"),
                          tier_stats(sat[(sat["tier"] == "A") & (sat["snv_accessible"].astype(str) == "True")], "saturation tier A, SNV"),  # noqa: E712
                          tier_stats(sat[sat["kif1a_same_substitution_pathogenic"].astype(str) == "True"], "KIF1A same substitution pathogenic"),  # noqa: E712
                          tier_stats(sat[(sat["kif1a_same_substitution_pathogenic"].astype(str) == "True") & (sat["snv_accessible"].astype(str) == "True")], "KIF1A same substitution pathogenic, SNV"),
                          tier_stats(sat[sat["domain"].astype(str).str.contains("motor")], "all motor substitutions")]
    dn = ["D65E", "R173Q", "Q237H", "I403T", "A1146V", "S1233L", "T1315A", "S1639L", "K1122E"]
    pe["KIF13A_denovo_and_known"] = {m: (float(pkm[m]) if m in pkm.index else None) for m in dn}
    # KIF1A patojenik seti (MANE numarasi) -> kanonik (NP_004312.2 = 1690 aa)
    cal = pd.read_csv(PRJ / "KIF1A_calibration_sets_20260915.tsv", sep="\t")

    def mane_to_canon(pos):
        pos = int(pos)
        if pos <= 393:
            return pos
        if 404 <= pos <= 857:
            return pos - 9
        if pos >= 950:
            return pos - 101
        return None

    p1m = p1.set_index("mutant")["popEVE"]

    def kif1a_score(pc):
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", str(pc))
        if not m:
            return None
        cp = mane_to_canon(m.group(2))
        if cp is None:
            return None
        key = f"{m.group(1)}{cp}{m.group(3)}"
        return float(p1m[key]) if key in p1m.index else None

    cal["popEVE"] = cal["protein_change"].map(kif1a_score)
    sets = {}
    for sname in ["missense_positive_core", "missense_negative_motor_population", "missense_negative_motor_population_ac2",
                  "missense_negative_clinvar_outside_motor", "heldout_pp3_dependent_positive"]:
        sc = cal.loc[cal["set"] == sname, "popEVE"]
        sets[sname] = {"n": int((cal["set"] == sname).sum()), "n_scored": int(sc.notna().sum()),
                       "frac_severe": float((sc <= POPEVE_SEVERE).mean()) if sc.notna().any() else None,
                       "frac_severe_or_moderate": float((sc <= t_mod).mean()) if sc.notna().any() else None,
                       "median": float(sc.median()) if sc.notna().any() else None}
    pos_ = cal.loc[cal["set"] == "missense_positive_core", "popEVE"].dropna()
    neg_ = cal.loc[cal["set"] == "missense_negative_motor_population", "popEVE"].dropna()
    if len(pos_) and len(neg_):
        u = stats.mannwhitneyu(-pos_, -neg_, alternative="two-sided")
        sets["auroc_core_vs_motor_population"] = float(u.statistic / (len(pos_) * len(neg_)))
    pe["KIF1A_sets"] = sets
    cal.to_csv(OUT / f"KIF1A_calibration_popEVE_{TAG}.tsv", sep="\t", index=False)
    S["popeve"] = pe
    print("popEVE:", json.dumps(pe, indent=1, default=float)[:6000], flush=True)

    # ---------------------------------------------------------------- 6b. asiri dagilima dayanikli ampirik test
    # Nakavt sayilari gen basina Poisson'dan cok daha daginik (kurucu alleller, akraba evliligi alt kohortlari).
    # Bu yuzden KIF13A'nin gozlemi, ayni frekans bandindaki genlerin GERCEK dagilimiyla karsilastirilir:
    # toleransli genler (H0 benzeri) ve DDG2P biallelik LoF gelisimsel bozukluk genleri (H_C benzeri).
    def emp(df, col, lo, hi, k, how):
        g = df[(df[col] >= lo) & (df[col] <= hi)]
        if not len(g):
            return {"band": [float(lo), float(hi)], "n_genes": 0, "p": None}
        x = (g["ko"] == k).mean() if how == "eq" else (g["ko"] <= k).mean()
        return {"band": [float(lo), float(hi)], "n_genes": int(len(g)), "p": float(x),
                "mean_ko": float(g["ko"].mean()), "frac_ko0": float((g["ko"] == 0).mean())}

    def pair(df_t, df_a, col, center, k, how="eq", w=0.3):
        lo, hi = (1 - w) * center, (1 + w) * center
        if col == "het":
            lo, hi = math.floor(lo), math.ceil(hi)
        t_, a_ = emp(df_t, col, lo, hi, k, how), emp(df_a, col, lo, hi, k, how)
        lr = (a_["p"] / t_["p"]) if (t_["p"] and a_["p"] is not None) else None
        return {"center": float(center), "k": k, "how": how, "tolerant": t_, "ardd_lof": a_, "LR_ardd_over_tolerant": lr}

    rk = S["rgc_target_genes"]["KIF13A"]
    gk = S["gh_target_genes"]["KIF13A"]
    share_sas = gk["true_lof_share_sas_vs_cls_hc_any"]
    EMP = {
        "rgc_KIF13A_def_obs2_eq": pair(tol_r, ardd_r, "n", rk["n_hc_snv_alleles"], 2, "eq"),
        "rgc_KIF13A_def_obs2_le": pair(tol_r, ardd_r, "n", rk["n_hc_snv_alleles"], 2, "le"),
        "rgc_KIF13A_true_obs0": pair(tol_r, ardd_r, "n", rk["q_true_lof_snv_rgc_scale"] * rk["N_total"], 0, "eq"),
        "gh_KIF13A_def_obs0": pair(tol, ardd, "het", gk["het"], 0, "eq"),
        "gh_KIF13A_true_obs0": pair(tol, ardd, "het", gk["het"] * share_sas, 0, "eq"),
        "rgc_KIF1A_obs0": pair(tol_r, ardd_r, "n", S["rgc_target_genes"]["KIF1A"]["n_hc_snv_alleles"], 0, "eq"),
        "gh_KIF1A_obs0": pair(tol, ardd, "het", S["gh_target_genes"]["KIF1A"]["het"], 0, "eq", w=0.5),
        "rgc_KIF13B_obs3_eq": pair(tol_r, ardd_r, "n", S["rgc_target_genes"]["KIF13B"]["n_hc_snv_alleles"], 3, "eq"),
        "gh_KIF13B_obs0": pair(tol, ardd, "het", S["gh_target_genes"]["KIF13B"]["het"], 0, "eq"),
    }
    # birlesik olabilirlik orani (G&H ve RGC ortusmez: G&H RGC-ME kohort listesinde yok)
    def comb(*keys):
        v = [EMP[k]["LR_ardd_over_tolerant"] for k in keys]
        return float(np.prod(v)) if all(x is not None for x in v) else None
    EMP["combined_LR_scenario_rgc_homs_are_true_lof_or_unknown"] = comb("rgc_KIF13A_def_obs2_eq", "gh_KIF13A_def_obs0")
    EMP["combined_LR_scenario_rgc_homs_are_nmd_escape"] = comb("rgc_KIF13A_true_obs0", "gh_KIF13A_true_obs0")
    S["empirical_overdispersion"] = EMP
    print("Ampirik (asiri dagilim):", json.dumps(EMP, indent=1, default=float), flush=True)
    # band icindeki genlerin dagilimi (sekil/tablo icin)
    bdf = tol_r[(tol_r["n"] >= 0.7 * rk["n_hc_snv_alleles"]) & (tol_r["n"] <= 1.3 * rk["n_hc_snv_alleles"])]
    adf = ardd_r[(ardd_r["n"] >= 0.7 * rk["n_hc_snv_alleles"]) & (ardd_r["n"] <= 1.3 * rk["n_hc_snv_alleles"])]
    pd.concat([bdf.assign(group="tolerant"), adf.assign(group="ardd_lof")])[
        ["gene", "group", "n", "q", "ko", "sum_homAA", "comp_het", "mean", "loeuf_v41"]].to_csv(
        OUT / f"rgc_band_genes_KIF13A_{TAG}.tsv", sep="	", index=False)

    # ---------------------------------------------------------------- 7. ozet tablo
    rows_out = []
    g = S["gh_target_genes"]["KIF13A"]
    rows_out.append({"kohort": "Genes & Health (44.028, yetiskin, otozigot)", "tanim": "G&H pLoF (ST16), q kohort ici",
                     "q": g["q"], "F": S["gh_fit_tolerant_het10_50"]["F_eff"], "beklenen": g["exp_F_eff_het10_50"],
                     "gozlenen": g["hom"] + g["chet"]})
    rows_out.append({"kohort": "Genes & Health", "tanim": "gercek LoF (SAS payi ile)", "q": g["q"] * g["true_lof_share_sas_vs_cls_hc_any"],
                     "F": S["gh_fit_tolerant_het10_50"]["F_eff"],
                     "beklenen": g["exp_true_lof_Feff_het10_50_vs_cls_hc_any"], "gozlenen": 0})
    r_ = S["rgc_target_genes"]["KIF13A"]
    rows_out.append({"kohort": "RGC-ME (983.578)", "tanim": "RGC pLOF AAF<%1 (ST6), q = HC SNV (ST2)", "q": r_["q"],
                     "F": S["rgc_fit_tolerant_n230_920"]["F_eff"], "beklenen": r_["exp_Feff_n230_920"],
                     "gozlenen": r_["ko_obs"]})
    rows_out.append({"kohort": "RGC-ME", "tanim": "gercek LoF (gnomAD payi, RGC olcegi)", "q": r_["q_true_lof_snv_rgc_scale"],
                     "F": S["rgc_fit_tolerant_n230_920"]["F_eff"], "beklenen": r_["exp_true_lof_Feff_n230_920"],
                     "gozlenen": None})
    if "KIF13A:cls_true_lof" in S["gnomad_hom"]:
        h = S["gnomad_hom"]["KIF13A:cls_true_lof"]
        rows_out.append({"kohort": "gnomAD v4.1 (yalniz homozigot)", "tanim": "gercek LoF", "q": gq["KIF13A"]["cls_true_lof"]["q_joint"],
                         "F": "populasyona gore", "beklenen": h["exp_hom_Ffit"], "gozlenen": h["obs_hom"]})
    c = S["chundru"]["target_genes"]["KIF13A"]
    rows_out.append({"kohort": "Chundru 2024 DD triolari (29.745)", "tanim": "LoF/LoF, notr beklenti (yazarlarin)",
                     "q": None, "F": "alt gruba gore", "beklenen": c["LoF/LoF"]["exp"], "gozlenen": c["LoF/LoF"]["obs"]})
    if "KIF13A_under_H_C" in S["chundru"]:
        hc = S["chundru"]["KIF13A_under_H_C"]
        rows_out.append({"kohort": "Chundru 2024", "tanim": "LoF/LoF, H_C altinda (ARDD carpani)", "q": None, "F": None,
                         "beklenen": hc["exp_H_C"], "gozlenen": c["LoF/LoF"]["obs"]})
    for r in rows_out:
        e, o = r["beklenen"], r["gozlenen"]
        if isinstance(e, float) and o is not None:
            r["P(X<=gozlenen)"] = pois_le(o, e)
            r["P(X>=gozlenen)"] = pois_ge(int(o), e)
            r["P(0)"] = math.exp(-e)
    pd.DataFrame(rows_out).to_csv(OUT / f"KIF13A_expected_observed_{TAG}.tsv", sep="\t", index=False)
    S["summary_table"] = rows_out
    jdump(f"nakavt_summary_{TAG}.json", S)
    print(pd.DataFrame(rows_out).to_string(), flush=True)


if __name__ == "__main__":
    main()
