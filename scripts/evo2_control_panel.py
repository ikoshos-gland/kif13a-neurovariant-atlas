"""Evo 2 control panel for KIF13A c.3545-45A>G (rs41267714).

Scores a labelled panel (clinical positives, in-silico positives, common-variant negatives, target) in
both orientations, decomposes each delta log-likelihood into the site term and the downstream term, and
builds a per-position landscape around the intron 29 / exon 30 acceptor from the target's reference windows.

    mcp\\atlas_server\\.venv\\Scripts\\python.exe scripts\\evo2_control_panel.py

Credits: 4 NIM calls per panel variant (ref+alt, two orientations). The landscape reuses the target's
reference logits, so it costs nothing extra.
"""

from __future__ import annotations

import asyncio
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core, genome_models as g2, variants  # noqa: E402

core.load_dotenv()
W = 8192
HALF = W // 2
BASES = "ACGT"
COMP = str.maketrans("ACGTN", "TGCAN")
REGION = "https://rest.ensembl.org/sequence/region/human"

# (label, class, gene, GRCh38 + strand variant, HGVS on MANE)
PANEL = [
    ("KIF1A c.798+1G>A (P/LP, 2 stars)", "positive_clinical", "KIF1A", "2-240783738-C-T", "c.798+1G>A"),
    ("KIF1A c.2582+1G>A (P/LP, 2 stars)", "positive_clinical", "KIF1A", "2-240758359-C-T", "c.2582+1G>A"),
    ("KIF1A c.4319-2A>G (LP)", "positive_clinical", "KIF1A", "2-240723560-T-C", "c.4319-2A>G"),
    ("KIF1A c.609-1G>A (LP)", "positive_clinical", "KIF1A", "2-240785101-C-T", "c.609-1G>A"),
    ("KIF13A c.3545-1G>A (SpliceAI 0.99)", "positive_insilico", "KIF13A", "6-17781302-C-T", "c.3545-1G>A"),
    ("KIF13A c.3545-2A>G (SpliceAI 0.99)", "positive_insilico", "KIF13A", "6-17781303-T-C", "c.3545-2A>G"),
    ("KIF13A c.3545-45A>G rs41267714 (target)", "target", "KIF13A", "6-17781346-T-C", "c.3545-45A>G"),
]


def revcomp(s: str) -> str:
    return s[::-1].translate(COMP)


async def common_intronic_negatives(n: int = 4) -> list[tuple]:
    """Common (joint AF > 1%) KIF13A intronic SNVs 10-100 nt from an exon, excluding the target."""
    result = await variants.gnomad_variants("KIF13A", "all", save=True)
    path = ROOT / result["saved"]["path"]
    import csv
    picks = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        m = re.search(r"c\.(\d+)([+-])(\d+)([ACGT])>([ACGT])$", r["hgvsc"] or "")
        if not m or r["consequence"] != "intron_variant" or not r["joint_af"]:
            continue
        off, af = int(m.group(3)), float(r["joint_af"])
        if 10 <= off <= 100 and af > 0.01 and r["variant_id"] != "6-17781346-T-C":
            picks.append((af, r))
    picks.sort(key=lambda x: -x[0])
    return [(f"KIF13A {r['hgvsc']} AF={af:.3f}", "negative_common", "KIF13A", r["variant_id"], r["hgvsc"])
            for af, r in picks[:n]]


async def window(chrom: str, pos: int) -> str:
    start = pos - HALF + 1
    url = f"{REGION}/{chrom}:{start}..{start + W - 1}:1"
    return (await core.request("GET", url, headers={"Accept": "text/plain"})).text.strip().upper()


def site_probs(logits: np.ndarray, idx: int) -> dict[str, float]:
    lp = g2.log_softmax(logits)[idx - 1]
    return {b: float(np.exp(lp[ord(b)])) for b in BASES}


async def score(label, klass, gene, var, hgvs, keep_logits=False):
    chrom, pos, ref, alt = var.split("-")
    pos = int(pos)
    fwd = await window(chrom, pos)
    idx_f = HALF - 1
    if fwd[idx_f] != ref:
        raise ValueError(f"{label}: + strand base at {chrom}:{pos} is {fwd[idx_f]}, not {ref}")
    out = {"label": label, "class": klass, "gene": gene, "variant": var, "hgvsc": hgvs}
    kept = {}
    # both KIF1A and KIF13A are minus-strand genes: gene-sense orientation = reverse complement of + strand
    for name, seq, idx, r, a in (("plus", fwd, idx_f, ref, alt),
                                 ("sense", revcomp(fwd), W - 1 - idx_f, ref.translate(COMP), alt.translate(COMP))):
        alt_seq = seq[:idx] + a + seq[idx + 1:]
        lg_ref, lg_alt = await g2.evo2_logits(seq), await g2.evo2_logits(alt_seq)
        lp_ref, lp_alt = g2.token_logprobs(seq, lg_ref), g2.token_logprobs(alt_seq, lg_alt)
        d = lp_alt - lp_ref
        probs = site_probs(lg_ref, idx)
        top = np.argsort(-np.abs(d[idx + 1:]))[:3]
        out.update({
            f"{name}_delta": round(float(d.sum()), 3),
            f"{name}_site_term": round(float(d[idx]), 3),
            f"{name}_downstream_term": round(float(d[idx + 1:].sum()), 3),
            f"{name}_max_abs_upstream_change": float(np.abs(d[:idx]).max()),
            f"{name}_P_ref_at_site": round(probs[r], 4),
            f"{name}_P_alt_at_site": round(probs[a], 4),
            f"{name}_top_downstream_offsets": [f"+{int(k) + 1}:{d[idx + 1 + k]:.2f}" for k in top],
        })
        if keep_logits:
            kept[name] = (seq, idx, lg_ref)
    return out, kept


def landscape(kept: dict, exon_len: int, span_up: int = 60, span_down: int = 20) -> list[dict]:
    """Per-position expectation around the acceptor from reference logits only.
    sense: upstream (5' in the pre-mRNA) context; plus: the opposite context. Offsets are relative to c.3545-45."""
    seq_s, idx_s, lg_s = kept["sense"]
    seq_p, _idx_p, lg_p = kept["plus"]
    lps, lpp = g2.log_softmax(lg_s), g2.log_softmax(lg_p)
    rows = []
    for k in range(-span_up + 45, 45 + span_down):  # k = offset from the target in sense coordinates
        j = idx_s + k
        off = k - 45  # relative to c.3545
        if off >= exon_len:
            break
        label = f"c.3545{off}" if off < 0 else f"c.{3545 + off}"
        base = seq_s[j]
        p_sense = float(np.exp(lps[j - 1][ord(base)]))
        jp = W - 1 - j  # same nucleotide on the + strand, predicted from the other side
        p_plus = float(np.exp(lpp[jp - 1][ord(base.translate(COMP))]))
        rows.append({"position": label, "region": "intron" if off < 0 else "exon", "ref_sense": base,
                     "P_ref_from_5prime_context": round(p_sense, 4), "P_ref_from_3prime_context": round(p_plus, 4),
                     "P_ref_mean": round((p_sense + p_plus) / 2, 4), "is_target": off == -45})
    return rows


async def exon30_length() -> int:
    data = await core.get_json("https://rest.ensembl.org/lookup/id/ENST00000259711", params={"expand": 1})
    for ex in data["Exon"]:
        if ex["end"] == 17781301:  # c.3545 is the first exonic base (minus strand: exon end coordinate)
            return ex["end"] - ex["start"] + 1
    raise RuntimeError("exon starting at c.3545 not found")


async def main() -> int:
    panel = PANEL[:6] + await common_intronic_negatives() + PANEL[6:]
    rows, kept = [], {}
    for item in panel:
        row, k = await score(*item, keep_logits=item[1] == "target")
        kept.update(k)
        rows.append(row)
        print(f"{row['class']:<18} {row['label']:<42} sense {row['sense_delta']:>8} (site {row['sense_site_term']:>7}, "
              f"down {row['sense_downstream_term']:>7})  plus {row['plus_delta']:>8}  "
              f"P_ref/P_alt sense {row['sense_P_ref_at_site']}/{row['sense_P_alt_at_site']}", flush=True)
    exon_len = await exon30_length()
    land = landscape(kept, exon_len)
    prov = dict(source="NVIDIA NIM arc/evo2-40b forward (unembed) on Ensembl GRCh38 windows",
                urls=[g2.NIM_FORWARD, REGION], tool="scripts/evo2_control_panel.py", params={"window": W})
    print(core.save_tsv("evo2", f"control_panel_rs41267714_w{W}_{core.today()}.tsv", rows, base=core.INTERIM, **prov))
    print(core.save_tsv("evo2", f"acceptor_landscape_intron29_{core.today()}.tsv", land, base=core.INTERIM, **prov))
    print("exon 30 length:", exon_len)
    for r in land:
        flag = "  <== target" if r["is_target"] else ""
        print(f"  {r['position']:>10} {r['ref_sense']} {r['P_ref_from_5prime_context']:.3f} {r['P_ref_from_3prime_context']:.3f}{flag}")
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
