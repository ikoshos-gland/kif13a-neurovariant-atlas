"""Place every KIF13A motor-domain residue in 3D: ATP site, microtubule interface, burial, pLDDT.

The KIF13A motor domain has no experimental structure, but its paralog KIF1A does, bound to a
microtubule: RCSB 8UTN and siblings, from Benoit et al., Nat Commun 2024;15:5530,
doi 10.1038/s41467-024-48720-4 (the primary citation is read back out of the deposited file).
So we borrow the geometry:

  1. parse the KIF1A[1-393] + alpha/beta tubulin complex, find the KIF1A chain by its UniProt
     accession (Q12756) in the mmCIF struct_ref block, plus its nucleotide ligand and the tubulin chains
  2. superpose the KIF13A AlphaFold model (AF-Q9H1H9-F1) onto that KIF1A chain with CA atoms of the
     residue pairs in the curated paralog map, trimming divergent loops iteratively
  3. for every KIF13A motor residue measure the minimum heavy-atom distance to the nucleotide and to
     tubulin, its relative solvent accessibility (Shrake-Rupley on the isolated motor domain) and its pLDDT
  4. join that onto the candidate list and give each variant a structural class

Nothing here is re-downloaded if it is already in data/raw (see data/raw/PROVENANCE.jsonl).

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\structure_axis.py
    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\structure_axis.py --pdb-id 8UTS --kif1a-chain K
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core  # noqa: E402

from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # noqa: E402
from Bio.PDB.MMCIFParser import MMCIFParser  # noqa: E402
from Bio.PDB.SASA import ShrakeRupley  # noqa: E402
from Bio.PDB.Superimposer import Superimposer  # noqa: E402

RCSB_FILES = "https://files.rcsb.org/download"
KIF1A_UNP = "Q12756"          # canonical; positions 1-393 are identical in MANE Q12756-3
KIF13A_UNP = "Q9H1H9"         # MANE, 1805 aa
AF_DIR = "alphafold/AF-Q9H1H9-F1_v6"
PARALOG_MAP = "paralog_maps/KIF1A_to_KIF13A_pairwise_MANE.tsv"

# KIF13A domains, MANE numbering (same table as scripts/build_kif13a_candidates.py)
MOTOR_START, MOTOR_END, NECK_END = 5, 352, 380
PAIR_MAX_KIF1A = 365          # the cryo-EM construct is KIF1A[1-393]; stop before the disordered end
NECK_CLASS = "boyun baglayici (yapi ekseni guvenilmez)"
OUTSIDE_CLASS = "motor disi (yapi ekseni yok)"

# Tien et al. 2013 (PLoS ONE 8:e80635) theoretical maximum accessible surface area, A^2
MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLU": 223, "GLN": 225, "GLY": 104,
           "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
           "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}

# KIF1A residues with well-known mechanistic roles, used as sanity checks (KIF1A numbering)
SANITY = {99: "P-loop / switch I bolgesi, ClinVar P/LP", 254: "ClinVar P/LP sicak nokta",
          305: "L11-alpha4, mikrotubul yakini, ClinVar P/LP", 350: "kinezin-tubulin tuz koprusu (literatur)",
          364: "boyun baglayici, ClinVar P/LP"}


def table(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def latest(pattern: str, base: Path) -> Path:
    return sorted(base.glob(pattern))[-1]


async def ensure_cif(pdb_id: str) -> Path:
    """Download <pdb_id>.cif into data/raw/pdb once, with a provenance sidecar."""
    path = core.RAW / "pdb" / f"{pdb_id}.cif"
    if path.exists():
        print(f"using cached {core.rel(path)} ({path.stat().st_size} bytes)")
        return path
    url = f"{RCSB_FILES}/{pdb_id}.cif"
    content = await core.get_bytes(url, timeout=600)
    rec = core.save_bytes("pdb", f"{pdb_id}.cif", content, source="RCSB PDB", urls=[url],
                          tool="scripts/structure_axis.py", params={"pdb_id": pdb_id})
    print(f"downloaded {rec['path']} ({rec['bytes']} bytes, sha256 {rec['sha256']})")
    return path


def describe_entities(cif: Path) -> dict:
    """Map chains to entities using the mmCIF metadata, so nothing about the complex is assumed."""
    d = MMCIF2Dict(str(cif))
    ent_desc = dict(zip(d["_entity.id"], d["_entity.pdbx_description"]))
    poly_strands = dict(zip(d["_entity_poly.entity_id"], d["_entity_poly.pdbx_strand_id"]))
    accession = dict(zip(d["_struct_ref.entity_id"], d["_struct_ref.pdbx_db_accession"]))
    chains: dict[str, dict] = {}
    for ent, strands in poly_strands.items():
        for ch in strands.split(","):
            chains[ch.strip()] = {"entity": ent, "description": ent_desc.get(ent, ""),
                                  "uniprot": accession.get(ent, "")}
    citation = {}
    for i, cid in enumerate(d.get("_citation.id", [])):
        if cid == "primary":
            citation = {"title": d["_citation.title"][i], "journal": d["_citation.journal_abbrev"][i],
                        "year": d["_citation.year"][i], "doi": d["_citation.pdbx_database_id_DOI"][i],
                        "pmid": d["_citation.pdbx_database_id_PubMed"][i]}
    return {"title": d["_struct.title"][0], "chains": chains, "citation": citation,
            "hetero": sorted({(a, c) for a, c, g in zip(d["_atom_site.auth_asym_id"], d["_atom_site.auth_comp_id"],
                                                        d["_atom_site.group_PDB"]) if g == "HETATM"})}


def heavy(residue) -> np.ndarray:
    return np.array([a.coord for a in residue if a.element != "H"], dtype=float)


def min_distance(coords: np.ndarray, target: np.ndarray, chunk: int = 4000) -> tuple[float, int]:
    """Minimum distance from a small atom set to a large one, and the index of the nearest target atom."""
    best, best_i = float("inf"), -1
    for start in range(0, len(target), chunk):
        block = target[start:start + chunk]
        d = np.linalg.norm(coords[:, None, :] - block[None, :, :], axis=2)
        i = int(np.argmin(d))
        if float(d.flat[i]) < best:
            best, best_i = float(d.flat[i]), start + i % len(block)
    return best, best_i


def superpose(pairs: list[tuple], trim: float, rounds: int = 6) -> tuple[Superimposer, list[tuple], list[float]]:
    """Least-squares CA fit, then drop pairs deviating more than `trim` and refit (divergent loops)."""
    def fit(items):
        sup = Superimposer()
        sup.set_atoms([p[2] for p in items], [p[3] for p in items])
        return sup

    current, history = pairs, []
    for _ in range(rounds):
        sup = fit(current)
        history.append(sup.rms)
        rot, tran = sup.rotran
        moved = np.array([p[3].coord for p in current]) @ rot + tran
        dev = np.linalg.norm(moved - np.array([p[2].coord for p in current]), axis=1)
        keep = [p for p, v in zip(current, dev) if v <= trim]
        if len(keep) == len(current) or len(keep) < 100:
            break
        current = keep
    return fit(current), current, history


def build_figure(residues: list[dict], candidates: list[dict], path: Path, atp_cut: float, mt_cut: float) -> None:
    """Two panels, hand-written SVG (this environment has no plotting library).

    a) position track 1-380: one bar per residue coloured by structural class, candidate residues as lollipops
    b) nucleotide distance vs tubulin distance: one dot per residue, candidate residues drawn on top
    """
    W, H = 1100, 780
    L, R = 70, 26
    PW = W - L - R
    BAR_Y, BAR_H = 120, 20                # track panel
    T2, H2 = 235, 400                     # scatter panel
    CAP = 40.0                            # distances are clamped for display
    COL = {"ATP cebi": "#d1495b", "mikrotubul arayuzu": "#2a78d6", "gomulu cekirdek": "#7d7c75",
           "yuzey": "#d8d7d0", NECK_CLASS: "#efd98a"}
    TIER = {"A": "#c81e1e", "B": "#eb6834", "C": "#2a78d6", "D": "#6250d6"}
    SHORT = {NECK_CLASS: "boyun baglayici (guvenilmez)"}

    def xp(pos: float) -> float:
        return L + PW * (pos - 1) / (NECK_END - 1)

    def xs(v: float) -> float:
        return L + PW * min(v, CAP) / CAP

    def ys(v: float) -> float:
        return T2 + H2 * (1 - min(v, CAP) / CAP)

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
         'font-family="Arial, Helvetica, sans-serif">', f'<rect width="{W}" height="{H}" fill="#ffffff"/>']
    o.append(f'<text x="{L}" y="26" font-size="14" fill="#0b0b0b">KIF13A motor alani yapi ekseni: '
             'AlphaFold modeli, KIF1A-mikrotubul kompleksine bindirildi</text>')

    # ---------------------------------------------------------------- a) position track
    o.append(f'<text x="{L}" y="52" font-size="11" fill="#52514e">a) kalinti sinifi boyunca MANE pozisyonu '
             f'1-{NECK_END}; cubuklarin ustundeki lolipoplar aday kalintilar (A ve B katmanlari)</text>')
    step = PW / (NECK_END - 1)
    for r in residues:
        o.append(f'<rect x="{xp(r["protein_pos"]):.2f}" y="{BAR_Y}" width="{step + 0.6:.2f}" '
                 f'height="{BAR_H}" fill="{COL[r["structural_class"]]}"/>')
    o.append(f'<rect x="{xp(MOTOR_START):.1f}" y="{BAR_Y}" width="{xp(MOTOR_END) - xp(MOTOR_START):.1f}" '
             f'height="{BAR_H}" fill="none" stroke="#0b0b0b" stroke-width="0.9"/>')
    o.append(f'<text x="{(xp(MOTOR_START) + xp(MOTOR_END)) / 2:.0f}" y="{BAR_Y + BAR_H + 15}" font-size="10" '
             f'fill="#52514e" text-anchor="middle">motor alan {MOTOR_START}-{MOTOR_END}</text>')
    o.append(f'<text x="{xp((MOTOR_END + NECK_END) / 2):.0f}" y="{BAR_Y + BAR_H + 15}" font-size="10" '
             'fill="#8a7a30" text-anchor="middle">boyun</text>')
    best: dict[int, str] = {}
    for c in candidates:
        t = c["tier"]
        if t in TIER and ("ABCD".index(t) < "ABCD".index(best.get(c["protein_pos"], "D")) or c["protein_pos"] not in best):
            best[c["protein_pos"]] = t
    for pos, tier in sorted(best.items()):
        if tier not in ("A", "B"):
            continue
        x, top = xp(pos), BAR_Y - (26 if tier == "A" else 14)
        o.append(f'<line x1="{x:.1f}" y1="{BAR_Y}" x2="{x:.1f}" y2="{top}" stroke="{TIER[tier]}" stroke-width="1"/>')
        o.append(f'<circle cx="{x:.1f}" cy="{top}" r="3" fill="{TIER[tier]}"/>')
    for pos in (1, 50, 100, 150, 200, 250, 300, 350, NECK_END):
        o.append(f'<text x="{xp(pos):.1f}" y="{BAR_Y + BAR_H + 30}" font-size="10" fill="#898781" '
                 f'text-anchor="middle">{pos}</text>')

    # ---------------------------------------------------------------- b) distance scatter
    o.append(f'<text x="{L}" y="{T2 - 14}" font-size="11" fill="#52514e">b) her kalinti icin nukleotide '
             f'(ANP + Mg) ve tubuline en kisa agir-atom uzakligi; {CAP:.0f} A ustu kirpildi, etiketler '
             'A katmani kalintilari</text>')
    o.append(f'<rect x="{L}" y="{T2}" width="{PW}" height="{H2}" fill="#fcfcfa" stroke="#dedcd4"/>')
    for v in (0, 5, 10, 20, 30, 40):
        o.append(f'<line x1="{L}" y1="{ys(v):.1f}" x2="{L + PW}" y2="{ys(v):.1f}" stroke="#efedE6"/>')
        o.append(f'<text x="{L - 8}" y="{ys(v) + 4:.1f}" font-size="10" fill="#898781" text-anchor="end">{v}</text>')
        o.append(f'<line x1="{xs(v):.1f}" y1="{T2}" x2="{xs(v):.1f}" y2="{T2 + H2}" stroke="#efede6"/>')
        o.append(f'<text x="{xs(v):.1f}" y="{T2 + H2 + 15}" font-size="10" fill="#898781" text-anchor="middle">{v}</text>')
    o.append(f'<line x1="{xs(atp_cut):.1f}" y1="{T2}" x2="{xs(atp_cut):.1f}" y2="{T2 + H2}" stroke="#d1495b" '
             'stroke-width="1.1" stroke-dasharray="5 4"/>')
    o.append(f'<line x1="{L}" y1="{ys(mt_cut):.1f}" x2="{L + PW}" y2="{ys(mt_cut):.1f}" stroke="#2a78d6" '
             'stroke-width="1.1" stroke-dasharray="5 4"/>')
    o.append(f'<text x="{xs(atp_cut) + 6:.1f}" y="{T2 + 14}" font-size="10" fill="#d1495b">ATP esigi {atp_cut:.0f} A</text>')
    o.append(f'<text x="{L + PW - 6}" y="{ys(mt_cut) - 6:.1f}" font-size="10" fill="#2a78d6" text-anchor="end">'
             f'mikrotubul esigi {mt_cut:.0f} A</text>')
    for r in residues:
        o.append(f'<circle cx="{xs(r["dist_nucleotide_A"]):.1f}" cy="{ys(r["dist_tubulin_A"]):.1f}" r="2" '
                 f'fill="{COL[r["structural_class"]]}" opacity="0.55"/>')
    by_pos = {r["protein_pos"]: r for r in residues}
    # one marker per residue (several substitutions share a position), strongest tier wins
    marks = []
    for pos, tier in sorted(best.items()):
        r = by_pos.get(pos)
        if r is None:
            continue
        marks.append((xs(r["dist_nucleotide_A"]), ys(r["dist_tubulin_A"]), tier, f'{r["aa"].capitalize()}{pos}'))
    for x, y, tier, _n in marks:
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{TIER[tier]}" stroke="#ffffff" stroke-width="0.9"/>')
    # greedy label placement: first free slot above or below, otherwise the label is dropped
    boxes: list[tuple[float, float, float, float]] = []
    for x, y, tier, name in sorted((m for m in marks if m[2] == "A"), key=lambda m: (m[1], m[0])):
        w = 6.3 * len(name)
        for dy in (-8, 12, -19, 23, -30, 34, -41, 45):
            bx, by = x - w / 2, y + dy - 9
            if all(bx > b[2] or bx + w < b[0] or by > b[3] or by + 11 < b[1] for b in boxes):
                boxes.append((bx, by, bx + w, by + 11))
                o.append(f'<text x="{x:.1f}" y="{y + dy:.1f}" font-size="10" fill="#0b0b0b" '
                         f'text-anchor="middle">{name}</text>')
                break
    o.append(f'<text x="{L + PW / 2:.0f}" y="{T2 + H2 + 34}" font-size="11" fill="#52514e" text-anchor="middle">'
             'nukleotide uzaklik, A</text>')
    o.append(f'<text x="18" y="{T2 + H2 / 2:.0f}" font-size="11" fill="#52514e" text-anchor="middle" '
             f'transform="rotate(-90 18 {T2 + H2 / 2:.0f})">tubuline uzaklik, A</text>')

    # ---------------------------------------------------------------- legends, wrapped
    lx, ly = L, H - 60
    for name, colour in COL.items():
        label = f"{SHORT.get(name, name)} ({sum(1 for r in residues if r['structural_class'] == name)})"
        width = 14 + 6.4 * len(label) + 16
        if lx + width > W - R:
            lx, ly = L, ly + 19
        o.append(f'<rect x="{lx}" y="{ly - 8}" width="9" height="9" fill="{colour}" stroke="#b6b5ae" stroke-width="0.5"/>')
        o.append(f'<text x="{lx + 14}" y="{ly}" font-size="11" fill="#52514e">{label}</text>')
        lx += width
    lx, ly = L, ly + 21
    for tier, colour in TIER.items():
        label = f"katman {tier} ({sum(1 for t in best.values() if t == tier)} kalinti)"
        width = 14 + 6.4 * len(label) + 16
        if lx + width > W - R:
            lx, ly = L, ly + 19
        o.append(f'<circle cx="{lx + 4}" cy="{ly - 4}" r="4" fill="{colour}"/>')
        o.append(f'<text x="{lx + 14}" y="{ly}" font-size="11" fill="#52514e">{label}</text>')
        lx += width
    o.append("</svg>")
    path.write_text("\n".join(o), encoding="utf-8")
    print(f"{core.rel(path)} {path.stat().st_size} bytes")


async def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-id", default="8UTN", help="KIF1A-microtubule cryo-EM entry (default 8UTN)")
    ap.add_argument("--kif1a-chain", default="auto", help="chain to superpose onto, or 'auto' for the best fit")
    ap.add_argument("--trim", type=float, default=3.0, help="CA deviation (A) above which a pair is dropped")
    ap.add_argument("--atp-cutoff", type=float, default=5.0, help="nucleotide contact distance (A)")
    ap.add_argument("--mt-cutoff", type=float, default=5.0, help="tubulin contact distance (A)")
    ap.add_argument("--buried-cutoff", type=float, default=0.15, help="relative SASA below which a residue is core")
    ap.add_argument("--max-rmsd", type=float, default=3.0, help="abort if the trimmed fit is worse than this")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()

    print(f"structure axis | {core.today()} | pdb {args.pdb_id} | trim {args.trim} A | "
          f"cutoffs ATP {args.atp_cutoff} A, MT {args.mt_cutoff} A, buried relSASA {args.buried_cutoff}")

    # ------------------------------------------------------------------ 1. the cryo-EM complex
    cif = await ensure_cif(args.pdb_id)
    meta = describe_entities(cif)
    citation = meta["citation"]
    print(f"title: {meta['title']}")
    if citation:
        print(f"primary citation (from the deposited file): {citation['journal']} {citation['year']}, "
              f"doi {citation['doi']}, PMID {citation['pmid']}")
    kif1a_chains = [c for c, v in meta["chains"].items() if core.unversioned(v["uniprot"]) == KIF1A_UNP]
    tubulin_chains = [c for c, v in meta["chains"].items() if "tubulin" in v["description"].lower()]
    if not kif1a_chains:
        print(f"ABORT: no chain in {args.pdb_id} maps to UniProt {KIF1A_UNP}; this entry is not KIF1A.")
        return 1
    for c, v in sorted(meta["chains"].items()):
        print(f"  chain {c}: entity {v['entity']} {v['description']} [{v['uniprot'] or 'no UNP'}]")
    print(f"  KIF1A chains: {kif1a_chains} | tubulin chains: {tubulin_chains}")
    print(f"  heteroatoms: {meta['hetero']}")

    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure(args.pdb_id, str(cif))[0]
    tub_coords, tub_label = [], []
    for cid in tubulin_chains:
        kind = "alpha" if "alpha" in meta["chains"][cid]["description"].lower() else "beta"
        for res in model[cid]:
            if res.id[0] != " ":
                continue
            for atom in res:
                if atom.element != "H":
                    tub_coords.append(atom.coord)
                    tub_label.append(f"{cid}/{kind}/{res.resname}{res.id[1]}")
    tub = np.array(tub_coords, dtype=float)
    print(f"  tubulin heavy atoms: {len(tub)}")

    # ------------------------------------------------------------------ 2. the AlphaFold model and the map
    af_cif = core.RAW / AF_DIR / "AF-Q9H1H9-F1-model_v6.cif"
    conf = json.loads((core.RAW / AF_DIR / "AF-Q9H1H9-F1-confidence_v6.json").read_text(encoding="utf-8"))
    plddt = dict(zip(conf["residueNumber"], conf["confidenceScore"]))
    af_model = parser.get_structure("kif13a", str(af_cif))[0]
    af_chain = list(af_model)[0]
    af_res = {r.id[1]: r for r in af_chain if r.id[0] == " "}
    print(f"KIF13A AlphaFold model: {core.rel(af_cif)}, {len(af_res)} residues, "
          f"pLDDT for {len(plddt)} residues")

    pmap = [(int(r["KIF1A_pos"]), int(r["KIF13A_pos"]), r["KIF1A_aa"], r["KIF13A_aa"])
            for r in table(core.INTERIM / PARALOG_MAP)]

    fits = {}
    for cid in kif1a_chains:
        kres = {r.id[1]: r for r in model[cid] if r.id[0] == " "}
        pairs = [(k1, k13, kres[k1]["CA"], af_res[k13]["CA"], a1 == a13)
                 for k1, k13, a1, a13 in pmap
                 if k1 <= PAIR_MAX_KIF1A and k1 in kres and k13 in af_res
                 and "CA" in kres[k1] and "CA" in af_res[k13]]
        sup, kept, history = superpose(pairs, args.trim)
        fits[cid] = {"sup": sup, "kept": kept, "pairs": pairs, "raw_rmsd": history[0]}
        print(f"  superposition on chain {cid}: {len(pairs)} CA pairs, RMSD {history[0]:.2f} A "
              f"-> trimmed {len(kept)} pairs, RMSD {sup.rms:.2f} A "
              f"({sum(1 for p in pairs if p[4])} identical residues)")

    chain_id = args.kif1a_chain if args.kif1a_chain != "auto" else min(fits, key=lambda c: fits[c]["sup"].rms)
    if chain_id not in fits:
        print(f"ABORT: chain {chain_id} is not a KIF1A chain in {args.pdb_id} (KIF1A chains: {kif1a_chains})")
        return 1
    fit = fits[chain_id]
    rot, tran = fit["sup"].rotran
    print(f"using KIF1A chain {chain_id} as the primary frame: trimmed RMSD {fit['sup'].rms:.2f} A over "
          f"{len(fit['kept'])} CA pairs (all pairs {fit['raw_rmsd']:.2f} A over {len(fit['pairs'])})")
    if fit["sup"].rms > args.max_rmsd:
        print(f"ABORT: trimmed RMSD {fit['sup'].rms:.2f} A is above --max-rmsd {args.max_rmsd}; "
              "the superposition is not trustworthy, so no distances are reported.")
        return 1
    # A two-heads-bound complex holds two genuinely microtubule-bound head states. They touch tubulin
    # with slightly different residues, so the tubulin distance is the minimum over every head frame
    # whose own fit is acceptable, and each frame is also kept as its own column.
    mt_frames = [c for c in fits if fits[c]["sup"].rms <= args.max_rmsd]
    frame_rmsd = ", ".join("%s:%.2f" % (c, fits[c]["sup"].rms) for c in mt_frames)
    print(f"  tubulin distances are taken over head frames {mt_frames} (trimmed RMSD {frame_rmsd} A)")

    # the nucleotide and Mg sitting in the chosen head
    lig_res = [r for r in model[chain_id] if r.id[0] != " "]
    lig = np.array([a.coord for r in lig_res for a in r if a.element != "H"], dtype=float)
    lig_names = "+".join(sorted({r.resname for r in lig_res}))
    if not len(lig):
        print(f"ABORT: chain {chain_id} carries no nucleotide ligand; pick a chain or entry that does.")
        return 1
    print(f"  nucleotide ligand in chain {chain_id}: {lig_names} ({len(lig)} heavy atoms)")

    # does this head actually touch the microtubule? (sanity: it must, this is an MT-bound complex)
    kif1a_contact = 0
    for res in model[chain_id]:
        if res.id[0] != " ":
            continue
        d, _ = min_distance(heavy(res), tub)
        kif1a_contact += d <= args.mt_cutoff
    print(f"  KIF1A chain {chain_id} residues within {args.mt_cutoff} A of tubulin: {kif1a_contact}")

    # ------------------------------------------------------------------ 3. per-residue geometry
    # SASA on the isolated motor domain: the full 1805 aa model would let the low-confidence tail
    # drape over the motor and fake burial.
    sasa_model = parser.get_structure("kif13a_motor", str(af_cif))[0]
    sasa_chain = list(sasa_model)[0]
    for res in [r for r in sasa_chain if r.id[0] != " " or r.id[1] > NECK_END]:
        sasa_chain.detach_child(res.id)
    ShrakeRupley().compute(sasa_model, level="R")
    sasa = {r.id[1]: r.sasa for r in sasa_chain}

    paired = {p[1]: p[0] for p in fit["pairs"]}
    kept_pairs = {p[1] for p in fit["kept"]}
    cb = {pos: (r["CB"] if "CB" in r else r["CA"]).coord for pos, r in af_res.items() if "CA" in r}
    cb_pos = sorted(p for p in cb if p <= NECK_END)
    cb_arr = np.array([cb[p] for p in cb_pos], dtype=float)

    residues: list[dict] = []
    for pos in range(1, NECK_END + 1):
        res = af_res.get(pos)
        if res is None:
            continue
        atoms = heavy(res)
        coords = atoms @ rot + tran
        d_lig, _ = min_distance(coords, lig)
        per_frame = {}
        for cid in mt_frames:
            r2, t2 = fits[cid]["sup"].rotran
            d, i = min_distance(atoms @ r2 + t2, tub)
            per_frame[cid] = (d, tub_label[i])
        best_frame = min(per_frame, key=lambda c: per_frame[c][0])
        d_tub, near_tub = per_frame[best_frame]
        rel = sasa[pos] / MAX_ASA.get(res.resname, 200.0)
        # The neck linker leaves the head: its AlphaFold path is not the cryo-EM docked/undocked path,
        # so its distances are printed but never turned into a claim.
        if pos > MOTOR_END:
            klass = NECK_CLASS
        elif d_lig <= args.atp_cutoff:
            klass = "ATP cebi"
        elif d_tub <= args.mt_cutoff:
            klass = "mikrotubul arayuzu"
        elif rel < args.buried_cutoff:
            klass = "gomulu cekirdek"
        else:
            klass = "yuzey"
        neighbours = int(np.sum(np.linalg.norm(cb_arr - cb[pos], axis=1) <= 10.0) - 1)
        row = {
            "protein_pos": pos, "aa": res.resname, "domain": "motor" if MOTOR_START <= pos <= MOTOR_END else
            ("N_term" if pos < MOTOR_START else "neck_linker"),
            "structural_class": klass,
            "dist_nucleotide_A": round(d_lig, 2), "dist_tubulin_A": round(d_tub, 2),
            "tubulin_frame": best_frame, "nearest_tubulin": near_tub,
            "sasa_A2": round(sasa[pos], 1), "rel_sasa": round(rel, 3), "neighbours_10A": neighbours,
            "plddt": plddt.get(pos), "kif1a_pos": paired.get(pos),
            "superpose_pair": pos in paired, "superpose_pair_kept": pos in kept_pairs,
        }
        for cid in mt_frames:
            row[f"dist_tubulin_chain{cid}_A"] = round(per_frame[cid][0], 2)
        residues.append(row)

    summary_counts = Counter(r["structural_class"] for r in residues if r["domain"] == "motor")
    print("\nmotor domain residues by class:", summary_counts.most_common())

    # ------------------------------------------------------------------ 4. sanity checks
    kif13a_of = {k1: k13 for k1, k13, _a, _b in pmap}
    checks = []
    clashes = [r["protein_pos"] for r in residues if r["dist_tubulin_A"] < 2.0]
    print(f"\nsuperposed KIF13A residues closer than 2 A to tubulin (rigid-body rotamer clashes): "
          f"{len(clashes)} {clashes}")
    print("\nsanity checks: the real KIF1A residue in the cryo-EM structure, then its mapped KIF13A residue")
    by_pos = {r["protein_pos"]: r for r in residues}
    for k1, note in SANITY.items():
        k13 = kif13a_of.get(k1)
        r = by_pos.get(k13) if k13 else None
        native: dict[str, float | str] = {}
        for cid in kif1a_chains:
            nres = {x.id[1]: x for x in model[cid] if x.id[0] == " "}.get(k1)
            if nres is None:
                continue
            n_lig, _ = min_distance(heavy(nres), lig)
            n_tub, i = min_distance(heavy(nres), tub)
            native[cid] = {"aa": nres.resname, "dist_nucleotide_A": round(n_lig, 2),
                           "dist_tubulin_A": round(n_tub, 2), "nearest_tubulin": tub_label[i]}
            print(f"  KIF1A {nres.resname}{k1} in head {cid}: ATP {n_lig:.1f} A, MT {n_tub:.1f} A ({tub_label[i]})")
        if r is None:
            print(f"    -> KIF13A {k13}: no structural row")
            checks.append({"kif1a_pos": k1, "kif13a_pos": k13, "note": note, "result": "eslesme yok",
                           "native_per_head": native})
            continue
        print(f"    -> KIF13A {r['aa']}{k13}: ATP {r['dist_nucleotide_A']:.1f} A, "
              f"MT {r['dist_tubulin_A']:.1f} A ({r['nearest_tubulin']}, head {r['tubulin_frame']}), "
              f"relSASA {r['rel_sasa']:.2f}, pLDDT {r['plddt']:.0f} -> {r['structural_class']}  [{note}]")
        checks.append({"kif1a_pos": k1, "kif13a_pos": k13, "aa": r["aa"], "note": note,
                       "native_per_head": native,
                       "dist_nucleotide_A": r["dist_nucleotide_A"], "dist_tubulin_A": r["dist_tubulin_A"],
                       "tubulin_frame": r["tubulin_frame"], "nearest_tubulin": r["nearest_tubulin"],
                       "rel_sasa": r["rel_sasa"], "plddt": r["plddt"], "result": r["structural_class"]})

    # ------------------------------------------------------------------ 5. join onto the candidates
    cand_path = latest("KIF13A_candidates_*.tsv", ROOT / "data" / "processed")
    cands = table(cand_path)
    print(f"\ncandidates: {core.rel(cand_path)} ({len(cands)} rows)")
    joined = []
    for c in cands:
        pos = int(c["protein_pos"])
        r = by_pos.get(pos)
        joined.append({
            "tier": c["tier"], "protein_change": c["protein_change"], "protein_pos": pos,
            "domain": c["domain"], "am_pathogenicity": c["am_pathogenicity"],
            "confidently_observed": c["confidently_observed"], "kif1a_pos": c["kif1a_pos"],
            "kif1a_pathogenic_at_residue": c["kif1a_pathogenic_at_residue"],
            "structural_class": r["structural_class"] if r else OUTSIDE_CLASS,
            "dist_nucleotide_A": r["dist_nucleotide_A"] if r else None,
            "dist_tubulin_A": r["dist_tubulin_A"] if r else None,
            "tubulin_frame": r["tubulin_frame"] if r else None,
            "nearest_tubulin": r["nearest_tubulin"] if r else None,
            "rel_sasa": r["rel_sasa"] if r else None, "neighbours_10A": r["neighbours_10A"] if r else None,
            "plddt": r["plddt"] if r else None,
            "structure_evidence_ok": bool(r and r["plddt"] and r["plddt"] >= 70 and r["superpose_pair"]),
            "variant_id": c["variant_id"], "hgvsc": c["hgvsc"],
        })
    joined.sort(key=lambda r: ("ABCD-".index(r["tier"]), -(float(r["am_pathogenicity"] or 0))))

    stamp = core.today()
    out_res = core.save_tsv("structure", f"KIF13A_motor_residue_structure_{stamp}.tsv", residues,
                            base=core.INTERIM,
                            source=f"derived: RCSB {args.pdb_id} + AlphaFold AF-{KIF13A_UNP}-F1_v6 + paralog map",
                            urls=[f"{RCSB_FILES}/{args.pdb_id}.cif"], tool="scripts/structure_axis.py",
                            params={"pdb_id": args.pdb_id, "kif1a_chain": chain_id, "trim_A": args.trim,
                                    "atp_cutoff_A": args.atp_cutoff, "mt_cutoff_A": args.mt_cutoff,
                                    "buried_rel_sasa": args.buried_cutoff,
                                    "trimmed_rmsd_A": round(fit["sup"].rms, 3), "n_pairs": len(fit["kept"])})
    out_cand = core.save_tsv("structure", f"KIF13A_candidates_structure_{stamp}.tsv", joined, base=core.INTERIM,
                             source=f"derived: {core.rel(cand_path)} + structural axis",
                             urls=[f"{RCSB_FILES}/{args.pdb_id}.cif"], tool="scripts/structure_axis.py",
                             params={"pdb_id": args.pdb_id, "kif1a_chain": chain_id,
                                     "trimmed_rmsd_A": round(fit["sup"].rms, 3)})
    out_sum = core.save_json("structure", f"KIF13A_structure_axis_summary_{stamp}.json", {
        "pdb_id": args.pdb_id, "pdb_title": meta["title"], "pdb_citation": citation,
        "kif1a_chain_primary": chain_id, "tubulin_head_frames": mt_frames,
        "per_frame_trimmed_rmsd_A": {c: round(fits[c]["sup"].rms, 3) for c in fits},
        "kif1a_chains": kif1a_chains, "tubulin_chains": tubulin_chains, "ligand": lig_names,
        "n_pairs_all": len(fit["pairs"]), "rmsd_all_A": round(fit["raw_rmsd"], 3),
        "n_pairs_trimmed": len(fit["kept"]), "rmsd_trimmed_A": round(fit["sup"].rms, 3),
        "kif1a_residues_touching_tubulin": kif1a_contact,
        "kif13a_residues_clashing_with_tubulin_under_2A": clashes,
        "cutoffs": {"atp_A": args.atp_cutoff, "mt_A": args.mt_cutoff, "buried_rel_sasa": args.buried_cutoff,
                    "trim_A": args.trim},
        "motor_class_counts": dict(summary_counts),
        "candidate_class_counts": {t: dict(Counter(r["structural_class"] for r in joined if r["tier"] == t))
                                   for t in "ABCD"},
        "sanity_checks": checks,
    }, base=core.INTERIM, source="derived", tool="scripts/structure_axis.py", urls=[])
    for rec in (out_res, out_cand, out_sum):
        print(f"  wrote {rec['path']} ({rec['bytes']} bytes)")

    for tier in "ABCD":
        rows = [r for r in joined if r["tier"] == tier]
        print(f"\ntier {tier} ({len(rows)} variants, {sum(1 for r in rows if r['confidently_observed'] == 'True')} "
              f"confidently observed): {Counter(r['structural_class'] for r in rows).most_common()}")
    print("\nTIER A with structural context")
    for r in [x for x in joined if x["tier"] == "A"]:
        conf_flag = "" if r["confidently_observed"] == "True" else "  [AC0]"
        if r["dist_nucleotide_A"] is None:
            print(f"  {r['protein_change']:<9} {r['structural_class']}{conf_flag}")
            continue
        print(f"  {r['protein_change']:<9} {r['structural_class']:<20} ATP {r['dist_nucleotide_A']:>6.1f} A  "
              f"MT {r['dist_tubulin_A']:>6.1f} A  relSASA {r['rel_sasa']:.2f}  pLDDT {r['plddt']:.0f}{conf_flag}")

    if not args.no_figure:
        build_figure(residues, joined, ROOT / "figures" / "KIF13A_structure_axis.svg",
                     args.atp_cutoff, args.mt_cutoff)
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
