"""Build tleap/Amber topologies for the Rab10(GTP, Mg2+):KIF13A(1085-1160) complex arms.

    python prepare_complex.py --af3-root AF3_DIR --out-root ../systems --all
    python prepare_complex.py --af3-root AF3_DIR --out-root ../systems --system WT_wt_s5

Needs AmberTools (tleap), OpenMM 8.3.1, PDBFixer 1.12, PROPKA 3.5.1, numpy and scipy in ONE
interpreter (the ambertools.sif image or the local ~/amber-env). Per system (systems/<arm>_<start>/):

 1. AF3 mmCIF (model_0) -> heavy-atom model, verified against the UniProt sequences, UniProt numbering
    kept (Rab10 P61026 1-175, KIF13A Q9H1H9 = MANE 1085-1160), OXT dropped.
 2. Paired arms: residue 1122 side chain stripped beyond CB and rebuilt by PDBFixer (K1122E, K1122R).
 3. Caps: Rab10 M1 free (NH3+), Rab10 P175-NME, KIF13A ACE-S1085 ... P1160-NME. The monomer build.py
    copies cap atoms from the neighbouring AF residue; the AF3 construct has no neighbours, so the
    neighbour is generated with ideal peptide geometry: NME N on the AF3 OXT direction (in the carboxyl
    plane) and CH3 trans (omega 180); ACE C/O/CH3 from an ideal residue 1084 with the clash-free phi.
 4. PROPKA 3 on the heavy-atom complex -> Rab10 H72/H110/H144 states (rule in complex_protocol.json),
    neutral tautomers chosen by heavy-atom H-bond partners.
 5. GTP: CCD atom names == Meagher/Redman/Carlson prep names (explicit map below, 32 heavy atoms);
    tleap adds the 12 hydrogens from the prep template.
 6. tleap (ff19SB, OPC, OPC ion set, GTP.prep + frcmod.phos) -> dry topology (mbondi3) ->
    restrained GBn2 minimisation (hydrogens and rebuilt side chain free) -> solvateOct OPC with a
    rotation-safe buffer -> Na+/Cl- (0.15 M + neutralising) -> complex.prmtop/.inpcrd.
 7. Mg2+ coordination shell from the AF3 model (<= 2.6 A) -> mg_restraints.json (prmtop indices).
 8. QC: charges, counts, OPC virtual sites, box, periodic water overlaps, Mg distances, energies,
    sha256 of every input -> qc.json; residue_map.tsv/json keyed by UniProt numbers.
"""

from __future__ import annotations

import argparse
import io
import itertools
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from complex_common import (DEFAULT_PROTOCOL, HERE, load_protocol, provenance, read_json, sha256,
                            system_by_id, system_manifest, utc_now, write_json, write_text)

ROOT = HERE.parent  # complex_md/
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
                "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
                "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
                "HID": "H", "HIE": "H", "HIP": "H"}
# CCD GTP (AF3 output) -> Meagher et al. GTP.prep heavy-atom names. Identical by design of the prep
# file; kept explicit so a renamed source is caught. Gamma (O1G/O2G/O3G), beta (O1B/O2B) and alpha
# (O1A/O2A) non-bridging oxygens are chemically equivalent within each set in the prep (same type and
# charge), so the CCD stereo labels need no swapping.
GTP_CCD_TO_PREP = {name: name for name in [
    "PG", "O1G", "O2G", "O3G", "O3B", "PB", "O1B", "O2B", "O3A", "PA", "O1A", "O2A", "O5'", "C5'",
    "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2",
    "N2", "N3", "C4"]}
CAPS = {"ACE": ["CH3", "C", "O"], "NME": ["N", "CH3"]}
POS_ATOMS = {"ARG": ["NE", "NH1", "NH2"], "LYS": ["NZ"], "HIP": ["ND1", "NE2"]}
NEG_ATOMS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}


# ----------------------------------------------------------------------------- geometry

def dihedral(a, b, c, d) -> float:
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def place(a, b, c, bond: float, angle_deg: float, torsion_deg: float) -> np.ndarray:
    """NeRF: atom d with |cd| = bond, angle(b, c, d) = angle, dihedral(a, b, c, d) = torsion."""
    angle, torsion = math.radians(angle_deg), math.radians(torsion_deg)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * math.cos(angle), bond * math.sin(angle) * math.cos(torsion),
                   bond * math.sin(angle) * math.sin(torsion)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def angle(a, b, c) -> float:
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


# ----------------------------------------------------------------------------- model container

class Residue:
    def __init__(self, molecule: str, chain: str, name: str, number: int | None, label: str):
        self.molecule = molecule      # Rab10, KIF13A, GTP, MG
        self.chain = chain            # PDB chain written for tleap/PROPKA
        self.name = name              # residue name (HIS variants set later)
        self.number = number          # UniProt/MANE number, None for caps/ligands
        self.label = label            # e.g. "Rab10 K22", "KIF13A ACE"
        self.atoms: dict[str, tuple[str, np.ndarray]] = {}

    def xyz(self, name: str) -> np.ndarray:
        return self.atoms[name][1]

    def heavy(self):
        return [(n, e, x) for n, (e, x) in self.atoms.items() if e != "H"]


def read_af3(cif: Path, protocol: dict) -> tuple[list[Residue], dict]:
    from openmm import app, unit
    pdbx = app.PDBxFile(str(cif))
    xyz = np.asarray(pdbx.positions.value_in_unit(unit.angstrom))
    con = protocol["construct"]
    chains = {c.id: c for c in pdbx.topology.chains()}
    expected_chains = [con["rab10"]["af3_chain"], con["kif13a"]["af3_chain"], con["gtp"]["af3_chain"], con["mg"]["af3_chain"]]
    if sorted(chains) != sorted(expected_chains):
        raise ValueError(f"AF3 chains {sorted(chains)} != {expected_chains}")
    residues: list[Residue] = []
    report: dict[str, Any] = {}
    for key, offset in [("rab10", con["rab10"]["first"] - 1), ("kif13a", con["kif13a"]["first"] - 1)]:
        spec = con[key]
        chain = chains[spec["af3_chain"]]
        seq = []
        for r in chain.residues():
            res = Residue(spec["label"], spec["pdb_chain"], r.name, int(r.id) + offset, f"{spec['label']} {r.name}{int(r.id) + offset}")
            for a in r.atoms():
                if a.element is None:
                    raise ValueError("Atom without element in AF3 model")
                if a.element.symbol == "H":
                    continue
                if a.name == "OXT":
                    res.atoms["OXT"] = ("O", xyz[a.index])
                    continue
                res.atoms[a.name] = (a.element.symbol, xyz[a.index])
            residues.append(res)
            seq.append(THREE_TO_ONE[r.name])
        numbers = [r.number for r in residues if r.molecule == spec["label"]]
        if numbers != list(range(spec["first"], spec["last"] + 1)):
            raise ValueError(f"{spec['label']} numbering mismatch")
        report[f"{key}_sequence"] = "".join(seq)
    gtp_chain = list(chains[con["gtp"]["af3_chain"]].residues())
    mg_chain = list(chains[con["mg"]["af3_chain"]].residues())
    if len(gtp_chain) != 1 or gtp_chain[0].name != "GTP" or len(mg_chain) != 1 or mg_chain[0].name != "MG":
        raise ValueError("AF3 ligand chains must hold one GTP and one MG")
    gtp = Residue("GTP", con["gtp"]["pdb_chain"], "GTP", None, "GTP")
    for a in gtp_chain[0].atoms():
        if a.element.symbol != "H":
            gtp.atoms[a.name] = (a.element.symbol, xyz[a.index])
    if set(gtp.atoms) != set(GTP_CCD_TO_PREP) or len(gtp.atoms) != con["gtp"]["heavy_atoms"]:
        raise ValueError(f"AF3 GTP heavy atoms differ from CCD naming: {sorted(gtp.atoms)}")
    gtp.atoms = {GTP_CCD_TO_PREP[n]: v for n, v in gtp.atoms.items()}
    mg = Residue("MG", con["mg"]["pdb_chain"], "MG", None, "MG")
    for a in mg_chain[0].atoms():
        mg.atoms["MG"] = ("Mg", xyz[a.index])
    residues += [gtp, mg]
    return residues, report


def find(residues: list[Residue], molecule: str, number: int) -> Residue:
    hits = [r for r in residues if r.molecule == molecule and r.number == number and r.name not in CAPS]
    if len(hits) != 1:
        raise KeyError(f"{molecule} {number}")
    return hits[0]


# ----------------------------------------------------------------------------- PDB writer

def atom_field(name: str, element: str) -> str:
    if len(name) >= 4 or len(element) == 2:
        return f"{name:<4s}"
    return f" {name:<3s}"


def pdb_text(residues: list[Residue], *, his_names: bool = True, het_ligands: bool = True, leap: bool = False) -> str:
    """leap=True: Amber library names (aminoct12.lib NME methyl carbon is "C", not "CH3")."""
    lines, serial = [], 1
    previous_chain = None
    for i, res in enumerate(residues):
        if previous_chain is not None and res.chain != previous_chain:
            lines.append("TER")
        previous_chain = res.chain
        resname = res.name if his_names or res.name not in {"HID", "HIE", "HIP"} else "HIS"
        record = "HETATM" if (het_ligands and res.molecule in {"GTP", "MG"}) else "ATOM  "
        number = res.number if res.number is not None else (1 if res.molecule in {"GTP", "MG"} else res.cap_number)
        for name, (element, xyz) in res.atoms.items():
            if leap and res.name == "NME" and name == "CH3":
                name = "C"
            lines.append(f"{record}{serial:5d} {atom_field(name, element)} {resname:>3s} {res.chain}{number:4d}    "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}")
            serial += 1
    lines += ["TER", "END"]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------- caps

def add_caps(residues: list[Residue], protocol: dict) -> tuple[list[Residue], list[dict]]:
    """ACE before KIF13A 1085, NME after Rab10 175 and KIF13A 1160 (see module docstring)."""
    checks = []
    rab_last = find(residues, "Rab10", protocol["construct"]["rab10"]["last"])
    kif_first = find(residues, "KIF13A", protocol["construct"]["kif13a"]["first"])
    kif_last = find(residues, "KIF13A", protocol["construct"]["kif13a"]["last"])

    def nme_after(res: Residue, chain: str, cap_number: int) -> Residue:
        n_ca, ca, c = res.xyz("N"), res.xyz("CA"), res.xyz("C")
        oxt = res.atoms.pop("OXT")[1]
        direction = (oxt - c) / np.linalg.norm(oxt - c)
        n = c + 1.335 * direction
        ch3 = place(ca, c, n, 1.449, 121.9, 180.0)
        cap = Residue(res.molecule, chain, "NME", None, f"{res.molecule} NME")
        cap.cap_number = cap_number
        cap.atoms = {"N": ("N", n), "CH3": ("C", ch3)}
        checks.append({"cap": f"{res.molecule} NME after {res.number}", "method": "N on AF3 OXT direction, omega 180",
                       "C-N_A": float(np.linalg.norm(n - c)), "CA-C-N_deg": angle(ca, c, n),
                       "O-C-N_deg": angle(res.xyz("O"), c, n),
                       "omega_deg": dihedral(ca, c, n, ch3), "psi_deg": dihedral(n_ca, ca, c, n)})
        return cap

    def ace_before(res: Residue, chain: str, cap_number: int) -> Residue:
        n, ca, c = res.xyz("N"), res.xyz("CA"), res.xyz("C")
        others = np.array([x for r in residues if r is not res for _, (e, x) in r.atoms.items() if e != "H"])
        trials = []
        for phi in [-65.0, -90.0, -120.0, -150.0, 180.0, 60.0, -30.0, 120.0]:  # preference order
            c_ace = place(c, ca, n, 1.335, 121.9, phi)
            ch3 = place(ca, n, c_ace, 1.522, 116.2, 180.0)
            o = place(ca, n, c_ace, 1.229, 122.7, 0.0)
            new = np.array([c_ace, ch3, o])
            clearance = float(np.min(np.linalg.norm(new[:, None] - others[None], axis=-1)))
            trials.append((phi, clearance, c_ace, ch3, o))
            if clearance >= 3.0:
                break
        best = trials[-1] if trials[-1][1] >= 3.0 else max(trials, key=lambda t: t[1])
        phi, clearance, c_ace, ch3, o = best
        cap = Residue(res.molecule, chain, "ACE", None, f"{res.molecule} ACE")
        cap.cap_number = cap_number
        cap.atoms = {"CH3": ("C", ch3), "C": ("C", c_ace), "O": ("O", o)}
        checks.append({"cap": f"{res.molecule} ACE before {res.number}", "method": "ideal residue, omega 180",
                       "phi_deg": phi, "clearance_A": clearance, "C-N_A": float(np.linalg.norm(c_ace - n)),
                       "C-N-CA_deg": angle(c_ace, n, ca), "omega_deg": dihedral(ch3, c_ace, n, ca)})
        return cap

    out = []
    for res in residues:
        if res is kif_first:
            out.append(ace_before(res, res.chain, res.number - 1))
        out.append(res)
        if res is rab_last:
            out.append(nme_after(res, res.chain, res.number + 1))
        if res is kif_last:
            out.append(nme_after(res, res.chain, res.number + 1))
    for check in checks:
        if not 1.30 <= check["C-N_A"] <= 1.37:
            raise ValueError(f"Cap peptide bond out of range: {check}")
        if abs(abs(check.get("omega_deg", 180.0)) - 180.0) > 1.0:
            raise ValueError(f"Cap omega not trans: {check}")
    leftover = [r.label for r in out if "OXT" in r.atoms]
    if leftover:
        raise ValueError(f"OXT left on {leftover}")
    return out, checks


# ----------------------------------------------------------------------------- mutation

def rebuild_beyond_cb(residues: list[Residue], mutation: dict, seed: int) -> dict:
    """Strip residue side chain beyond CB, rename, let PDBFixer rebuild heavy atoms (template + local fit)."""
    from openmm import app, unit
    from pdbfixer import PDBFixer
    target = find(residues, "KIF13A", mutation["position"])
    if target.name != mutation["from"]:
        raise ValueError(f"Residue {mutation['position']} is {target.name}, expected {mutation['from']}")
    protein = [r for r in residues if r.molecule in {"Rab10", "KIF13A"} and r.name not in CAPS]
    text_lines, serial, prev = [], 1, None
    for r in protein:
        if prev is not None and r.chain != prev:
            text_lines.append("TER")
        prev = r.chain
        for name, (element, xyz) in r.atoms.items():
            if name == "OXT":
                continue
            text_lines.append(f"ATOM  {serial:5d} {atom_field(name, element)} {r.name:>3s} {r.chain}{r.number:4d}    "
                              f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}")
            serial += 1
    text_lines += ["TER", "END"]
    fixer = PDBFixer(pdbfile=io.StringIO("\n".join(text_lines) + "\n"))
    modeller = app.Modeller(fixer.topology, fixer.positions)
    hit = [r for r in modeller.topology.residues() if r.chain.id == target.chain and int(r.id) == mutation["position"]]
    if len(hit) != 1:
        raise ValueError("Mutation residue not found in PDBFixer topology")
    removed = [a.name for a in hit[0].atoms() if a.name not in {"N", "CA", "C", "O", "CB"}]
    modeller.delete([a for a in hit[0].atoms() if a.name not in {"N", "CA", "C", "O", "CB"}])
    hit = [r for r in modeller.topology.residues() if r.chain.id == target.chain and int(r.id) == mutation["position"]]
    hit[0].name = mutation["to"]
    before = np.asarray(modeller.positions.value_in_unit(unit.angstrom))
    fixer.topology, fixer.positions = modeller.topology, modeller.positions
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.missingTerminals = {}
    touched = {(r.chain.id, int(r.id)) for r in fixer.missingAtoms}
    if touched != {(target.chain, mutation["position"])}:
        raise ValueError(f"PDBFixer wants to add atoms outside the mutation site: {touched}")
    fixer.addMissingAtoms(seed=seed)
    after = np.asarray(fixer.positions.value_in_unit(unit.angstrom))
    new_atoms, kept_shift = {}, 0.0
    old_names = [(a.residue.chain.id, int(a.residue.id), a.name) for a in modeller.topology.atoms()]
    old_xyz = dict(zip(old_names, before))
    for a in fixer.topology.atoms():
        key = (a.residue.chain.id, int(a.residue.id), a.name)
        if key in old_xyz:
            kept_shift = max(kept_shift, float(np.linalg.norm(after[a.index] - old_xyz[key])))
        elif key[:2] == (target.chain, mutation["position"]) and a.element.symbol != "H":
            new_atoms[a.name] = (a.element.symbol, after[a.index])
    if kept_shift > 1e-3:
        raise ValueError(f"PDBFixer moved existing atoms by {kept_shift:.4f} A")
    original_side_chain = {n: [float(v) for v in target.atoms[n][1]] for n in removed}
    kept = {n: target.atoms[n] for n in ["N", "CA", "C", "O", "CB"]}
    target.atoms = {**kept, **new_atoms}
    target.name = mutation["to"]
    target.label = f"KIF13A {mutation['to']}{mutation['position']}"
    expected = {"GLU": {"CG", "CD", "OE1", "OE2"}, "ARG": {"CG", "CD", "NE", "CZ", "NH1", "NH2"}}[mutation["to"]]
    if set(new_atoms) != expected:
        raise ValueError(f"Rebuilt atoms {sorted(new_atoms)} != {sorted(expected)}")
    return {"position": mutation["position"], "from": mutation["from"], "to": mutation["to"],
            "removed_beyond_cb": removed, "original_side_chain_xyz_A": original_side_chain,
            "rebuilt": sorted(new_atoms), "max_shift_existing_atoms_A": kept_shift,
            "tool": "PDBFixer.addMissingAtoms", "seed": seed}


# ----------------------------------------------------------------------------- PROPKA and His

def run_propka(residues: list[Residue], out: Path, python: str) -> dict:
    pdb = out / "propka_input.pdb"
    write_text(pdb, pdb_text(residues, his_names=False))
    for stale in [out / "propka_input.pka"]:
        if stale.exists():
            stale.unlink()
    proc = subprocess.run([python, "-m", "propka", "-q", pdb.name], cwd=out, capture_output=True, text=True)
    write_text(out / "propka_stdout.txt", proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if proc.returncode != 0 or not (out / "propka_input.pka").exists():
        raise RuntimeError(f"PROPKA failed: {proc.stderr[-2000:]}")
    text = (out / "propka_input.pka").read_text()
    summary = text.split("SUMMARY OF THIS PREDICTION")[1].split("------")[0]
    groups = []
    for line in summary.splitlines():
        m = re.match(r"^\s+(\S+)\s+(\d+)\s+(\S)\s+(-?[\d.]+)\s+(-?[\d.]+)", line)
        if m:
            groups.append({"group": m.group(1), "number": int(m.group(2)), "chain": m.group(3),
                           "pka": float(m.group(4)), "model_pka": float(m.group(5))})
    version = subprocess.run([python, "-c", "import propka; print(propka.__version__)"], capture_output=True, text=True).stdout.strip()
    ligand_lines = [l for l in proc.stderr.splitlines() + proc.stdout.splitlines() if "GTP" in l or "ligand" in l.lower()]
    determinants = {}
    body = text.split("SUMMARY OF THIS PREDICTION")[0]
    for line in body.splitlines():
        m = re.match(r"^(HIS)\s*(\d+) (\S)\s", line)
        if m:
            determinants.setdefault(f"{m.group(3)}:{int(m.group(2))}", []).append(" ".join(line.split()))
    return {"version": version, "input": pdb.name, "input_sha256": sha256(pdb), "groups": groups,
            "his_determinant_lines": determinants,
            "ligand_note": "PROPKA 3.5.1 reads the GTP HETATM record (GTP oxygens appear as determinants, e.g. of R70); Mg2+ is not used",
            "notes": ligand_lines[:20]}


def his_states(residues: list[Residue], propka: dict, protocol: dict) -> list[dict]:
    rule = protocol["his"]
    cutoff = rule["hbond_cutoff_A"]
    acceptors, donors, both = [], [], []
    for r in residues:
        for name, (element, xyz) in r.atoms.items():
            entry = (r, name, xyz)
            if r.molecule == "GTP":
                if element == "O":
                    (both if name in {"O2'", "O3'"} else acceptors).append(entry)
                elif name in {"N1", "N2"}:
                    donors.append(entry)
                elif name in {"N3", "N7"}:
                    acceptors.append(entry)
                continue
            if r.molecule == "MG":
                continue
            if name == "O" or (r.name, name) in {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
                                                   ("ASN", "OD1"), ("GLN", "OE1")}:
                acceptors.append(entry)
            elif (name == "N" and r.name not in {"PRO", "ACE"}) or (r.name, name) in {
                    ("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("ASN", "ND2"),
                    ("GLN", "NE2"), ("TRP", "NE1")}:
                donors.append(entry)
            elif (r.name, name) in {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH")}:
                both.append(entry)
    pkas = {(g["chain"], g["number"]): g["pka"] for g in propka["groups"] if g["group"] == "HIS"}
    states = []
    for number in rule["rab10_residues"]:
        res = find(residues, "Rab10", number)
        if res.name != "HIS":
            raise ValueError(f"Rab10 {number} is {res.name}, not HIS")
        pka = pkas[(res.chain, number)]

        def partners(atom: str) -> dict:
            xyz = res.xyz(atom)
            near = {"acceptor": [], "donor": [], "both": []}
            for kind, pool in [("acceptor", acceptors), ("donor", donors), ("both", both)]:
                for r, name, x in pool:
                    if r is res:
                        continue
                    d = float(np.linalg.norm(x - xyz))
                    if d <= cutoff:
                        near[kind].append(f"{r.label}:{name} {d:.2f}")
            return near

        nd1, ne2 = partners("ND1"), partners("NE2")
        score_hid = len(nd1["acceptor"]) + len(ne2["donor"]) + 0.5 * (len(nd1["both"]) + len(ne2["both"]))
        score_hie = len(ne2["acceptor"]) + len(nd1["donor"]) + 0.5 * (len(ne2["both"]) + len(nd1["both"]))
        tautomer = "HID" if score_hid > score_hie else "HIE"
        if pka > rule["protonated_above_pka"]:
            state, flag = "HIP", None
        elif pka < rule["neutral_below_pka"]:
            state, flag = tautomer, None
        else:
            state = "HIP" if pka >= rule["ph"] else tautomer
            flag = "candidate_sensitivity_arm"
        res.name = state
        states.append({"residue": f"Rab10 H{number}", "propka_pka": pka, "state": state, "flag": flag,
                       "propka_determinants": propka["his_determinant_lines"].get(f"{res.chain}:{number}", []),
                       "neutral_tautomer_by_hbond": tautomer, "score_HID": score_hid, "score_HIE": score_hie,
                       "ND1_partners": nd1, "NE2_partners": ne2})
    return states


# ----------------------------------------------------------------------------- tleap

def leap_header(protocol: dict, params_dir: Path) -> str:
    p = protocol["parameters"]
    lines = [f"source {name}" for name in p["leaprc"]]
    lines.append("addAtomTypes {")
    lines += [f'  {{ "{t}" "{e}" "{h}" }}' for t, e, h in p["extra_atom_types"]]
    lines.append("}")
    lines.append(f"loadAmberPrep {(params_dir / Path(p['gtp_prep']).name).as_posix()}")
    lines.append(f"loadAmberParams {(params_dir / Path(p['gtp_frcmod']).name).as_posix()}")
    lines.append(f"GTP = {p['gtp_prep_unit']}")
    lines.append(f"set default PBRadii {p['pb_radii']}")
    return "\n".join(lines) + "\n"


def run_leap(tleap: str, out: Path, name: str, body: str, header: str) -> str:
    script = out / f"leap_{name}.in"
    write_text(script, header + body + "quit\n")
    proc = subprocess.run([tleap, "-f", script.name], cwd=out, capture_output=True, text=True)
    log = proc.stdout + ("\n--- stderr ---\n" + proc.stderr if proc.stderr else "")
    write_text(out / f"leap_{name}.log", log)
    if (out / "leap.log").exists():
        (out / "leap.log").unlink()
    m = re.search(r"Exiting LEaP: Errors = (\d+); Warnings = (\d+); Notes = (\d+)", log)
    if proc.returncode != 0 or not m or int(m.group(1)) != 0:
        raise RuntimeError(f"tleap {name} failed, see {out / f'leap_{name}.log'}")
    return log


def leap_charges(log: str) -> list[float]:
    return [float(x) for x in re.findall(r"Total (?:unperturbed )?charge:\s+(-?[\d.]+)", log)]


def read_prmtop(path: Path) -> dict[str, list]:
    flags, current, fmt = {}, None, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("%FLAG"):
            current = line.split()[1]
            flags[current] = []
        elif line.startswith("%FORMAT"):
            fmt = re.match(r"%FORMAT\((\d+)([aIE])(\d+)", line.replace("%FORMAT(", "%FORMAT(").strip())
            fmt = (int(fmt.group(1)), fmt.group(2), int(fmt.group(3))) if fmt else None
        elif current and fmt and not line.startswith("%"):
            count, kind, width = fmt
            for i in range(0, len(line), width):
                token = line[i:i + width]
                if not token.strip() and kind != "a":
                    continue
                flags[current].append(token.strip() if kind == "a" else int(token) if kind == "I" else float(token))
    return flags


def read_inpcrd_box(path: Path) -> list[float]:
    last = Path(path).read_text().rstrip("\n").splitlines()[-1]
    return [float(last[i:i + 12]) for i in range(0, 72, 12)]


BOUNDARY_OO_CUTOFF_A = 2.3


def boundary_overlap_waters(prmtop_path: Path, inpcrd_path: Path, cutoff_A: float = BOUNDARY_OO_CUTOFF_A) -> list[dict]:
    """Waters whose O lies within cutoff of another water O across the truncated-octahedron faces.

    tleap clips the replicated OPC box to the octahedron; a few waters on opposite faces then overlap
    through the periodic boundary (seen here: ~19 pairs, O-O down to 1.45 A, forces ~4e6 kJ/mol/nm).
    One water of each pair (the higher residue number) is removed in the final tleap session.
    """
    from openmm import app, unit
    from scipy.spatial import cKDTree
    inpcrd = app.AmberInpcrdFile(str(inpcrd_path))
    vec = np.array([[v.x, v.y, v.z] for v in inpcrd.boxVectors.value_in_unit(unit.angstrom)])
    xyz = read_inpcrd_xyz(inpcrd_path)
    top = read_prmtop(prmtop_path)
    labels, pointers = top["RESIDUE_LABEL"], top["RESIDUE_POINTER"]
    waters = [(r + 1, pointers[r] - 1) for r, label in enumerate(labels) if label == "WAT"]
    o_xyz = np.array([xyz[first] for _, first in waters])
    tree = cKDTree(o_xyz)
    remove = {}
    for cell in itertools.product((-1, 0, 1), repeat=3):
        if cell == (0, 0, 0):
            continue
        shifted = o_xyz + np.array(cell, dtype=float) @ vec
        for i, hits in enumerate(tree.query_ball_point(shifted, cutoff_A)):
            for j in hits:
                a, b = waters[i][0], waters[j][0]
                d = float(np.linalg.norm(shifted[i] - o_xyz[j]))
                victim = max(a, b)
                if victim not in remove or d < remove[victim]["distance_A"]:
                    remove[victim] = {"leap_residue": victim, "partner": min(a, b), "distance_A": d, "cell": list(cell)}
    # keep one victim per pair: drop entries whose partner is itself a victim of the same pair
    return sorted(remove.values(), key=lambda r: r["leap_residue"])


def read_inpcrd_xyz(path: Path) -> np.ndarray:
    lines = Path(path).read_text().splitlines()
    natom = int(lines[1].split()[0])
    values = []
    for line in lines[2:]:
        values += [float(line[i:i + 12]) for i in range(0, len(line), 12) if line[i:i + 12].strip()]
        if len(values) >= 3 * natom:
            break
    return np.asarray(values[:3 * natom]).reshape(natom, 3)


def replace_pdb_coords(template: Path, xyz_angstrom: np.ndarray, target: Path) -> None:
    out, k = [], 0
    for line in template.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x, y, z = xyz_angstrom[k]
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            k += 1
        out.append(line)
    if k != len(xyz_angstrom):
        raise ValueError(f"PDB template has {k} atoms, coordinates {len(xyz_angstrom)}")
    write_text(target, "\n".join(out) + "\n")


def split_pdb_by_residue_names(source: Path, keep, target: Path) -> None:
    """Write the lines of selected leap residues (keep(resname, resseq) -> bool), keeping TER cards."""
    out, last_kept = [], False
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            last_kept = keep(line[17:21].strip(), int(line[22:26]))
            if last_kept:
                out.append(line)
        elif line.startswith("TER"):
            if out and out[-1] != "TER" and last_kept:
                out.append("TER")
    out.append("END")
    write_text(target, "\n".join(out) + "\n")


# ----------------------------------------------------------------------------- residue map

def residue_map(prmtop: dict, residues: list[Residue], protocol: dict) -> list[dict]:
    labels = prmtop["RESIDUE_LABEL"]
    pointers = prmtop["RESIDUE_POINTER"] + [prmtop["POINTERS"][0] + 1]
    solute = [r for r in residues]
    rows = []
    for i, res in enumerate(solute):
        label = labels[i]
        expected = {"GTP": "gtp", "MG": "MG"}.get(res.molecule, res.name)
        if label != expected:
            raise ValueError(f"prmtop residue {i} {label} != expected {expected} ({res.label})")
        accession = {"Rab10": protocol["construct"]["rab10"]["uniprot"],
                     "KIF13A": protocol["construct"]["kif13a"]["uniprot"]}.get(res.molecule)
        rows.append({"prmtop_residue_index": i, "leap_residue_number": i + 1, "molecule": res.molecule,
                     "pdb_chain": res.chain, "uniprot": accession if res.name not in CAPS else None,
                     "uniprot_number": res.number if res.name not in CAPS else None,
                     "resname_prmtop": label, "one_letter": THREE_TO_ONE.get(res.name, ""),
                     "cap": res.name in CAPS, "first_atom": pointers[i] - 1, "end_atom": pointers[i + 1] - 1})
    return rows


# ----------------------------------------------------------------------------- relaxation

def tethered_relaxation(out: Path, protocol: dict, free_residue_index: int, mg_pairs: list[dict], tether: dict) -> dict:
    """Restrained GBn2 relaxation with temporary flat-bottom pulls that seat a rebuilt R1122 like WT K1122.

    Same system, positional restraints and minimiser settings as restrained_relaxation. Stage "pull_position":
    centroid of R1122 NE/NH1/NH2 pulled to within flat_bottom_A of the WT K1122 NZ position of the same AF3 start.
    Stage "pull_salt_bridge" (only if the closest R1122 N to E1130 OE1/OE2 is still > accept_A): adds a pull of
    that closest N-O distance to <= target_A. Stage "release": every pull switched off, re-minimised, so the
    reported geometry is a pull-free local minimum of the side chain.
    """
    import openmm as mm
    from openmm import app, unit
    settings = protocol["mutant_relax"]
    spec = tether["spec"]
    prmtop = app.AmberPrmtopFile(str(out / "dry_prerelax.prmtop"))
    inpcrd = app.AmberInpcrdFile(str(out / "dry_prerelax.inpcrd"))
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=app.GBn2,
                                 implicitSolventSaltConc=protocol["solvation"]["salt_molar"] * unit.molar,
                                 temperature=protocol["simulation"]["temperature_k"] * unit.kelvin,
                                 soluteDielectric=1.0, solventDielectric=78.5)
    atoms = list(prmtop.topology.atoms())
    xyz0 = np.asarray(inpcrd.positions.value_in_unit(unit.nanometer))
    free = set()
    for a in atoms:
        if a.element is not None and a.element.symbol == "H":
            free.add(a.index)
        elif a.residue.index == free_residue_index and a.name not in {"N", "CA", "C", "O", "CB"}:
            free.add(a.index)
    restraint = mm.CustomExternalForce("0.5*kpos*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("kpos", settings["restraint_k_kj_mol_nm2"])
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    for a in atoms:
        if a.index not in free:
            restraint.addParticle(a.index, xyz0[a.index])
    restraint.setForceGroup(30)
    system.addForce(restraint)
    system.addForce(mg_force(mg_pairs, protocol["mg_restraints"]["k_kj_mol_nm2"], periodic=False))
    n_idx, o_idx = tether["guanidinium_n"], tether["e1130_o"]
    if len(n_idx) != 3 or len(o_idx) != 2:
        raise ValueError("Tether needs R1122 NE/NH1/NH2 and E1130 OE1/OE2")
    pos = mm.CustomCompoundBondForce(3, "0.5*kpull_pos*step(d-tol)*(d-tol)^2;"
                                        "d=sqrt((xc-tx)^2+(yc-ty)^2+(zc-tz)^2);"
                                        "xc=(x1+x2+x3)/3; yc=(y1+y2+y3)/3; zc=(z1+z2+z3)/3")
    pos.addGlobalParameter("kpull_pos", 0.0)
    for name in ("tx", "ty", "tz", "tol"):
        pos.addPerBondParameter(name)
    target = np.asarray(tether["target_nm"], dtype=float)
    pos.addBond(n_idx, [*target, spec["position_pull"]["flat_bottom_A"] / 10])
    pos.setForceGroup(28)
    system.addForce(pos)
    terms = [f"distance(p{i + 1},p{j + 4})" for i in range(3) for j in range(2)]
    dmin = terms[0]
    for term in terms[1:]:
        dmin = f"min({dmin},{term})"
    sb = mm.CustomCompoundBondForce(5, f"0.5*kpull_sb*step(d-r0)*(d-r0)^2; d={dmin}")
    sb.addGlobalParameter("kpull_sb", 0.0)
    sb.addPerBondParameter("r0")
    sb.addBond(list(n_idx) + list(o_idx), [spec["salt_bridge_pull"]["target_A"] / 10])
    sb.setForceGroup(27)
    system.addForce(sb)
    context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("CPU"))
    context.setPositions(inpcrd.positions)
    phys = set(range(27))

    def measure(stage: str) -> dict:
        state = context.getState(getEnergy=True, getPositions=True, groups=phys)
        x = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        centroid = x[n_idx].mean(axis=0)
        no = min(float(np.linalg.norm(x[i] - x[j])) for i in n_idx for j in o_idx)
        pull = sum(context.getState(getEnergy=True, groups={g}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                   for g in (27, 28))
        return {"stage": stage, "potential_kj_mol": state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                "centroid_to_wt_nz_A": float(np.linalg.norm(centroid - target)) * 10,
                "closest_n_to_e1130_o_A": no * 10, "pull_energy_kj_mol": pull}

    t0 = time.time()
    stages = [measure("start")]
    context.setParameter("kpull_pos", spec["position_pull"]["k_kj_mol_nm2"])
    mm.LocalEnergyMinimizer.minimize(context, settings["tolerance_kj_mol_nm"], settings["max_iterations"])
    stages.append(measure("pull_position"))
    variant = "position"
    if stages[-1]["closest_n_to_e1130_o_A"] > spec["salt_bridge_pull"]["accept_A"]:
        context.setParameter("kpull_sb", spec["salt_bridge_pull"]["k_kj_mol_nm2"])
        mm.LocalEnergyMinimizer.minimize(context, settings["tolerance_kj_mol_nm"], settings["max_iterations"])
        stages.append(measure("pull_salt_bridge"))
        variant = "position+salt_bridge"
    context.setParameter("kpull_pos", 0.0)
    context.setParameter("kpull_sb", 0.0)
    mm.LocalEnergyMinimizer.minimize(context, settings["tolerance_kj_mol_nm"], settings["max_iterations"])
    stages.append(measure("release"))
    state = context.getState(getEnergy=True, getPositions=True, groups=phys)
    e1 = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    xyz1 = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    if not np.isfinite(e1):
        raise ValueError("Non-finite energy after tethered relaxation")
    heavy = [a.index for a in atoms if a.element.symbol != "H"]
    restrained_heavy = [i for i in heavy if i not in free]
    disp = np.linalg.norm(xyz1 - xyz0, axis=1) * 10
    return {"xyz_angstrom": xyz1 * 10, "atoms": atoms,
            "summary": {"method": settings["method"] + "; plus temporary flat-bottom pulls (tether), then released",
                        "k_kj_mol_nm2": settings["restraint_k_kj_mol_nm2"],
                        "free_atoms": len(free), "free_heavy_atoms": len([i for i in free if i in heavy]),
                        "potential_before_kj_mol": stages[0]["potential_kj_mol"], "potential_after_kj_mol": e1,
                        "restrained_heavy_rmsd_A": float(np.sqrt(np.mean(disp[restrained_heavy] ** 2))),
                        "restrained_heavy_max_shift_A": float(disp[restrained_heavy].max()),
                        "free_heavy_max_shift_A": float(max([disp[i] for i in free if i in heavy], default=0.0)),
                        "wall_s": time.time() - t0, "openmm": mm.__version__,
                        "energy_note": "GBn2 potential of the dry complex, restraint and pull terms excluded; diagnostic only"},
            "tether": {"variant": variant, "stages": stages, "target_wt_k1122_nz_A": (target * 10).tolist(),
                       "spec": spec}}


def restrained_relaxation(out: Path, protocol: dict, free_residue_index: int | None, mg_pairs: list[dict]) -> dict:
    import openmm as mm
    from openmm import app, unit
    settings = protocol["mutant_relax"]
    prmtop = app.AmberPrmtopFile(str(out / "dry_prerelax.prmtop"))
    inpcrd = app.AmberInpcrdFile(str(out / "dry_prerelax.inpcrd"))
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=app.GBn2,
                                 implicitSolventSaltConc=protocol["solvation"]["salt_molar"] * unit.molar,
                                 temperature=protocol["simulation"]["temperature_k"] * unit.kelvin,
                                 soluteDielectric=1.0, solventDielectric=78.5)
    atoms = list(prmtop.topology.atoms())
    xyz0 = np.asarray(inpcrd.positions.value_in_unit(unit.nanometer))
    free = set()
    for a in atoms:
        if a.element is not None and a.element.symbol == "H":
            free.add(a.index)
        elif free_residue_index is not None and a.residue.index == free_residue_index and a.name not in {"N", "CA", "C", "O", "CB"}:
            free.add(a.index)
    restraint = mm.CustomExternalForce("0.5*kpos*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("kpos", settings["restraint_k_kj_mol_nm2"])
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    for a in atoms:
        if a.index not in free:
            restraint.addParticle(a.index, xyz0[a.index])
    restraint.setForceGroup(30)
    system.addForce(restraint)
    mgf = mg_force(mg_pairs, protocol["mg_restraints"]["k_kj_mol_nm2"], periodic=False)
    system.addForce(mgf)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("CPU"))
    context.setPositions(inpcrd.positions)
    groups_phys = set(range(29))
    e0 = context.getState(getEnergy=True, groups=groups_phys).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    t0 = time.time()
    mm.LocalEnergyMinimizer.minimize(context, settings["tolerance_kj_mol_nm"], settings["max_iterations"])
    state = context.getState(getEnergy=True, getPositions=True, groups=groups_phys)
    e1 = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    xyz1 = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    heavy = [a.index for a in atoms if a.element.symbol != "H"]
    restrained_heavy = [i for i in heavy if i not in free]
    disp = np.linalg.norm(xyz1 - xyz0, axis=1) * 10
    if not np.isfinite(e1):
        raise ValueError("Non-finite energy after restrained relaxation")
    return {"xyz_angstrom": xyz1 * 10, "atoms": atoms,
            "summary": {"method": settings["method"], "k_kj_mol_nm2": settings["restraint_k_kj_mol_nm2"],
                        "free_atoms": len(free), "free_heavy_atoms": len([i for i in free if i in heavy]),
                        "potential_before_kj_mol": e0, "potential_after_kj_mol": e1,
                        "restrained_heavy_rmsd_A": float(np.sqrt(np.mean(disp[restrained_heavy] ** 2))),
                        "restrained_heavy_max_shift_A": float(disp[restrained_heavy].max()),
                        "free_heavy_max_shift_A": float(max([disp[i] for i in free if i in heavy], default=0.0)),
                        "wall_s": time.time() - t0, "openmm": mm.__version__,
                        "energy_note": "GBn2 potential of the dry complex, restraint terms excluded; diagnostic only"}}


def mg_force(pairs: list[dict], k: float, periodic: bool):
    import openmm as mm
    force = mm.CustomBondForce("step(r-rfb)*0.5*kmg*(r-rfb)^2")
    force.addGlobalParameter("kmg", k)
    force.addPerBondParameter("rfb")
    for pair in pairs:
        force.addBond(pair["mg_index"], pair["ligand_index"], [pair["flat_bottom_nm"]])
    force.setUsesPeriodicBoundaryConditions(periodic)
    force.setForceGroup(29)
    return force


def clash_report(xyz: np.ndarray, atoms, focus_residue_index: int | None, bonds: set) -> dict:
    """Non-bonded heavy-heavy contacts < 2.5 A (Mg ligands excluded) across the solute and around 1122."""
    from scipy.spatial import cKDTree
    heavy = np.array([a.index for a in atoms if a.element.symbol != "H"])
    tree = cKDTree(xyz[heavy])
    pairs = tree.query_pairs(2.5)
    clashes, focus = [], []
    for i, j in pairs:
        a, b = atoms[heavy[i]], atoms[heavy[j]]
        if a.residue.index == b.residue.index or (a.index, b.index) in bonds or (b.index, a.index) in bonds:
            continue
        if "MG" in {a.residue.name, b.residue.name}:
            continue
        # skip 1-3 pairs across peptide bonds
        if abs(a.residue.index - b.residue.index) == 1 and {a.name, b.name} & {"C", "N", "CA", "O", "CH3"}:
            continue
        d = float(np.linalg.norm(xyz[a.index] - xyz[b.index]))
        item = f"{a.residue.name}{a.residue.index + 1}:{a.name}-{b.residue.name}{b.residue.index + 1}:{b.name} {d:.2f}"
        clashes.append(item)
        if focus_residue_index is not None and focus_residue_index in {a.residue.index, b.residue.index}:
            focus.append(item)
    out = {"heavy_pairs_below_2.5A": clashes, "count": len(clashes)}
    if focus_residue_index is not None:
        own = [a.index for a in atoms if a.residue.index == focus_residue_index and a.element.symbol != "H"
               and a.name not in {"N", "CA", "C", "O", "CB"}]
        others = np.array([a.index for a in atoms if a.residue.index != focus_residue_index])
        d = np.linalg.norm(xyz[own][:, None] - xyz[others][None], axis=-1)
        k = np.unravel_index(np.argmin(d), d.shape)
        nearest = atoms[others[k[1]]]
        out["residue_1122"] = {"clashes": focus, "min_distance_to_other_residue_A": float(d.min()),
                               "nearest": f"{nearest.residue.name}{nearest.residue.index + 1}:{nearest.name}"}
    return out


# ----------------------------------------------------------------------------- QC of the solvated system

def solvated_qc(out: Path, protocol: dict, mg_pairs: list[dict], n_solute_atoms: int) -> dict:
    import openmm as mm
    from openmm import app, unit
    from scipy.spatial import cKDTree
    prmtop = app.AmberPrmtopFile(str(out / "complex.prmtop"))
    inpcrd = app.AmberInpcrdFile(str(out / "complex.inpcrd"))
    sim = protocol["simulation"]
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=sim["nonbonded_cutoff_nm"] * unit.nanometer,
                                 constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=sim["ewald_error_tolerance"])
    atoms = list(prmtop.topology.atoms())
    n = system.getNumParticles()
    waters = [r for r in prmtop.topology.residues() if r.name in {"WAT", "HOH"}]
    vsites = [i for i in range(n) if system.isVirtualSite(i)]
    ep_mass = {float(system.getParticleMass(i).value_in_unit(unit.dalton)) for i in vsites}
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    charges = np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(n)])
    opc_ok = True
    for w in waters[:50] + waters[-50:]:
        names = {a.name: a.index for a in w.atoms()}
        if set(names) != {"O", "H1", "H2", "EPW"} or not system.isVirtualSite(names["EPW"]):
            opc_ok = False
        if abs(charges[names["O"]]) > 1e-6 or abs(charges[names["EPW"]] + 1.358284) > 1e-3:
            opc_ok = False
    box = inpcrd.boxVectors
    vec = np.array([[v.x, v.y, v.z] for v in box.value_in_unit(unit.nanometer)]) if box is not None else None
    xyz = np.asarray(inpcrd.positions.value_in_unit(unit.nanometer))
    # periodic overlap / void check on water oxygens (all 27 images; catches a box/lattice mismatch)
    o_idx = np.array([next(a.index for a in w.atoms() if a.name == "O") for w in waters])
    shifts = np.array([i * vec[0] + j * vec[1] + k * vec[2] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    images = (xyz[o_idx][None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    tree = cKDTree(images)
    d, _ = tree.query(xyz[o_idx], k=2)
    nearest_oo = d[:, 1]
    solute_heavy = [a.index for a in atoms[:n_solute_atoms] if a.element is not None and a.element.symbol != "H"]
    solute_tree = cKDTree((xyz[solute_heavy][None] + shifts[:, None, :]).reshape(-1, 3))
    ds, _ = solute_tree.query(xyz[o_idx], k=1)
    ion_names = {"Na+": 0, "Cl-": 0}
    for r in prmtop.topology.residues():
        if r.name in ion_names:
            ion_names[r.name] += 1
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("CPU"))
    context.setPositions(inpcrd.positions)
    if box is not None:
        context.setPeriodicBoxVectors(*box)
    state = context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = np.linalg.norm(np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)), axis=1)
    worst = int(np.argmax(forces))
    mg_now = [{**p, "distance_in_complex_inpcrd_A": float(np.linalg.norm(xyz[p["mg_index"]] - xyz[p["ligand_index"]]) * 10)}
              for p in mg_pairs]
    mg_index = mg_pairs[0]["mg_index"]
    mg_tree_d, _ = cKDTree(images).query(xyz[mg_index], k=8)
    return {
        "particles_openmm": n, "atoms_prmtop": len(atoms), "particles_match": n == len(atoms),
        "water_molecules": len(waters), "virtual_sites": len(vsites), "virtual_sites_equal_waters": len(vsites) == len(waters),
        "virtual_site_masses_dalton": sorted(ep_mass), "opc_spot_check_ok": opc_ok,
        "net_charge_e": float(charges.sum()), "solute_charge_e": float(charges[:n_solute_atoms].sum()),
        "ions": ion_names,
        "box_vectors_nm": vec.tolist() if vec is not None else None,
        "box_volume_nm3": float(abs(np.linalg.det(vec))) if vec is not None else None,
        "water_oo_nearest_min_A": float(nearest_oo.min() * 10), "water_oo_pairs_below_2.2A": int((nearest_oo < 0.22).sum()),
        "water_without_O_neighbour_within_3.5A": int((nearest_oo > 0.35).sum()),
        "water_o_to_solute_heavy_min_A": float(ds.min() * 10),
        "single_point_energy_kj_mol": energy, "energy_finite": bool(np.isfinite(energy)),
        "max_force_kj_mol_nm": float(forces[worst]),
        "max_force_atom": f"{atoms[worst].residue.name}{atoms[worst].residue.index + 1}:{atoms[worst].name}",
        "mg_ligands": mg_now,
        "water_O_near_mg_A": [float(x * 10) for x in mg_tree_d[:6]],
        "openmm": mm.__version__,
    }


# ----------------------------------------------------------------------------- main build

MG_WATER_O_A = 2.08          # Mg-O(water) placement distance (crystal 5LPN/5SZJ: 1.94-2.44 A)
OPC_OH_A, OPC_HOH_DEG = 0.8724331, 103.6


def octahedral_vacancies(mg: np.ndarray, ligands: np.ndarray) -> tuple[list[np.ndarray], float]:
    """Best-fit octahedron (Kabsch over all ligand-to-vertex assignments); return the empty vertices."""
    unit_vectors = (ligands - mg) / np.linalg.norm(ligands - mg, axis=1)[:, None]
    vertices = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
    best = None
    for perm in itertools.permutations(range(6), len(unit_vectors)):
        ideal = vertices[list(perm)]
        u, _, vt = np.linalg.svd(ideal.T @ unit_vectors)
        d = np.sign(np.linalg.det(vt.T @ u.T))
        rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
        rmsd = float(np.sqrt(np.mean(np.sum((ideal @ rot.T - unit_vectors) ** 2, axis=1))))
        if best is None or rmsd < best[0]:
            best = (rmsd, rot, perm)
    rmsd, rot, perm = best
    return [vertices[i] @ rot.T for i in range(6) if i not in perm], rmsd


def mg_shell_waters(relaxed_pdb: Path, mg_pairs: list[dict], out_pdb: Path) -> dict:
    """Two OPC waters on the empty octahedral sites of Mg2+ (Rab10 crystal structures have exactly two)."""
    xyz, names, heavy = [], [], []
    for line in relaxed_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            names.append(f"{line[17:21].strip()}{int(line[22:26])}:{line[12:16].strip()}")
            heavy.append(not line[12:16].strip().startswith("H"))
    xyz = np.asarray(xyz)
    mg = xyz[mg_pairs[0]["mg_index"]]
    ligands = np.array([xyz[p["ligand_index"]] for p in mg_pairs])
    vacancies, rmsd = octahedral_vacancies(mg, ligands)
    if len(vacancies) != 2:
        raise ValueError(f"Expected two empty Mg sites, found {len(vacancies)} (ligands {len(mg_pairs)})")
    lines, report, placed = [], [], []
    half = math.radians(OPC_HOH_DEG / 2)
    others = np.array([i for i in range(len(xyz)) if heavy[i] and i != mg_pairs[0]["mg_index"]])
    for k, direction in enumerate(vacancies):
        o = mg + MG_WATER_O_A * direction
        helper = vacancies[1 - k] if abs(np.dot(vacancies[1 - k], direction)) < 0.9 else np.array([1.0, 0.0, 0.0])
        p = np.cross(direction, helper)
        p /= np.linalg.norm(p)
        h1 = o + OPC_OH_A * (math.cos(half) * direction + math.sin(half) * p)
        h2 = o + OPC_OH_A * (math.cos(half) * direction - math.sin(half) * p)
        d = np.linalg.norm(xyz[others] - o, axis=1)
        j = int(np.argmin(d))
        report.append({"water": k + 1, "mg_o_A": MG_WATER_O_A, "nearest_solute_heavy": names[others[j]],
                       "nearest_solute_heavy_A": float(d[j]),
                       "angles_to_ligands_deg": [angle(o, mg, l) for l in ligands]})
        placed.append(o)
        for name, x in [("O", o), ("H1", h1), ("H2", h2)]:
            lines.append(f"ATOM  {len(lines) + 1:5d} {atom_field(name, name[0])} WAT W{k + 1:4d}    "
                         f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          {name[0]:>2s}")
        lines.append("TER")
    if min(r["nearest_solute_heavy_A"] for r in report) < 2.2:  # 2.3-2.4 A = H-bond to D64 OD / GTP O2A, as in 5LPN
        raise ValueError(f"Mg shell water clashes with the solute: {report}")
    write_text(out_pdb, "\n".join(lines + ["END"]) + "\n")
    return {"octahedron_fit_rmsd_unit_vectors": rmsd, "waters": report,
            "water_water_O_A": float(np.linalg.norm(placed[0] - placed[1]))}


def mg_shell(residues: list[Residue], protocol: dict) -> list[dict]:
    spec = protocol["mg_restraints"]
    mg = next(r for r in residues if r.molecule == "MG").xyz("MG")
    shell = []
    for r in residues:
        if r.molecule == "MG":
            continue
        for name, (element, xyz) in r.atoms.items():
            d = float(np.linalg.norm(xyz - mg))
            if d <= spec["coordination_cutoff_A"]:
                shell.append({"residue": r.label, "molecule": r.molecule, "number": r.number, "atom": name,
                              "af3_distance_A": round(d, 3),
                              "flat_bottom_A": round(d + spec["flat_bottom_margin_A"], 3)})
    return shell


def ambertools_version() -> str | None:
    meta = Path(sys.prefix) / "conda-meta"
    hits = sorted(meta.glob("ambertools-*.json")) if meta.exists() else []
    return hits[0].name.removesuffix(".json") if hits else os.environ.get("AMBERTOOLS_VERSION")


def prepare(system: dict, protocol: dict, protocol_path: Path, af3_root: Path, out_root: Path,
            params_dir: Path, tleap: str, propka_python: str, force: bool) -> Path:
    out = out_root / system["system_id"]
    if (out / "qc.json").exists() and not force:
        raise FileExistsError(f"{out} already built (qc.json); use --force to rebuild")
    if out.exists() and force:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    rng_seed = system["seed"]
    cif = af3_root / system["af3_job"] / f"fold_{system['af3_job']}_model_{system['af3_model']}.cif"
    residues, seq_report = read_af3(cif, protocol)
    con = protocol["construct"]
    native = protocol["starts"][system["start"]]["residue_1122"]
    expected_kif = list(con["kif13a"]["sequence_wt"])
    if native == "GLU":
        expected_kif[1122 - con["kif13a"]["first"]] = "E"
    if seq_report["rab10_sequence"] != con["rab10"]["sequence"]:
        raise ValueError("Rab10 sequence differs from P61026 1-175")
    if seq_report["kif13a_sequence"] != "".join(expected_kif):
        raise ValueError("KIF13A sequence differs from the expected start sequence")
    if find(residues, "KIF13A", 1122).name != native:
        raise ValueError("Residue 1122 identity differs from the start definition")
    for position, name in [(1130, "GLU"), (1124, "GLU"), (1125, "ASP"), (1126, "ASP")]:
        if find(residues, "KIF13A", position).name != name:
            raise ValueError(f"KIF13A {position} is not {name}")
    if find(residues, "Rab10", 70).name != "ARG":
        raise ValueError("Rab10 70 is not ARG")
    mg_af3 = mg_shell(residues, protocol)
    mutation_report = None
    if system["mutation"]:
        mutation_report = rebuild_beyond_cb(residues, system["mutation"], rng_seed)
    residues, cap_checks = add_caps(residues, protocol)
    sg = [(r.label, r.xyz("SG")) for r in residues if r.name == "CYS"]
    disulfide_candidates = [f"{a}-{b} {np.linalg.norm(x - y):.2f}" for i, (a, x) in enumerate(sg)
                            for b, y in sg[i + 1:] if np.linalg.norm(x - y) < 2.5]
    if disulfide_candidates:
        raise ValueError(f"Unexpected disulfide geometry: {disulfide_candidates}")
    propka = run_propka(residues, out, propka_python)
    his = his_states(residues, propka, protocol)
    write_text(out / "complex_heavy.pdb", pdb_text(residues, his_names=True, leap=True))

    # --- dry topology for the restrained relaxation
    header = leap_header(protocol, params_dir)
    log_dry = run_leap(tleap, out, "dry_prerelax", "mol = loadPdb complex_heavy.pdb\ncheck mol\ncharge mol\n"
                       "savePdb mol dry_prerelax.pdb\nsaveAmberParm mol dry_prerelax.prmtop dry_prerelax.inpcrd\n", header)
    dry_charge = leap_charges(log_dry)
    dry_top = read_prmtop(out / "dry_prerelax.prmtop")
    rmap = residue_map(dry_top, residues, protocol)
    index_of = {(row["molecule"], row["uniprot_number"]): row["prmtop_residue_index"] for row in rmap if not row["cap"]}
    names = dry_top["ATOM_NAME"]
    gtp_row = next(row for row in rmap if row["molecule"] == "GTP")
    mg_row = next(row for row in rmap if row["molecule"] == "MG")
    gtp_atoms = names[gtp_row["first_atom"]:gtp_row["end_atom"]]
    gtp_heavy_after = [n for n in gtp_atoms if not n.startswith("H")]
    if sorted(gtp_heavy_after) != sorted(GTP_CCD_TO_PREP.values()) or len(gtp_atoms) != 44:
        raise ValueError(f"tleap GTP atoms unexpected: {gtp_atoms}")
    mg_pairs = []
    for item in mg_af3:
        if item["molecule"] == "GTP":
            row = gtp_row
        else:
            row = rmap[index_of[(item["molecule"], item["number"])]]
        ligand_index = row["first_atom"] + names[row["first_atom"]:row["end_atom"]].index(item["atom"])
        mg_pairs.append({**item, "mg_index": mg_row["first_atom"], "ligand_index": ligand_index,
                         "flat_bottom_nm": item["flat_bottom_A"] / 10.0})
    free_residue = index_of[("KIF13A", 1122)] if system["mutation"] else None
    tether_spec = protocol["arms"][system["arm"]].get("tether") or {}
    tether_report = None
    if tether_spec.get("enabled"):
        if not system["mutation"] or system["mutation"]["to"] != "ARG":
            raise ValueError("The tether option is defined for the K1122R paired rebuild only")
        r_row, e_row = rmap[free_residue], rmap[index_of[("KIF13A", 1130)]]
        n_idx = [i for i in range(r_row["first_atom"], r_row["end_atom"]) if names[i] in {"NE", "NH1", "NH2"}]
        o_idx = [i for i in range(e_row["first_atom"], e_row["end_atom"]) if names[i] in {"OE1", "OE2"}]
        wt_nz = np.asarray(mutation_report["original_side_chain_xyz_A"]["NZ"]) / 10.0
        relax = tethered_relaxation(out, protocol, free_residue, mg_pairs,
                                    {"spec": tether_spec, "guanidinium_n": n_idx, "e1130_o": o_idx, "target_nm": wt_nz})
        tether_report = relax["tether"]
    else:
        relax = restrained_relaxation(out, protocol, free_residue, mg_pairs)
    replace_pdb_coords(out / "dry_prerelax.pdb", relax["xyz_angstrom"], out / "complex_relaxed.pdb")
    bonds = set()
    bond_list = dry_top["BONDS_INC_HYDROGEN"] + dry_top["BONDS_WITHOUT_HYDROGEN"]
    for i in range(0, len(bond_list), 3):
        bonds.add((bond_list[i] // 3, bond_list[i + 1] // 3))
    clashes = clash_report(relax["xyz_angstrom"], relax["atoms"], index_of[("KIF13A", 1122)], bonds)
    site = {}
    res1122 = rmap[index_of[("KIF13A", 1122)]]
    local = {names[i]: relax["xyz_angstrom"][i] for i in range(res1122["first_atom"], res1122["end_atom"])}
    if "OE1" in local:
        site["carboxylate_OE1_OE2_A"] = float(np.linalg.norm(local["OE1"] - local["OE2"]))
        low, high = protocol["mutant_relax"]["carboxylate_oo_A"]
        if not low <= site["carboxylate_OE1_OE2_A"] <= high:
            raise ValueError(f"E1122 carboxylate O-O {site['carboxylate_OE1_OE2_A']:.3f} A outside [{low}, {high}]")
    if "CZ" in local:
        site["guanidinium_CZ_N_A"] = [float(np.linalg.norm(local["CZ"] - local[n])) for n in ["NE", "NH1", "NH2"]]

    def charged_xyz(molecule: str, number: int) -> np.ndarray:
        row = rmap[index_of[(molecule, number)]]
        atom_names = {**POS_ATOMS, **NEG_ATOMS}[row["resname_prmtop"]]
        return np.array([relax["xyz_angstrom"][i] for i in range(row["first_atom"], row["end_atom"]) if names[i] in atom_names])

    def min_d(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.min(np.linalg.norm(a[:, None] - b[None], axis=-1)))

    r70, e1130, r1122 = charged_xyz("Rab10", 70), charged_xyz("KIF13A", 1130), charged_xyz("KIF13A", 1122)
    site.update({"res1122": res1122["resname_prmtop"], "res1122_charged_to_R70_A": min_d(r1122, r70),
                 "res1122_charged_to_E1130_A": min_d(r1122, e1130), "E1130_to_R70_A": min_d(e1130, r70),
                 "note": "minimum distances between charged-group atoms after the restrained relaxation"})
    if tether_report is not None:
        r1122_n = np.array([relax["xyz_angstrom"][i] for i in range(rmap[free_residue]["first_atom"], rmap[free_residue]["end_atom"])
                            if names[i] in {"NE", "NH1", "NH2"}])
        wt_nz_A = np.asarray(tether_report["target_wt_k1122_nz_A"])
        tether_report["final"] = {
            "r1122_n_to_e1130_o_A": site["res1122_charged_to_E1130_A"], "r1122_n_to_r70_A": site["res1122_charged_to_R70_A"],
            "e1130_to_r70_A": site["E1130_to_R70_A"],
            "guanidinium_centroid_to_wt_k1122_nz_A": float(np.linalg.norm(r1122_n.mean(axis=0) - wt_nz_A)),
            "closest_r1122_n_to_wt_k1122_nz_A": float(np.min(np.linalg.norm(r1122_n - wt_nz_A, axis=1))),
            "guanidinium_CZ_N_A": site.get("guanidinium_CZ_N_A"),
            "residue_1122_clashes": clashes["residue_1122"]["clashes"],
            "residue_1122_min_distance_to_other_residue_A": clashes["residue_1122"]["min_distance_to_other_residue_A"]}
        if site["res1122_charged_to_E1130_A"] > tether_spec["accept"]["r1122_e1130_max_A"]:
            raise ValueError(f"Tether failed: R1122-E1130 {site['res1122_charged_to_E1130_A']:.2f} A after release")
    if system["mutation"] and clashes["residue_1122"]["clashes"]:
        raise ValueError(f"Clash at rebuilt residue 1122: {clashes['residue_1122']}")

    # --- rotation-safe truncated octahedron
    receptor_ids = {i for i, row in enumerate(rmap) if row["molecule"] in {"Rab10", "GTP", "MG"} and row["pdb_chain"] in {"A", "C", "D"}}
    ligand_ids = {i for i, row in enumerate(rmap) if row["molecule"] == "KIF13A"}
    leap_resseq = {}
    for line in (out / "dry_prerelax.pdb").read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            leap_resseq.setdefault(int(line[22:26]), None)
    seq_numbers = sorted(leap_resseq)
    if len(seq_numbers) != len(rmap):
        raise ValueError("tleap PDB residue count differs from the residue map")
    rec_set = {seq_numbers[i] for i in receptor_ids}
    lig_set = {seq_numbers[i] for i in ligand_ids}
    split_pdb_by_residue_names(out / "complex_relaxed.pdb", lambda name, seq: seq in rec_set, out / "receptor_relaxed.pdb")
    split_pdb_by_residue_names(out / "complex_relaxed.pdb", lambda name, seq: seq in lig_set, out / "ligand_relaxed.pdb")
    shell_report = mg_shell_waters(out / "complex_relaxed.pdb", mg_pairs, out / "mg_waters.pdb")
    from scipy.spatial.distance import pdist
    diameter = float(pdist(relax["xyz_angstrom"]).max())
    solv = protocol["solvation"]
    target_a = diameter + 2 * solv["buffer_A"]
    buffer = solv["buffer_A"]
    probes = []
    for attempt in range(8):
        log_probe = run_leap(tleap, out, "probe", "mol = loadPdb complex_relaxed.pdb\nw = loadPdb mg_waters.pdb\n"
                             f"mol = combine {{ mol w }}\nsolvateOct mol OPCBOX {buffer:.3f}\n"
                             "saveAmberParm mol probe.prmtop probe.inpcrd\n", header)
        box = read_inpcrd_box(out / "probe.inpcrd")
        n_wat = read_prmtop(out / "probe.prmtop")["RESIDUE_LABEL"].count("WAT")
        probes.append({"buffer_A": buffer, "box": box, "waters": n_wat})
        if not solv["rotation_safe"] or box[0] >= target_a:
            break
        buffer += max(0.5, (target_a - box[0]) * 0.75)
    else:
        raise RuntimeError("Could not reach a rotation-safe box")
    overlap_removed = boundary_overlap_waters(out / "probe.prmtop", out / "probe.inpcrd")
    for stale in ["probe.prmtop", "probe.inpcrd"]:
        (out / stale).unlink()
    solute_charge = round(dry_charge[-1])
    if abs(dry_charge[-1] - solute_charge) > 1e-3:
        raise ValueError(f"Non-integral solute charge {dry_charge[-1]}")
    n_neutral = abs(solute_charge)
    n_waters_kept = probes[-1]["waters"] - len(overlap_removed)
    n_pairs = int(math.floor((n_waters_kept - n_neutral) * solv["salt_molar"] / solv["water_molarity"] + 0.5))
    n_na = n_pairs + (n_neutral if solute_charge < 0 else 0)
    n_cl = n_pairs + (n_neutral if solute_charge > 0 else 0)
    body = ("mol = loadPdb complex_relaxed.pdb\nrec = loadPdb receptor_relaxed.pdb\nlig = loadPdb ligand_relaxed.pdb\n"
            "charge mol\ncharge rec\ncharge lig\n"
            "saveAmberParm mol complex_dry.prmtop complex_dry.inpcrd\n"
            "saveAmberParm rec receptor_dry.prmtop receptor_dry.inpcrd\n"
            "saveAmberParm lig ligand_dry.prmtop ligand_dry.inpcrd\n"
            "w = loadPdb mg_waters.pdb\nmol = combine { mol w }\n"
            f"solvateOct mol OPCBOX {buffer:.3f}\n"
            + "".join(f"remove mol mol.{r['leap_residue']}\n" for r in sorted(overlap_removed, key=lambda r: -r["leap_residue"]))
            + f"addIonsRand mol {solv['cation']} {n_na} {solv['anion']} {n_cl}\n"
            "charge mol\n"
            "saveAmberParm mol complex.prmtop complex.inpcrd\nsavePdb mol complex_solvated.pdb\n")
    log_final = run_leap(tleap, out, "final", body, header)
    final_charges = leap_charges(log_final)
    final_top = read_prmtop(out / "complex.prmtop")
    rmap_final = residue_map(final_top, residues, protocol)
    if [r["first_atom"] for r in rmap_final] != [r["first_atom"] for r in rmap]:
        raise ValueError("Solute atom layout changed between dry and solvated topologies")
    n_solute = rmap[-1]["end_atom"]
    dry_names = read_prmtop(out / "complex_dry.prmtop")["ATOM_NAME"]
    if final_top["ATOM_NAME"][:n_solute] != dry_names or len(dry_names) != n_solute:
        raise ValueError("Dry complex atoms do not match the solute block of the solvated topology")
    rec_top, lig_top = read_prmtop(out / "receptor_dry.prmtop"), read_prmtop(out / "ligand_dry.prmtop")
    rec_atoms = [i for row in rmap if row["prmtop_residue_index"] in receptor_ids for i in range(row["first_atom"], row["end_atom"])]
    lig_atoms = [i for row in rmap if row["prmtop_residue_index"] in ligand_ids for i in range(row["first_atom"], row["end_atom"])]
    if [dry_names[i] for i in rec_atoms] != rec_top["ATOM_NAME"] or [dry_names[i] for i in lig_atoms] != lig_top["ATOM_NAME"]:
        raise ValueError("Receptor/ligand dry topologies do not match the complex atom order")
    labels = final_top["RESIDUE_LABEL"]
    # tleap writes solute, then ions, then waters; find the two placed shell waters by their distance to Mg
    final_xyz = read_inpcrd_xyz(out / "complex.inpcrd")
    mg_atom = mg_pairs[0]["mg_index"]
    shell_waters = []
    for r, label in enumerate(labels):
        if label == "WAT":
            first = final_top["RESIDUE_POINTER"][r] - 1
            if np.linalg.norm(final_xyz[first] - final_xyz[mg_atom]) < MG_WATER_O_A + 0.05:
                shell_waters.append({final_top["ATOM_NAME"][i]: i for i in range(first, first + 4)})
    if len(shell_waters) != 2 or any(set(w) != {"O", "H1", "H2", "EPW"} for w in shell_waters):
        raise ValueError(f"Expected the two placed Mg shell waters in the solvated topology, found {shell_waters}")
    counts = {"WAT": labels.count("WAT"), "Na+": labels.count("Na+"), "Cl-": labels.count("Cl-")}
    if counts["Na+"] != n_na or counts["Cl-"] != n_cl:
        raise ValueError(f"Ion counts {counts} differ from requested Na+ {n_na}, Cl- {n_cl}")
    write_json(out / "mg_restraints.json", {
        "system_id": system["system_id"], "k_kj_mol_nm2": protocol["mg_restraints"]["k_kj_mol_nm2"],
        "form": protocol["mg_restraints"]["form"], "pairs": mg_pairs, "shell_waters_unrestrained": shell_waters,
        "ligand_atom_names_prmtop": [final_top["ATOM_NAME"][p["ligand_index"]] for p in mg_pairs]})
    for p in mg_pairs:
        if final_top["ATOM_NAME"][p["ligand_index"]] != p["atom"] or final_top["ATOM_NAME"][p["mg_index"]] != "MG":
            raise ValueError("Mg restraint atom mapping failed")
    qc_solv = solvated_qc(out, protocol, mg_pairs, n_solute)
    if qc_solv["water_oo_pairs_below_2.2A"] != 0:
        raise ValueError(f"Boundary water overlaps remain after removal: {qc_solv['water_oo_pairs_below_2.2A']}")
    if not (qc_solv["particles_match"] and qc_solv["virtual_sites_equal_waters"] and qc_solv["opc_spot_check_ok"]
            and abs(qc_solv["net_charge_e"]) < 1e-3 and qc_solv["energy_finite"]):
        raise ValueError(f"Solvated QC failed: {qc_solv}")
    rows_tsv = ["\t".join(["system_id"] + list(rmap_final[0]))]
    rows_tsv += ["\t".join([system["system_id"]] + ["" if v is None else str(v) for v in row.values()]) for row in rmap_final]
    write_text(out / "residue_map.tsv", "\n".join(rows_tsv) + "\n")
    receptor_rows = sorted(receptor_ids)
    write_json(out / "residue_map.json", {"system_id": system["system_id"], "residues": rmap_final,
                                          "n_solute_atoms": n_solute,
                                          "receptor_residue_indices": receptor_rows,
                                          "ligand_residue_indices": sorted(ligand_ids),
                                          "receptor_atoms": len(rec_atoms), "ligand_atoms": len(lig_atoms)})
    inputs = {"af3_cif": cif, "gtp_prep": params_dir / "GTP.prep", "frcmod_phos": params_dir / "frcmod.phos",
              "protocol": protocol_path, "prepare_complex.py": Path(__file__), "complex_common.py": HERE / "complex_common.py"}
    for key in ["rab10", "kif13a"]:
        fasta = ROOT.parent / con[key]["fasta"]
        if fasta.exists():
            inputs[f"{key}_fasta"] = fasta
    qc = {
        "system": system, "utc": utc_now(), "wall_s": None,
        "inputs_sha256": {k: sha256(v) for k, v in inputs.items()},
        "inputs_path": {k: str(v) for k, v in inputs.items()},
        "ambertools": ambertools_version(), "python": sys.version.split()[0],
        "af3_sequences": seq_report,
        "gtp_atom_map_ccd_to_prep": GTP_CCD_TO_PREP,
        "gtp_heavy_atoms_af3": 32, "gtp_atoms_after_tleap": len(gtp_atoms), "gtp_hydrogens_added_by_tleap": len(gtp_atoms) - 32,
        "mutation": mutation_report, "caps": cap_checks,
        "propka": propka, "his_states": his,
        "mg_shell_af3": mg_af3, "mg_shell_waters": shell_report,
        "relaxation": relax["summary"], "tether": tether_report, "clashes_after_relaxation": clashes, "site_1122_geometry": site,
        "charges": {"dry_solute_tleap_e": dry_charge[-1] if dry_charge else None,
                    "solvated_tleap_e": final_charges[-1] if final_charges else None,
                    "expected_components": "Rab10 1-175 (free N-term, NME) + KIF13A 1085-1160 (ACE/NME) + GTP -4 + Mg +2"},
        "solvation": {"solute_diameter_A": diameter, "target_box_a_A": target_a, "probes": probes,
                      "final_buffer_A": buffer, "box_inpcrd": read_inpcrd_box(out / "complex.inpcrd"),
                      "boundary_overlap_waters_removed": len(overlap_removed), "boundary_overlap_detail": overlap_removed,
                      "boundary_overlap_cutoff_A": BOUNDARY_OO_CUTOFF_A,
                      "salt_pairs": n_pairs, "neutralising": n_neutral, "Na+": n_na, "Cl-": n_cl,
                      "ion_placement": "tleap addIonsRand (random water replacement; seeded by tleap per invocation, not reproducible bit-for-bit)"},
        "counts": {"total_atoms_prmtop": final_top["POINTERS"][0], "solute_atoms": n_solute,
                   "protein_and_caps_atoms": gtp_row["first_atom"], "gtp_atoms": len(gtp_atoms), "mg": 1,
                   "waters": counts["WAT"], "water_sites": 4 * counts["WAT"], "Na+": counts["Na+"], "Cl-": counts["Cl-"]},
        "openmm_solvated_qc": qc_solv,
        "leap_warnings": {name: re.search(r"Warnings = (\d+)", (out / f"leap_{name}.log").read_text()).group(1)
                          for name in ["dry_prerelax", "final"]},
    }
    qc["wall_s"] = time.time() - t_start
    write_json(out / "qc.json", qc)
    for name in ["complex.prmtop", "complex.inpcrd", "complex_dry.prmtop", "receptor_dry.prmtop", "ligand_dry.prmtop",
                 "mg_restraints.json", "residue_map.json", "qc.json"]:
        provenance(out / name, list(inputs.values()), params={"system": system})
    print(f"{system['system_id']}: atoms {qc['counts']['total_atoms_prmtop']} waters {counts['WAT']} "
          f"Na+ {n_na} Cl- {n_cl} box_a {qc['solvation']['box_inpcrd'][0]:.2f} A "
          f"His {[h['state'] for h in his]} wall {qc['wall_s']:.0f} s", flush=True)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--af3-root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, default=ROOT / "systems")
    parser.add_argument("--params-dir", type=Path, default=ROOT / "params")
    parser.add_argument("--tleap", default=shutil.which("tleap") or "tleap")
    parser.add_argument("--propka-python", default=sys.executable)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--system", action="append")
    group.add_argument("--all", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    protocol = load_protocol(args.protocol)
    systems = (system_manifest(protocol, args.protocol) if args.all else
               [system_by_id(protocol, args.protocol, s) for s in args.system])
    failures = []
    for system in systems:
        try:
            prepare(system, protocol, args.protocol, args.af3_root, args.out_root, args.params_dir,
                    args.tleap, args.propka_python, args.force)
        except Exception as error:  # keep building the others, report at the end
            failures.append((system["system_id"], repr(error)))
            print(f"FAILED {system['system_id']}: {error!r}", file=sys.stderr, flush=True)
    if failures:
        raise SystemExit(f"{len(failures)} system(s) failed: {failures}")


if __name__ == "__main__":
    main()
