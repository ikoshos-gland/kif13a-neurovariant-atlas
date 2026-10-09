"""Build one independently seeded, ACE/NME-capped KIF13A monomer MD task.

    python scripts/k1122e_md/build.py --task-index 0 --output-root data/interim/k1122e_md

Uses the cached AlphaFold mmCIF. Caps copy adjacent backbone geometry and use
the standard Amber ACE/NME residue templates; connectivity and bond distances
are validated before force-field assignment. No binding or folding ddG is made.
"""

from __future__ import annotations

import argparse
import io
import json
import math
import random
from pathlib import Path
from typing import Any

from common import DEFAULT_CONFIG, HERE, ROOT, provenance, sha256, task_at, write_json


def engine():
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    from pdbfixer import PDBFixer
    return np, mm, app, unit, PDBFixer


def cap_fragment(fixer, source, start: int, end: int):
    """Use neighboring AF peptide geometry to construct standard neutral caps.

    ACE carbonyl C/O come from residue start-1 and CH3 from its CA. NME N
    comes from end+1 and CH3 from its CA. Internal coordinates and peptide
    bonds are inherited, rather than invented as arbitrary Cartesian offsets.
    """
    np, _, app, unit, _ = engine()
    source_residues = {int(r.id): r for r in source.topology.residues()}
    source_xyz = source.positions.value_in_unit(unit.nanometer)
    source_atoms = {i: {a.name: a for a in r.atoms()} for i, r in source_residues.items()}
    fixed_xyz = fixer.positions.value_in_unit(unit.nanometer)
    residues = list(fixer.topology.residues())
    if [int(r.id) for r in residues] != list(range(start, end + 1)):
        raise ValueError("Fragment numbering changed during heavy-atom repair")
    top = app.Topology()
    chain = top.addChain("A")
    atom_map = {}
    coords = []
    cap_atoms = {}
    ace = top.addResidue("ACE", chain, id=str(start - 1))
    for name, original in [("CH3", "CA"), ("C", "C"), ("O", "O")]:
        src = source_atoms[start - 1][original]
        cap_atoms[("ACE", name)] = top.addAtom(name, src.element, ace)
        coords.append(source_xyz[src.index])
    for residue in residues:
        new_residue = top.addResidue(residue.name, chain, id=residue.id)
        for atom in residue.atoms():
            if atom.name == "OXT" or atom.element.symbol == "H":
                continue
            atom_map[atom] = top.addAtom(atom.name, atom.element, new_residue)
            coords.append(fixed_xyz[atom.index])
    nme = top.addResidue("NME", chain, id=str(end + 1))
    for name, original in [("N", "N"), ("CH3", "CA")]:
        src = source_atoms[end + 1][original]
        cap_atoms[("NME", name)] = top.addAtom(name, src.element, nme)
        coords.append(source_xyz[src.index])
    for a, b in fixer.topology.bonds():
        if a in atom_map and b in atom_map:
            top.addBond(atom_map[a], atom_map[b])
    first_n = next(atom_map[a] for a in residues[0].atoms() if a.name == "N")
    last_c = next(atom_map[a] for a in residues[-1].atoms() if a.name == "C")
    top.addBond(cap_atoms[("ACE", "CH3")], cap_atoms[("ACE", "C")])
    top.addBond(cap_atoms[("ACE", "C")], cap_atoms[("ACE", "O")])
    top.addBond(cap_atoms[("ACE", "C")], first_n)
    top.addBond(last_c, cap_atoms[("NME", "N")])
    top.addBond(cap_atoms[("NME", "N")], cap_atoms[("NME", "CH3")])
    xyz = np.asarray(coords)
    checks = []
    for a, b, low, high in [
        (cap_atoms[("ACE", "C")], first_n, 0.12, 0.15),
        (last_c, cap_atoms[("NME", "N")], 0.12, 0.15),
        (cap_atoms[("ACE", "CH3")], cap_atoms[("ACE", "C")], 0.14, 0.17),
        (cap_atoms[("NME", "N")], cap_atoms[("NME", "CH3")], 0.13, 0.17),
    ]:
        distance = float(np.linalg.norm(xyz[a.index] - xyz[b.index]))
        checks.append({"atoms": [f"{a.residue.name}:{a.name}", f"{b.residue.name}:{b.name}"],
                       "distance_nm": distance, "accepted_nm": [low, high]})
        if not low <= distance <= high:
            raise ValueError(f"Invalid cap bond geometry: {checks[-1]}")
    return app.Modeller(top, xyz * unit.nanometer), checks


def rotate_to_principal_axes(modeller) -> None:
    """Rigidly align the fragment before measuring its common cubic box edge."""
    np, _, _, unit, _ = engine()
    xyz = np.array(modeller.positions.value_in_unit(unit.nanometer))
    indices = [a.index for a in modeller.topology.atoms() if a.name == "CA"]
    centered = xyz - xyz[indices].mean(axis=0)
    _, _, axes = np.linalg.svd(centered[indices], full_matrices=False)
    if np.linalg.det(axes) < 0:
        axes[-1] *= -1
    modeller.positions = (centered @ axes.T) * unit.nanometer


def cubic_box_lengths(xyz, settings: dict):
    """Give all orientations the same initial extent plus solvent padding."""
    import numpy as np
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError("Expected finite Cartesian atom coordinates")
    if settings["box_policy"] != "cubic_max_initial_atom_extent_plus_padding":
        raise ValueError("Only the frozen cubic solvent box policy is supported")
    edge = max(float(np.ptp(xyz, axis=0).max()) + 2 * settings["solvent_padding_nm"],
               2 * settings["nonbonded_cutoff_nm"] + 0.2)
    return np.full(3, edge)


def add_amber_cap_hydrogens(modeller, forcefield, ph: float, platform) -> None:
    """Initialize standard chemistry before overriding the Amber NME methyl."""
    from openmm import app
    # Directly loading hydrogens.xml leaves OpenMM's lazy-load flag false.
    # Its later automatic load would overwrite NME's custom CH3 parent.
    app.Modeller._loadStandardHydrogenDefinitions()
    app.Modeller.loadHydrogenDefinitions(str(HERE / "cap_hydrogens.xml"))
    modeller.addHydrogens(forcefield, pH=ph, platform=platform)


def prepare_capped_fragment(source, task: dict):
    """Repair only heavy atoms and apply a verified mutation before capping."""
    _, _, app, _, PDBFixer = engine()
    start, end = task["residues"]
    crop = app.Modeller(source.topology, source.positions)
    crop.delete([r for r in crop.topology.residues() if not start <= int(r.id) <= end])
    buffer = io.StringIO()
    app.PDBFile.writeFile(crop.topology, crop.positions, buffer, keepIds=True)
    buffer.seek(0)
    fixer = PDBFixer(pdbfile=buffer)
    mutation = task["mutation"]
    if mutation:
        actual = next(r.name for r in fixer.topology.residues() if int(r.id) == mutation["position"])
        if actual != mutation["from"]:
            raise ValueError(f"Reference residue mismatch: expected {mutation['from']}, found {actual}")
        fixer.applyMutations([f"{mutation['from']}-{mutation['position']}-{mutation['to']}"], "A")
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=task["seed"])
    modeller, checks = cap_fragment(fixer, source, start, end)
    rotate_to_principal_axes(modeller)
    return modeller, checks


def add_restraints(system, topology, positions, k: float, *, periodic: bool = True):
    _, mm, _, unit, _ = engine()
    expression = ("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2" if periodic else
                  "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint = mm.CustomExternalForce(expression)
    restraint.addGlobalParameter("k", k)
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    xyz = positions.value_in_unit(unit.nanometer)
    for atom in topology.atoms():
        if atom.residue.name not in {"HOH", "NA", "CL"} and atom.name in {"N", "CA", "C", "O"}:
            restraint.addParticle(atom.index, xyz[atom.index])
    restraint.setForceGroup(31)
    system.addForce(restraint)
    return restraint


def static_proxy(modeller, config: dict, platform_name: str) -> tuple[Any, dict]:
    """Evaluate an actual implicit-solvent local potential, never a folding ddG."""
    np, mm, app, unit, _ = engine()
    settings = config["static_proxy"]
    ionic = config["simulation"]["salt_molar"]
    temperature = config["simulation"]["temperature_k"]
    kappa = 367.434915 * math.sqrt(ionic / (settings["solvent_dielectric"] * temperature))
    forcefield = app.ForceField(*settings["forcefield"])
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
        soluteDielectric=settings["solute_dielectric"],
        solventDielectric=settings["solvent_dielectric"],
        implicitSolventKappa=kappa / unit.nanometer,
    )
    add_restraints(system, modeller.topology, modeller.positions, settings["restraint_k_kj_mol_nm2"], periodic=False)
    integ = mm.VerletIntegrator(0.001 * unit.picosecond)
    context = mm.Context(system, integ, mm.Platform.getPlatformByName(platform_name))
    context.setPositions(modeller.positions)
    groups = set(range(31))
    initial = context.getState(getEnergy=True, groups=groups).getPotentialEnergy()
    mm.LocalEnergyMinimizer.minimize(context, maxIterations=settings["minimization_max_iterations"])
    state = context.getState(getPositions=True, getEnergy=True, groups=groups)
    final = state.getPotentialEnergy()
    positions = state.getPositions(asNumpy=True)
    if not np.isfinite(final.value_in_unit(unit.kilojoule_per_mole)):
        raise ValueError("Nonfinite monomer potential")
    summary = {
        "initial_unrestrained_potential_kj_mol": float(initial.value_in_unit(unit.kilojoule_per_mole)),
        "final_unrestrained_potential_kj_mol": float(final.value_in_unit(unit.kilojoule_per_mole)),
        "minimized_with_backbone_restraints": True,
        "gb_screening_per_nm": kappa,
        "interpretation": settings["interpretation"],
        "is_folding_ddg": False,
        "is_binding_energy": False,
    }
    del context, integ
    return positions, summary


def verify_water_model(topology, system, settings: dict) -> dict:
    """Reject a mislabeled OPC system by inspecting its actual particles."""
    import openmm as mm
    from openmm import unit
    waters = [r for r in topology.residues() if r.name == "HOH"]
    if not waters:
        raise ValueError("No solvent waters found")
    result = {"physical_model": settings["water_model"],
              "coordinate_template": settings.get("solvent_template_model", settings["water_model"]),
              "water_molecules": len(waters),
              "virtual_sites": sum(system.isVirtualSite(i) for i in range(system.getNumParticles()))}
    if settings["water_model"] != "opc":
        return result
    nonbonded = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    expected_charges = {"O": 0.0, "H1": 0.679142, "H2": 0.679142, "M": -1.358284}
    weights = [0.7045587810585758, 0.1477206094707121, 0.1477206094707121]
    for water in waters:
        atoms = {a.name: a.index for a in water.atoms()}
        if set(atoms) != set(expected_charges):
            raise ValueError("OPC needs O, H1, H2 and M in every water")
        for name, charge in expected_charges.items():
            actual = nonbonded.getParticleParameters(atoms[name])[0].value_in_unit(unit.elementary_charge)
            if abs(actual - charge) > 1e-8:
                raise ValueError(f"Incorrect OPC charge on {name}: {actual}")
        index = atoms["M"]
        if not system.isVirtualSite(index) or system.getParticleMass(index).value_in_unit(unit.dalton) != 0:
            raise ValueError("OPC M must be a massless virtual site")
        site = system.getVirtualSite(index)
        if not isinstance(site, mm.ThreeParticleAverageSite):
            raise ValueError("Incorrect OPC virtual-site type")
        for j, name in enumerate(["O", "H1", "H2"]):
            if site.getParticle(j) != atoms[name] or abs(site.getWeight(j) - weights[j]) > 1e-10:
                raise ValueError("Incorrect OPC virtual-site geometry")
    if result["virtual_sites"] != len(waters):
        raise ValueError("OPC virtual-site count does not match water count")
    result.update({"sites_per_water": 4, "charges_e": expected_charges,
                   "virtual_site_weights": weights, "all_waters_passed": True})
    return result


def build(config_path: Path, task_index: int, output_root: Path,
          input_root: Path, platform_name: str = "CPU") -> Path:
    np, mm, app, unit, PDBFixer = engine()
    config, task = task_at(config_path, task_index)
    out = output_root / task["task_id"]
    if (out / "build.json").exists():
        raise FileExistsError(f"A completed build already exists: {out}")
    out.mkdir(parents=True, exist_ok=True)
    random.seed(task["seed"])
    np.random.seed(task["seed"])
    source_path = input_root / config["source"]["structure"]
    source = app.PDBxFile(str(source_path))
    chains = list(source.topology.chains())
    if len(chains) != 1 or chains[0].id != config["source"]["chain"]:
        raise ValueError("Expected exactly the verified AlphaFold chain")
    residues = list(source.topology.residues())
    if len(residues) != config["source"]["protein_length"]:
        raise ValueError("Full-length AlphaFold numbering must be verified before cropping")
    start, end = task["residues"]
    modeller, cap_checks = prepare_capped_fragment(source, task)
    sim = config["simulation"]
    forcefield = app.ForceField(*sim["forcefield"])
    add_amber_cap_hydrogens(modeller, forcefield, sim["ph"], mm.Platform.getPlatformByName(platform_name))
    proxy = {"enabled": False, "reason": config["static_proxy"].get("disabled_reason")}
    output_names = ["system.xml", "system.pdb", "protein.pdb", "build.json"]
    if config["static_proxy"].get("enabled", True):
        minimized, proxy = static_proxy(modeller, config, platform_name)
        # Separate diagnostic; the common AF backbone remains the MD start.
        with (out / "monomer_minimized.pdb").open("w") as handle:
            app.PDBFile.writeFile(modeller.topology, minimized, handle, keepIds=True)
        write_json(out / "static_proxy.json", proxy)
        output_names += ["monomer_minimized.pdb", "static_proxy.json"]
    with (out / "protein.pdb").open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
    residue_map = [{"residue_index": r.index, "residue_name": r.name,
                    "mane_position": int(r.id) if r.name not in {"ACE", "NME"} else None,
                    "cap": r.name in {"ACE", "NME"}}
                   for r in modeller.topology.residues()]
    protein_atoms = [a.index for a in modeller.topology.atoms()]
    box = cubic_box_lengths(modeller.positions.value_in_unit(unit.nanometer), sim)
    modeller.addSolvent(forcefield, model=sim.get("solvent_template_model", sim["water_model"]),
                       boxSize=mm.Vec3(*box) * unit.nanometer,
                       ionicStrength=sim["salt_molar"] * unit.molar,
                       positiveIon="Na+", negativeIon="Cl-", neutralize=True)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,
                                    nonbondedCutoff=sim["nonbonded_cutoff_nm"] * unit.nanometer,
                                    constraints=app.HBonds, rigidWater=True,
                                    ewaldErrorTolerance=sim["ewald_error_tolerance"])
    solvent_validation = verify_water_model(modeller.topology, system, sim)
    cmap_terms = sum(f.getNumTorsions() for f in system.getForces() if isinstance(f, mm.CMAPTorsionForce))
    if "amber19/protein.ff19SB.xml" in sim["forcefield"] and cmap_terms == 0:
        raise ValueError("Missing ff19SB backbone CMAP torsions")
    barostat = mm.MonteCarloBarostat(sim["pressure_bar"] * unit.bar,
                                   sim["temperature_k"] * unit.kelvin, 25)
    barostat.setRandomNumberSeed(task["seed"] + 2)
    barostat.setFrequency(0)
    system.addForce(barostat)
    add_restraints(system, modeller.topology, modeller.positions, sim["restraint_k_kj_mol_nm2"])
    (out / "system.xml").write_text(mm.XmlSerializer.serialize(system), encoding="utf-8")
    with (out / "system.pdb").open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
    metadata = {
        "task": task,
        "source_sha256": sha256(source_path),
        "config": config,
        "cap_bond_checks": cap_checks,
        "protein_atom_indices": protein_atoms,
        "protein_residue_map": residue_map,
        "protein_atoms": len(protein_atoms),
        "system_atoms": system.getNumParticles(),
        "box_nm": box.tolist(),
        "openmm_version": mm.__version__,
        "forcefield_files_sha256": {name: sha256(Path(app.__file__).parent / "data" / name)
                                    for name in sim["forcefield"]},
        "ff19sb_cmap_torsions": cmap_terms,
        "solvent_validation": solvent_validation,
        "cap_hydrogen_definition_sha256": sha256(HERE / "cap_hydrogens.xml"),
        "cap_hydrogen_standard_loader": "Modeller._loadStandardHydrogenDefinitions",
        "static_proxy": proxy,
    }
    confidence_path = input_root / config["source"]["confidence"]
    confidence = json.loads(confidence_path.read_text(encoding="utf-8"))
    entries = confidence if isinstance(confidence, list) else [confidence]
    values = entries[0].get("confidenceScore", entries[0].get("confidence_score", []))
    if isinstance(values, list) and len(values) == len(residues):
        selected = np.asarray(values[start - 1:end], dtype=float)
        core_start, core_end = config["analysis"]["core"]
        core = np.asarray(values[core_start - 1:core_end], dtype=float)
        metadata["af_confidence"] = {"fragment_min_plddt": float(selected.min()),
                                     "fragment_mean_plddt": float(selected.mean()),
                                     "core_min_plddt": float(core.min()),
                                     "core_mean_plddt": float(core.mean()),
                                     "site_plddt": float(values[1121])}
    else:
        raise ValueError("AlphaFold confidence schema or residue count unexpected")
    write_json(out / "build.json", metadata)
    for name in output_names:
        provenance(out / name, [config_path, source_path, confidence_path, HERE / "cap_hydrogens.xml"], params=task)
    print(json.dumps({"built": str(out), "atoms": system.getNumParticles(), "task": task["task_id"]}))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, default=ROOT)
    parser.add_argument("--platform", choices=["CPU", "Reference", "CUDA", "OpenCL"], default="CPU")
    args = parser.parse_args()
    build(args.config, args.task_index, args.output_root, args.input_root, args.platform)


if __name__ == "__main__":
    main()
