"""Parameterize each variant and boundary before any production MD is launched.

    python scripts/k1122e_md/validate_caps.py --output cap_validation.json

This engine check creates no production trajectory. It validates exact cached
AF-derived neutral caps, mutation chemistry and finite configured potentials.
The unsolvated potential checks parameterization only, not stability.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from common import DEFAULT_CONFIG, HERE, ROOT, load_config, provenance, sha256, task_manifest, write_json
from build import add_amber_cap_hydrogens, engine, prepare_capped_fragment


def validate(config_path: Path, input_root: Path, platform_name: str = "CPU") -> dict:
    np, mm, app, unit, _ = engine()
    config = load_config(config_path)
    source_path = input_root / config["source"]["structure"]
    source = app.PDBxFile(str(source_path))
    rows = []
    for task in task_manifest(config, config_path):
        if task["replica"] != 1:
            continue
        random.seed(task["seed"])
        np.random.seed(task["seed"])
        fragment, checks = prepare_capped_fragment(source, task)
        # Force the initial-load state to catch the original lazy-load bug.
        app.Modeller._hasLoadedStandardHydrogens = False
        app.Modeller._residueHydrogens = {}
        ff = app.ForceField(*config["simulation"]["forcefield"])
        platform = mm.Platform.getPlatformByName(platform_name)
        add_amber_cap_hydrogens(fragment, ff, config["simulation"]["ph"], platform)
        system = ff.createSystem(fragment.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        integ = mm.VerletIntegrator(0.001 * unit.picosecond)
        properties = {"Threads": "1"} if platform_name == "CPU" else {}
        context = mm.Context(system, integ, platform, properties)
        context.setPositions(fragment.positions)
        potential = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(potential):
            raise ValueError(f"Nonfinite cap validation potential: {task['task_id']}")
        mutated = task["mutation"]
        names = {int(r.id): r.name for r in fragment.topology.residues()}
        if mutated and names[mutated["position"]] != mutated["to"]:
            raise ValueError("Prepared mutation does not match frozen task")
        caps = [r for r in fragment.topology.residues() if r.name in {"ACE", "NME"}]
        row = {"task_id": task["task_id"], "boundary": task["boundary"],
               "variant": task["variant"], "atoms": system.getNumParticles(),
               "initial_potential_kj_mol": float(potential),
               "caps": {r.name: sorted(a.name for a in r.atoms()) for r in caps},
               "cap_bond_checks": checks, "passed": True}
        rows.append(row)
        print(json.dumps({"validated": task["task_id"], "atoms": system.getNumParticles(),
                          "finite_potential_kj_mol": float(potential)}), flush=True)
        del context, integ
    return {"openmm_version": mm.__version__, "source_sha256": sha256(source_path),
            "forcefield": config["simulation"]["forcefield"],
            "forcefield_files_sha256": {name: sha256(Path(app.__file__).parent / "data" / name)
                                        for name in config["simulation"]["forcefield"]},
            "cap_definition_sha256": sha256(HERE / "cap_hydrogens.xml"),
            "config_sha256": sha256(config_path), "cases": rows,
            "scientific_inference": "none; parameterization check only"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--input-root", type=Path, default=ROOT)
    parser.add_argument("--platform", choices=["CPU", "Reference", "CUDA", "OpenCL"], default="CPU")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = validate(args.config, args.input_root, args.platform)
    if args.output:
        write_json(args.output, result)
        provenance(args.output, [args.config, args.input_root / load_config(args.config)["source"]["structure"],
                                 HERE / "cap_hydrogens.xml"], params={"validation_cases": len(result["cases"])})
    print(json.dumps({"validation_cases": len(result["cases"]), "passed": True}))


if __name__ == "__main__":
    main()
