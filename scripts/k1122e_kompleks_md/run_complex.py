"""Run one prepared Rab10(GTP, Mg):KIF13A system: minimise, restrained heating/NVT/NPT, production.

    python run_complex.py --system-dir systems/WT_wt_s5 --run-dir runs/WT_wt_s5 --platform CUDA --auto
    python run_complex.py --system-dir systems/WT_wt_s5 --run-dir smoke/WT_wt_s5 --platform CPU --smoke

Conventions follow the monomer main arm (k1122e_md v3.1 run.py): OpenMM 8.3.1, PME 1.0 nm, HBonds, rigid
OPC, LangevinMiddle 310 K 1/ps, 2 fs, MC barostat 1 bar every 25 steps, CUDA mixed precision with
deterministic forces, 1 ns checkpoints, exact resume (DCD frame count checked against the checkpoint),
full-system frames every 100 ps written to a sibling "<run-root>_solvent/<id>" tree.

Complex-specific: the topology comes from tleap (complex.prmtop/.inpcrd, truncated octahedron, OPC with a
massless virtual site per water). Equilibration releases harmonic restraints stepwise (all solute heavy
atoms, then backbone + GTP + Mg only, then weak backbone, then none; stage list in complex_protocol.json).
Mg2+ flat-bottom restraints (mg_restraints.json) stay on for the whole run.

Trajectory: trajectory.dcd holds solute (protein, caps, GTP, Mg) + the two placed Mg-shell waters +
all Na+/Cl- every 10 ps; trajectory_atoms.json and trajectory_top.pdb describe it.

--auto (task files): no run.json -> fresh start; finished -> exit 0; production checkpoint -> resume;
interrupted equilibration -> partial files moved to aborted_<utc>/ and a fresh start.
A hard kill can leave DCD frames / state rows after the last checkpoint: on resume they are cut back to the
checkpoint (recorded in run.json "tail_repairs") and regenerated, so frame times stay exact.
Exit codes: 0 complete, 3 paused before the deadline (resubmit to continue), other = error.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import socket
import struct
import sys
import time
from pathlib import Path

from complex_common import DEFAULT_PROTOCOL, load_protocol, provenance, read_json, sha256, utc_now, write_json

PAUSED = 3
KEEP_ON_ARCHIVE = {"task.out", "task.err", "EXIT", "DONE", "run_complex.log"}


def steps_for(ns: float, timestep_fs: float) -> int:
    steps = ns * 1_000_000 / timestep_fs
    if abs(steps - round(steps)) > 1e-6:
        raise ValueError(f"{ns} ns is not an integral number of {timestep_fs} fs steps")
    return int(round(steps))


# ----------------------------------------------------------------------------- DCD / CSV tail repair

def dcd_frames(path: Path) -> tuple[int, int, int, int]:
    """(frames from file size, frames in header, atoms, frame bytes) for an OpenMM-written DCD."""
    with path.open("rb") as handle:
        head = handle.read(276)
    if struct.unpack("<i", head[:4])[0] != 84 or head[4:8] != b"CORD":
        raise ValueError(f"{path} is not a DCD file")
    nset = struct.unpack("<i", head[8:12])[0]
    box_flag = struct.unpack("<i", head[48:52])[0]
    comments = struct.unpack("<i", head[92:96])[0]
    header = 104 + comments + 8
    natoms = struct.unpack("<i", head[104 + comments:108 + comments])[0]
    frame = (56 if box_flag else 0) + 3 * (8 + 4 * natoms)
    size = path.stat().st_size
    return (size - header) // frame, nset, natoms, frame


def truncate_dcd(path: Path, keep: int) -> dict:
    with path.open("rb") as handle:
        head = handle.read(276)
    comments = struct.unpack("<i", head[92:96])[0]
    header = 104 + comments + 8
    first_step, interval = struct.unpack("<ii", head[12:20])
    full, nset, natoms, frame = dcd_frames(path)
    before = path.stat().st_size
    with path.open("r+b") as handle:
        handle.truncate(header + keep * frame)
        handle.seek(8)
        handle.write(struct.pack("<i", keep))
        handle.seek(20)
        handle.write(struct.pack("<i", first_step + max(keep - 1, 0) * interval))
    return {"file": path.name, "frames_before": full, "header_frames_before": nset, "frames_after": keep,
            "bytes_before": before, "bytes_after": path.stat().st_size}


def truncate_states(path: Path, last_step: int) -> dict:
    lines = path.read_text().splitlines(keepends=True)
    kept, dropped = [], 0
    for line in lines:
        if line.startswith("#") or not line.strip():
            kept.append(line)
            continue
        if int(float(line.split(",")[0])) <= last_step:
            kept.append(line)
        else:
            dropped += 1
    with path.open("w", newline="") as handle:
        handle.writelines(kept)
    return {"file": path.name, "rows_dropped": dropped}


# ----------------------------------------------------------------------------- system

def restraint_sets(topology, residue_map: dict) -> tuple[list[int], list[int]]:
    backbone, sidechain = [], []
    rows = residue_map["residues"]
    atoms = list(topology.atoms())
    for row in rows:
        for i in range(row["first_atom"], row["end_atom"]):
            atom = atoms[i]
            if atom.element is None or atom.element.symbol == "H":
                continue
            if row["molecule"] in {"GTP", "MG"} or row["cap"] or atom.name in {"N", "CA", "C", "O"}:
                backbone.append(i)
            else:
                sidechain.append(i)
    return backbone, sidechain


def build_system(system_dir: Path, protocol: dict, residue_map: dict, mg: dict, seed: int = 0):
    import openmm as mm
    from openmm import app, unit
    sim = protocol["simulation"]
    prmtop = app.AmberPrmtopFile(str(system_dir / "complex.prmtop"))
    inpcrd = app.AmberInpcrdFile(str(system_dir / "complex.inpcrd"))
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=sim["nonbonded_cutoff_nm"] * unit.nanometer,
                                 constraints=app.HBonds, rigidWater=sim["rigid_water"],
                                 ewaldErrorTolerance=sim["ewald_error_tolerance"])
    if system.getNumParticles() != prmtop.topology.getNumAtoms():
        raise ValueError("Particle count differs from the prmtop atom count")
    n_vs = sum(system.isVirtualSite(i) for i in range(system.getNumParticles()))
    n_wat = sum(1 for r in prmtop.topology.residues() if r.name in {"HOH", "WAT"})
    if n_vs != n_wat:
        raise ValueError(f"OPC virtual sites {n_vs} != waters {n_wat}")
    xyz = inpcrd.positions.value_in_unit(unit.nanometer)
    backbone, sidechain = restraint_sets(prmtop.topology, residue_map)
    for name, atoms in [("k_backbone", backbone), ("k_sidechain", sidechain)]:
        force = mm.CustomExternalForce(f"0.5*{name}*periodicdistance(x,y,z,x0,y0,z0)^2")
        force.addGlobalParameter(name, 0.0)
        for p in ("x0", "y0", "z0"):
            force.addPerParticleParameter(p)
        for i in atoms:
            force.addParticle(i, xyz[i])
        force.setForceGroup(30)
        system.addForce(force)
    mg_force = mm.CustomBondForce("step(r-rfb)*0.5*kmg*(r-rfb)^2")
    mg_force.addGlobalParameter("kmg", mg["k_kj_mol_nm2"])
    mg_force.addPerBondParameter("rfb")
    for pair in mg["pairs"]:
        mg_force.addBond(pair["mg_index"], pair["ligand_index"], [pair["flat_bottom_nm"]])
    mg_force.setUsesPeriodicBoundaryConditions(True)
    mg_force.setForceGroup(29)
    system.addForce(mg_force)
    barostat = mm.MonteCarloBarostat(sim["pressure_bar"] * unit.bar, sim["temperature_k"] * unit.kelvin, 0)
    barostat.setRandomNumberSeed(seed + 2 if seed else 0)  # stored in system.xml, as in the main arm build.py
    system.addForce(barostat)
    return system, prmtop, inpcrd, {"backbone": len(backbone), "sidechain": len(sidechain),
                                    "virtual_sites": n_vs, "waters": n_wat}


def trajectory_selection(prmtop, residue_map: dict, mg: dict) -> list[int]:
    solute = list(range(residue_map["n_solute_atoms"]))
    shell = sorted(i for water in mg["shell_waters_unrestrained"] for i in water.values())
    ions = [a.index for a in prmtop.topology.atoms() if a.residue.name in {"Na+", "Cl-", "NA", "CL"}]
    selection = sorted(set(solute + shell + ions))
    if len(selection) != len(solute) + len(shell) + len(ions):
        raise ValueError("Overlapping trajectory selections")
    return selection


def barostat_of(system):
    import openmm as mm
    return next(f for f in system.getForces() if isinstance(f, mm.MonteCarloBarostat))


def snapshot(sim, mg: dict) -> dict:
    import numpy as np
    from openmm import unit
    state = sim.context.getState(getEnergy=True, getPositions=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    restr = sim.context.getState(getEnergy=True, groups={30}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    mge = sim.context.getState(getEnergy=True, groups={29}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    inv = np.linalg.inv(np.asarray(box))
    dists = []
    for pair in mg["pairs"]:
        d = xyz[pair["ligand_index"]] - xyz[pair["mg_index"]]
        frac = d @ inv
        d = (frac - np.round(frac)) @ np.asarray(box)
        dists.append(round(float(np.linalg.norm(d)) * 10, 3))
    dof = sim.system_dof
    return {"step": sim.currentStep, "potential_kj_mol": pe, "kinetic_kj_mol": ke,
            "temperature_k": 2 * ke / (dof * 0.00831446261815324),
            "volume_nm3": float(abs(np.linalg.det(box))), "restraint_kj_mol": restr, "mg_restraint_kj_mol": mge,
            "mg_ligand_distances_A": dists, "finite": bool(np.isfinite(pe))}


def degrees_of_freedom(system) -> int:
    import openmm as mm
    dof = sum(3 for i in range(system.getNumParticles()) if system.getParticleMass(i)._value > 0)
    dof -= system.getNumConstraints()
    if any(isinstance(f, mm.CMMotionRemover) for f in system.getForces()):
        dof -= 3
    return dof


def archive_partial(run_dir: Path, solvent_dir: Path) -> Path:
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    target = run_dir / f"aborted_{stamp}"
    target.mkdir()
    for item in list(run_dir.iterdir()):
        if item.name in KEEP_ON_ARCHIVE or item.name.startswith("aborted_") or item == target:
            continue
        shutil.move(str(item), str(target / item.name))
    if solvent_dir.exists():
        for item in list(solvent_dir.iterdir()):
            shutil.move(str(item), str(target / f"solvent_{item.name}"))
    return target


def run(args) -> int:
    import numpy as np
    import openmm as mm
    from openmm import app, unit

    protocol = load_protocol(args.protocol)
    sim_cfg = protocol["simulation"]
    system_dir, run_dir = args.system_dir.resolve(), args.run_dir.resolve()
    qc = read_json(system_dir / "qc.json")
    meta = qc["system"]
    residue_map = read_json(system_dir / "residue_map.json")
    mg = read_json(system_dir / "mg_restraints.json")
    solvent_dir = run_dir.parent.parent / f"{run_dir.parent.name}_solvent" / run_dir.name
    run_path = run_dir / "run.json"
    run_dir.mkdir(parents=True, exist_ok=True)

    resume = args.resume
    if args.auto:
        if run_path.exists():
            previous = read_json(run_path)
            if previous.get("completed"):
                print(json.dumps({"system": meta["system_id"], "status": "already complete"}), flush=True)
                return 0
            if previous.get("phase") == "production":
                resume = True
            else:
                target = archive_partial(run_dir, solvent_dir)
                print(json.dumps({"system": meta["system_id"], "status": "interrupted equilibration archived",
                                  "to": str(target)}), flush=True)
    if run_path.exists() and not resume:
        raise FileExistsError("A run exists here; use --resume / --auto, or a fresh run directory")
    solvent_dir.mkdir(parents=True, exist_ok=True)

    timestep = sim_cfg["timestep_fs"]
    smoke = args.smoke
    sm = sim_cfg["smoke"]
    if smoke:
        production_steps = sm["production_steps"] if args.duration_ns is None else steps_for(args.duration_ns, timestep)
        interval, solvent_interval, checkpoint_interval = sm["report_steps"], sm["solvent_report_steps"], sm["checkpoint_steps"]
    else:
        production_steps = steps_for(args.duration_ns if args.duration_ns else sim_cfg["production_ns"], timestep)
        interval = steps_for(sim_cfg["report_ps"] / 1000, timestep)
        solvent_interval = steps_for(sim_cfg["solvent_report_ps"] / 1000, timestep)
        checkpoint_interval = steps_for(sim_cfg["checkpoint_ps"] / 1000, timestep)
    if solvent_interval % interval or checkpoint_interval % solvent_interval:
        raise ValueError("Report, solvent and checkpoint intervals must nest")

    platform = mm.Platform.getPlatformByName(args.platform)
    properties = dict(sim_cfg["cuda_properties"]) if args.platform == "CUDA" else (
        {"Threads": str(args.cpu_threads)} if args.platform == "CPU" else {})
    system_xml = run_dir / "system.xml"
    if resume:
        previous = read_json(run_path)
        system = mm.XmlSerializer.deserialize(system_xml.read_text(encoding="utf-8"))
        if sha256(system_xml) != previous["system_sha256"]:
            raise ValueError("system.xml changed since the run started")
        prmtop = app.AmberPrmtopFile(str(system_dir / "complex.prmtop"))
        inpcrd = app.AmberInpcrdFile(str(system_dir / "complex.inpcrd"))
        counts = previous["restraint_atom_counts"]
        barostat_of(system).setFrequency(sim_cfg["barostat_frequency"])
    else:
        system, prmtop, inpcrd, counts = build_system(system_dir, protocol, residue_map, mg, meta["seed"])
        system_xml.write_text(mm.XmlSerializer.serialize(system), encoding="utf-8")
    selection = trajectory_selection(prmtop, residue_map, mg)
    integrator = mm.LangevinMiddleIntegrator(sim_cfg["temperature_k"] * unit.kelvin,
                                             sim_cfg["friction_per_ps"] / unit.picosecond,
                                             timestep * unit.femtosecond)
    integrator.setRandomNumberSeed(meta["seed"] + 3)
    sim = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    sim.system_dof = degrees_of_freedom(system)

    runtime = {
        "system_id": meta["system_id"], "arm": meta["arm"], "start": meta["start"], "replica": meta["replica"],
        "pair_key": meta["pair_key"], "residue_1122": meta["residue_1122"], "seed": meta["seed"],
        "protocol_id": protocol["protocol_id"], "protocol_sha256": sha256(args.protocol),
        "prmtop_sha256": sha256(system_dir / "complex.prmtop"), "inpcrd_sha256": sha256(system_dir / "complex.inpcrd"),
        "system_sha256": sha256(system_xml), "platform": args.platform, "platform_properties": properties,
        "openmm_version": mm.__version__, "smoke": smoke, "duration_override_ns": args.duration_ns,
        "timestep_fs": timestep, "target_production_steps": production_steps,
        "report_interval_steps": interval, "solvent_report_interval_steps": solvent_interval,
        "checkpoint_interval_steps": checkpoint_interval,
        "trajectory_atoms": len(selection), "system_particles": system.getNumParticles(),
        "restraint_atom_counts": counts, "mg_restraints": {"k_kj_mol_nm2": mg["k_kj_mol_nm2"], "pairs": len(mg["pairs"])},
        "solvent_dir": str(solvent_dir), "completed": False,
    }
    keys_must_match = ["system_sha256", "protocol_sha256", "prmtop_sha256", "platform", "openmm_version",
                       "report_interval_steps", "solvent_report_interval_steps", "target_production_steps",
                       "duration_override_ns", "smoke", "trajectory_atoms"]
    segment = {"start_utc": utc_now(), "host": socket.gethostname(), "slurm_job": os.environ.get("SLURM_JOB_ID"),
               "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES")}
    wall0 = time.monotonic()
    trajectory_path, solvent_path, states_path = run_dir / "trajectory.dcd", solvent_dir / "system_solvent.dcd", run_dir / "states.csv"
    dcd_outputs = [(trajectory_path, interval, selection), (solvent_path, solvent_interval, None)]

    if resume:
        previous = read_json(run_path)
        if previous["phase"] not in {"production", "complete"}:
            raise ValueError("Only a production checkpoint can resume")
        for key in keys_must_match:
            if previous.get(key) != runtime[key]:
                raise ValueError(f"Resume environment or protocol mismatch: {key}: {previous.get(key)} != {runtime[key]}")
        sim.loadCheckpoint(str(run_dir / "checkpoint.bin"))
        production_start = previous["production_start_step"]
        completed = sim.currentStep - production_start
        repairs = previous.get("tail_repairs", [])
        for path, every, subset in dcd_outputs:
            expected = sim.currentStep // every - production_start // every
            frames, header_frames, natoms, frame_bytes = dcd_frames(path)
            if natoms != (len(subset) if subset else system.getNumParticles()):
                raise ValueError(f"{path.name}: atom count mismatch")
            if frames < expected:
                raise ValueError(f"{path.name} has {frames} frames, checkpoint needs {expected}; cannot resume exactly")
            if header_frames != expected or path.stat().st_size != 276 + expected * frame_bytes:
                repairs.append({"utc": utc_now(), "checkpoint_step": sim.currentStep, **truncate_dcd(path, expected)})
        rows = truncate_states(states_path, sim.currentStep)
        if rows["rows_dropped"]:
            repairs.append({"utc": utc_now(), "checkpoint_step": sim.currentStep, **rows})
        for key in ["started_utc", "production_start_step", "equilibration", "segments", "wall_seconds"]:
            runtime[key] = previous.get(key)
        runtime["tail_repairs"] = repairs
        runtime["prior_wall_seconds"] = previous.get("wall_seconds", 0.0)
    else:
        runtime["started_utc"] = utc_now()
        runtime["segments"] = []
        sim.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        runtime["phase"] = "equilibration"
        write_json(run_path, runtime)
        stages_log = []
        barostat_on = False
        for stage in sim_cfg["stages"]:
            t0 = time.monotonic()
            sim.context.setParameter("k_backbone", stage["k_backbone"])
            sim.context.setParameter("k_sidechain", stage["k_sidechain"])
            entry = {"name": stage["name"], "kind": stage["kind"], "k_backbone": stage["k_backbone"],
                     "k_sidechain": stage["k_sidechain"], "before": snapshot(sim, mg)}
            if stage["kind"] == "minimize":
                iterations = sm["minimize_iterations"] if smoke else stage["max_iterations"]
                mm.LocalEnergyMinimizer.minimize(sim.context, 10.0, iterations)
                entry["iterations"] = iterations
                first_md = next(s for s in sim_cfg["stages"] if s["kind"] != "minimize")
                sim.context.setVelocitiesToTemperature(first_md.get("t_start_k", sim_cfg["temperature_k"]) * unit.kelvin,
                                                       meta["seed"] + 4)
            else:
                steps = sm["stage_steps"] if smoke else steps_for(stage["ns"], timestep)
                if stage["kind"] == "npt" and not barostat_on:
                    barostat = barostat_of(system)
                    barostat.setFrequency(sim_cfg["barostat_frequency"])
                    sim.context.reinitialize(preserveState=True)
                    barostat_on = True
                if "t_start_k" in stage:
                    chunks = min(stage["ramp_chunks"], steps)
                    base, extra = divmod(steps, chunks)
                    for c in range(chunks):
                        temperature = stage["t_start_k"] + (stage["t_end_k"] - stage["t_start_k"]) * (c + 1) / chunks
                        integrator.setTemperature(temperature * unit.kelvin)
                        sim.step(base + (1 if c < extra else 0))
                    integrator.setTemperature(sim_cfg["temperature_k"] * unit.kelvin)
                else:
                    sim.step(steps)
                entry["steps"] = steps
            entry["after"] = snapshot(sim, mg)
            entry["wall_s"] = time.monotonic() - t0
            if not entry["after"]["finite"]:
                raise ValueError(f"Non-finite energy after stage {stage['name']}")
            stages_log.append(entry)
            print(json.dumps({"system": meta["system_id"], "stage": stage["name"], "step": sim.currentStep,
                              "T_K": round(entry["after"]["temperature_k"], 1),
                              "V_nm3": round(entry["after"]["volume_nm3"], 2),
                              "mg_A": entry["after"]["mg_ligand_distances_A"],
                              "wall_s": round(entry["wall_s"], 1)}), flush=True)
        if not barostat_on:
            raise ValueError("Protocol has no NPT stage before production")
        sim.context.setParameter("k_backbone", 0.0)
        sim.context.setParameter("k_sidechain", 0.0)
        production_start = sim.currentStep
        completed = 0
        runtime["equilibration"] = stages_log
        runtime["production_start_step"] = production_start
        sim.saveCheckpoint(str(run_dir / "checkpoint.bin"))
        state = sim.context.getState(getPositions=True)
        modeller = app.Modeller(prmtop.topology, state.getPositions())
        keep = set(selection)
        modeller.delete([a for a in modeller.topology.atoms() if a.index not in keep])
        with (run_dir / "trajectory_top.pdb").open("w") as handle:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=False)
        write_json(run_dir / "trajectory_atoms.json", {"system_id": meta["system_id"], "atom_indices": selection,
                                                      "prmtop": str(system_dir / "complex.prmtop"),
                                                      "content": "solute (protein, caps, GTP, Mg) + 2 Mg-shell waters + all Na+/Cl-"})
        runtime["tail_repairs"] = []
        runtime["prior_wall_seconds"] = 0.0

    runtime["phase"] = "production"
    runtime["completed_production_steps"] = completed
    runtime["segments"] = runtime.get("segments") or []
    runtime["segments"].append(segment)
    write_json(run_path, runtime)
    for path, every, subset in dcd_outputs:
        sim.reporters.append(app.DCDReporter(str(path), every, append=resume, enforcePeriodicBox=False, atomSubset=subset))
    state_handle = states_path.open("a" if resume else "w", newline="")
    sim.reporters.append(app.StateDataReporter(state_handle, interval, step=True, time=True, potentialEnergy=True,
                                               kineticEnergy=True, temperature=True, volume=True, speed=True,
                                               separator=",", append=resume))
    deadline = args.deadline_epoch or float(os.environ.get("COMPLEX_DEADLINE_EPOCH", "0") or 0)
    blocks_done, last_block_wall = 0, None
    status = 0
    while completed < production_steps:
        block = min(checkpoint_interval - (sim.currentStep - production_start) % checkpoint_interval,
                    production_steps - completed)
        if deadline and last_block_wall and time.time() + 1.3 * last_block_wall * block / checkpoint_interval > deadline:
            status = PAUSED
            break
        t0 = time.monotonic()
        sim.step(block)
        state = sim.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(energy):
            raise ValueError("Non-finite production potential energy")
        completed += block
        state_handle.flush()
        tmp = run_dir / "checkpoint.bin.tmp"
        sim.saveCheckpoint(str(tmp))
        tmp.replace(run_dir / "checkpoint.bin")
        last_block_wall = time.monotonic() - t0
        blocks_done += 1
        runtime["completed_production_steps"] = completed
        runtime["wall_seconds"] = time.monotonic() - wall0 + runtime.get("prior_wall_seconds", 0.0)
        runtime["segments"][-1].update({"end_utc": utc_now(), "last_step": sim.currentStep})
        ns_day = block * timestep * 1e-6 / (last_block_wall / 86400)
        write_json(run_path, runtime)
        print(json.dumps({"system": meta["system_id"], "production_ns": completed * timestep * 1e-6,
                          "target_ns": production_steps * timestep * 1e-6, "ns_per_day": round(ns_day, 1),
                          "mg_A": snapshot(sim, mg)["mg_ligand_distances_A"]}), flush=True)
        if args.abort_after_steps is not None and blocks_done == 1:
            sim.step(args.abort_after_steps)  # test hook: frames written past the checkpoint, then a hard kill
            state_handle.flush()
            for reporter in sim.reporters:
                if isinstance(reporter, app.DCDReporter) and reporter._dcd is not None:
                    reporter._out.flush()
            os._exit(17)
        if args.stop_after_blocks and blocks_done >= args.stop_after_blocks and completed < production_steps:
            status = PAUSED
            break
    state_handle.close()
    sim.reporters.clear()
    if status == PAUSED:
        runtime["segments"][-1]["paused"] = True
        write_json(run_path, runtime)
        print(json.dumps({"system": meta["system_id"], "status": "paused", "completed_steps": completed}), flush=True)
        return PAUSED
    state = sim.context.getState(getPositions=True, getVelocities=True)
    (run_dir / "final_state.xml").write_text(mm.XmlSerializer.serialize(state), encoding="utf-8")
    runtime["phase"] = "complete"
    runtime["completed"] = True
    runtime["completed_utc"] = utc_now()
    write_json(run_path, runtime)
    inputs = [system_dir / "complex.prmtop", system_dir / "complex.inpcrd", system_dir / "mg_restraints.json",
              system_dir / "residue_map.json", args.protocol, system_xml]
    for path in [trajectory_path, solvent_path, states_path, run_dir / "checkpoint.bin", run_dir / "final_state.xml",
                 run_dir / "trajectory_atoms.json", run_path]:
        provenance(path, inputs, params={k: runtime[k] for k in ["system_id", "seed", "platform", "smoke",
                                                                   "target_production_steps", "duration_override_ns"]})
    print(json.dumps({"system": meta["system_id"], "status": "complete", "steps": completed}), flush=True)
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--system-dir", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--platform", choices=["CUDA", "OpenCL", "CPU", "Reference"], default="CUDA")
    parser.add_argument("--cpu-threads", type=int, default=int(os.environ.get("OPENMM_CPU_THREADS", "1")))
    parser.add_argument("--smoke", action="store_true", help="engine test; excluded from scientific analysis")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--auto", action="store_true", help="fresh, resume or skip as appropriate (task files)")
    parser.add_argument("--duration-ns", type=float, help="non-standard production length (excluded from inference)")
    parser.add_argument("--deadline-epoch", type=float, help="pause cleanly before this UNIX time (default env COMPLEX_DEADLINE_EPOCH)")
    parser.add_argument("--stop-after-blocks", type=int, help="test hook: pause after N checkpoint blocks")
    parser.add_argument("--abort-after-steps", type=int, help="test hook: hard exit N steps after the first checkpoint")
    args = parser.parse_args()
    args.protocol = args.protocol.resolve()
    sys.exit(run(args))


if __name__ == "__main__":
    main()
