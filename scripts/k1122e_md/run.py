"""Run one built KIF13A task with restrained equilibration and unrestrained MD.

    python scripts/k1122e_md/run.py --task-dir RUN/primary_WT_r01 --platform CUDA
    python scripts/k1122e_md/run.py --task-dir RUN/primary_WT_r01 --platform CPU --smoke

The smoke mode is a short engine test and cannot enter scientific summaries.
Production checkpoints support exact continuation only when DCD frame count
matches the saved checkpoint. An inconsistent interrupted tail is rejected.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import time
from pathlib import Path

from common import provenance, read_json, sha256, write_json

SOLVENT_EVERY_REPORTS = 10  # full-system frame every 10 reports (100 ps at the 10 ps protocol interval)


def steps_ns(duration: float, timestep_fs: float) -> int:
    steps = duration * 1_000_000 / timestep_fs
    if abs(steps - round(steps)) > 1e-5:
        raise ValueError("Duration must be an integral number of integration steps")
    return int(round(steps))


def run(task_dir: Path, platform_name: str, *, smoke: bool = False,
        resume: bool = False, duration_ns: float | None = None,
        cpu_threads: int = 1) -> None:
    import numpy as np
    import openmm as mm
    from openmm import app, unit

    build_path = task_dir / "build.json"
    build = read_json(build_path)
    settings = build["config"]["simulation"]
    task = build["task"]
    run_path = task_dir / "run.json"
    if smoke and resume:
        raise ValueError("Smoke jobs cannot resume")
    if run_path.exists() and not resume:
        raise FileExistsError("A run exists. Use --resume for a production checkpoint, or a fresh task directory.")
    system_path = task_dir / "system.xml"
    system = mm.XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
    if resume:
        next(force for force in system.getForces()
             if isinstance(force, mm.MonteCarloBarostat)).setFrequency(25)
    pdb = app.PDBFile(str(task_dir / "system.pdb"))
    integ = mm.LangevinMiddleIntegrator(settings["temperature_k"] * unit.kelvin,
                                       settings["friction_per_ps"] / unit.picosecond,
                                       settings["timestep_fs"] * unit.femtosecond)
    integ.setRandomNumberSeed(task["seed"] + 3)
    platform = mm.Platform.getPlatformByName(platform_name)
    properties = {}
    if platform_name == "CUDA":
        properties = {"Precision": "mixed", "DeterministicForces": "true"}
    elif platform_name == "CPU":
        properties = {"Threads": str(cpu_threads)}
    sim = app.Simulation(pdb.topology, system, integ, platform, properties)
    target_ns = task["production_ns"] if duration_ns is None else duration_ns
    if target_ns <= 0:
        raise ValueError("Production duration must be positive")
    production_steps = 100 if smoke else steps_ns(target_ns, settings["timestep_fs"])
    interval = 10 if smoke else int(round(settings["report_ps"] * 1000 / settings["timestep_fs"]))
    # v3.1: also keep ions every report and the full solvated system every 10th report (100 ps),
    # so Na+/Cl- pairing and water-mediated bridges can be measured after the run.
    solvent_interval = interval * SOLVENT_EVERY_REPORTS
    ion_atoms = [atom.index for atom in pdb.topology.atoms() if atom.residue.name in {"NA", "CL"}]
    # The full-system DCD (~0.7 GB per 100 ns) goes to a sibling "<run-root>_solvent" tree so the frozen
    # summary job, which tars the whole run root into /arf/home, does not copy it.
    solvent_dir = task_dir.parent.parent / f"{task_dir.parent.name}_solvent" / task_dir.name
    solvent_dir.mkdir(parents=True, exist_ok=True)
    dcd_outputs = [(task_dir / "trajectory.dcd", interval, build["protein_atom_indices"]),
                   (task_dir / "protein_ions.dcd", interval, build["protein_atom_indices"] + ion_atoms),
                   (solvent_dir / "system_solvent.dcd", solvent_interval, None)]
    checkpoint_interval = max(interval, int(round(settings["checkpoint_ps"] * 1000 / settings["timestep_fs"])))
    runtime = {
        "task_id": task["task_id"], "platform": platform_name,
        "system_sha256": sha256(system_path), "config_sha256": task["config_sha256"],
        "seed": task["seed"], "openmm_version": mm.__version__,
        "smoke": smoke, "duration_override_ns": duration_ns,
        "target_production_steps": production_steps,
        "report_interval_steps": interval,
        "solvent_report_interval_steps": solvent_interval,
        "dcd_outputs": {str(path): {"interval_steps": every,
                                    "atoms": len(subset) if subset else system.getNumParticles()}
                        for path, every, subset in dcd_outputs},
        "timestep_fs": settings["timestep_fs"],
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "completed": False,
    }
    start_wall = time.monotonic()
    if resume:
        previous = read_json(run_path)
        if previous["smoke"] or previous["phase"] not in {"production", "complete"}:
            raise ValueError("Only a production checkpoint can resume")
        for key in ["system_sha256", "config_sha256", "platform", "openmm_version", "report_interval_steps",
                    "solvent_report_interval_steps", "target_production_steps", "duration_override_ns"]:
            if previous.get(key) != runtime[key]:
                raise ValueError(f"Resume environment or protocol mismatch: {key}")
        sim.loadCheckpoint(str(task_dir / "checkpoint.bin"))
        production_start = previous["production_start_step"]
        completed_steps = sim.currentStep - production_start
        from mdtraj.formats import DCDTrajectoryFile
        for path, every, _ in dcd_outputs:
            expected_frames = sim.currentStep // every - production_start // every
            with DCDTrajectoryFile(str(path), "r") as trajectory:
                if len(trajectory) != expected_frames:
                    raise ValueError(f"{path.name} tail does not match checkpoint; unsafe append rejected")
        runtime["started_utc"] = previous["started_utc"]
        runtime["prior_wall_seconds"] = previous.get("wall_seconds", 0)
    else:
        sim.context.setPositions(pdb.positions)
        sim.minimizeEnergy(maxIterations=20 if smoke else settings["minimization_max_iterations"])
        sim.context.setVelocitiesToTemperature(settings["temperature_k"] * unit.kelvin, task["seed"] + 4)
        runtime["phase"] = "equilibration"
        write_json(run_path, runtime)
        sim.step(10 if smoke else steps_ns(settings["nvt_ns"], settings["timestep_fs"]))
        barostat = next(force for force in system.getForces() if isinstance(force, mm.MonteCarloBarostat))
        barostat.setFrequency(25)
        sim.context.reinitialize(preserveState=True)
        sim.step(10 if smoke else steps_ns(settings["restrained_npt_ns"], settings["timestep_fs"]))
        sim.context.setParameter("k", 0.0)
        sim.step(10 if smoke else steps_ns(settings["unrestrained_npt_ns"], settings["timestep_fs"]))
        production_start = sim.currentStep
        completed_steps = 0
        sim.saveCheckpoint(str(task_dir / "checkpoint.bin"))
    runtime["production_start_step"] = production_start
    runtime["phase"] = "production"
    runtime["completed_production_steps"] = completed_steps
    write_json(run_path, runtime)
    for path, every, subset in dcd_outputs:
        sim.reporters.append(app.DCDReporter(str(path), every, append=resume,
                                            enforcePeriodicBox=False, atomSubset=subset))
    state_handle = (task_dir / "states.csv").open("a" if resume else "w")
    sim.reporters.append(app.StateDataReporter(state_handle, interval, step=True, time=True,
                                              potentialEnergy=True, kineticEnergy=True,
                                              temperature=True, volume=True, speed=True,
                                              separator=","))
    while completed_steps < production_steps:
        block = min(checkpoint_interval, production_steps - completed_steps)
        sim.step(block)
        state = sim.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(energy):
            raise ValueError("Nonfinite production potential energy")
        completed_steps += block
        checkpoint_temp = task_dir / "checkpoint.bin.tmp"
        sim.saveCheckpoint(str(checkpoint_temp))
        checkpoint_temp.replace(task_dir / "checkpoint.bin")
        runtime["completed_production_steps"] = completed_steps
        runtime["wall_seconds"] = time.monotonic() - start_wall + runtime.get("prior_wall_seconds", 0)
        write_json(run_path, runtime)
        print(json.dumps({"task": task["task_id"],
                          "production_ns": completed_steps * settings["timestep_fs"] / 1_000_000,
                          "target_ns": production_steps * settings["timestep_fs"] / 1_000_000,
                          "wall_seconds": runtime["wall_seconds"]}), flush=True)
    state_handle.close()
    sim.reporters.clear()
    state = sim.context.getState(getPositions=True, getVelocities=True)
    with (task_dir / "final_state.xml").open("w") as handle:
        handle.write(mm.XmlSerializer.serialize(state))
    runtime["phase"] = "complete"
    runtime["completed"] = True
    runtime["completed_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    write_json(run_path, runtime)
    inputs = [build_path, system_path]
    for path in [task_dir / "trajectory.dcd", task_dir / "protein_ions.dcd", solvent_dir / "system_solvent.dcd",
                 task_dir / "states.csv", task_dir / "checkpoint.bin", task_dir / "final_state.xml",
                 task_dir / "run.json"]:
        provenance(path, inputs, params=runtime)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", type=Path, required=True)
    parser.add_argument("--platform", choices=["CPU", "Reference", "CUDA", "OpenCL"], default="CUDA")
    parser.add_argument("--cpu-threads", type=int, default=1)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--duration-ns", type=float)
    args = parser.parse_args()
    run(args.task_dir, args.platform, smoke=args.smoke, resume=args.resume,
        duration_ns=args.duration_ns, cpu_threads=args.cpu_threads)


if __name__ == "__main__":
    main()
