"""Write multigpu_run.sh task files for the 12 complex runs (one task = one run_complex.py --auto call).

    python make_complex_tasks.py --root /arf/scratch/avural/so3lr/complex_md --jobs 6
    python make_complex_tasks.py --root ... --jobs 1 --systems WT_wt_s5 K1122E_paired_wt_s5 --tag pilot --duration-ns 1

Task line format (bin/multigpu_run.sh): "<workdir> <command>", the command runs via bash -c inside <workdir>
with CUDA_VISIBLE_DEVICES set to one GPU. <workdir> is the run directory itself, so multigpu_run.sh's
DONE marker sits next to run.json and a finished run is skipped on resubmission. ${USER} and
${SLURM_JOB_ID} are expanded by bash at run time; COMPLEX_DEADLINE_EPOCH is exported by the SLURM script.
Runs are sorted by system size and cut into consecutive chunks, so the GPUs of one job finish at similar times.
"""

from __future__ import annotations

import argparse
import json
import shlex
from pathlib import Path

from complex_common import DEFAULT_PROTOCOL, load_protocol, system_manifest, write_text

SIF_DEFAULT = "/arf/home/avural/containers/mdana-openmm831-pyscf.sif"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--root", required=True, help="complex_md root on TRUBA (POSIX path)")
    parser.add_argument("--sif", default=SIF_DEFAULT)
    parser.add_argument("--runs-name", default="runs", help="run root under --root; full frames go to <runs-name>_solvent")
    parser.add_argument("--jobs", type=int, default=6, help="number of task files (one SLURM job each)")
    parser.add_argument("--systems", nargs="*", help="subset of system ids (default: all 12)")
    parser.add_argument("--tag", default="prod")
    parser.add_argument("--duration-ns", type=float, help="non-standard production length (pilot; excluded from inference)")
    parser.add_argument("--local-systems", type=Path, help="local copy of systems/ to read atom counts for balancing")
    args = parser.parse_args()
    protocol = load_protocol(args.protocol)
    manifest = system_manifest(protocol, args.protocol)
    if args.systems:
        unknown = set(args.systems) - {m["system_id"] for m in manifest}
        if unknown:
            raise SystemExit(f"unknown systems: {sorted(unknown)}")
        manifest = [m for m in manifest if m["system_id"] in args.systems]
    root = args.root.rstrip("/")
    scripts, systems, runs = f"{root}/scripts", f"{root}/systems", f"{root}/{args.runs_name}"
    sizes = {}
    for m in manifest:
        qc = (args.local_systems or Path(__file__).resolve().parent.parent / "systems") / m["system_id"] / "qc.json"
        sizes[m["system_id"]] = json.loads(qc.read_text())["counts"]["total_atoms_prmtop"] if qc.exists() else 0
    ordered = sorted(manifest, key=lambda m: -sizes[m["system_id"]])
    n_jobs = max(1, min(args.jobs, len(ordered)))
    per_job = -(-len(ordered) // n_jobs)
    # consecutive chunks of the size-sorted list: the runs of one job (one per GPU, in parallel) have similar
    # sizes, so neither GPU idles long while the other finishes
    jobs = [ordered[k:k + per_job] for k in range(0, len(ordered), per_job)]
    out_dir = Path(__file__).resolve().parent.parent / "tasks"
    written = []
    for j, members in enumerate(jobs, start=1):
        lines = [f"# {args.tag} job {j}/{len(jobs)}: {len(members)} runs; resubmit the same file to continue",
                 f"# protocol {protocol['protocol_id']}"]
        for m in members:
            run_dir = f"{runs}/{m['system_id']}"
            extra = f" --duration-ns {args.duration_ns:g}" if args.duration_ns else ""
            command = (f"apptainer exec --nv -B /arf --home /tmp/${{USER}}_home_${{SLURM_JOB_ID}} {shlex.quote(args.sif)} "
                       f"python {scripts}/run_complex.py --protocol {scripts}/complex_protocol.json "
                       f"--system-dir {systems}/{m['system_id']} --run-dir {run_dir} --platform CUDA --auto{extra} "
                       f">> {run_dir}/run_complex.log 2>&1")
            lines.append(f"{run_dir} {command}")
        path = out_dir / f"complex_{args.tag}_job{j:02d}.txt"
        write_text(path, "\n".join(lines) + "\n")
        written.append({"file": str(path), "systems": [m["system_id"] for m in members],
                        "atoms": [sizes[m["system_id"]] for m in members]})
    print(json.dumps(written, indent=1))


if __name__ == "__main__":
    main()
