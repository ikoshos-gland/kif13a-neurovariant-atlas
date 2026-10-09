"""Monomer v3 onizleme: suren bir gorevin KOPYASINI kesme anina kadar kirpar ve 'tamamlandi' diye isaretler.

Gercek kosuya dokunmaz. Cikti yalniz analyze.py'nin calisabilmesi icindir; preview_only alani ve
eligible_for_inference=False (hedef adim sayisi degismedigi icin) cikarim disi oldugunu kaydeder.
Kullanim: python onizleme_monomer.py <kaynak_gorev_dizini> <hedef_dizin> <kesme_ns>
"""
import json
import shutil
import sys
from pathlib import Path

import mdtraj as md

src, dst, cut_ns = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
dst.mkdir(parents=True, exist_ok=True)
for name in ("build.json", "protein.pdb"):
    shutil.copy2(src / name, dst / name)
run = json.loads((src / "run.json").read_text())
start, iv, dt_fs = run["production_start_step"], run["report_interval_steps"], run["timestep_fs"]
cut_steps = int(round(cut_ns * 1_000_000 / dt_fs))
n_need = (start + cut_steps) // iv - start // iv
frames = []
n = 0
for chunk in md.iterload(str(src / "trajectory.dcd"), top=str(src / "protein.pdb"), chunk=500):
    take = min(len(chunk), n_need - n)
    if take > 0:
        frames.append(chunk[:take])
        n += take
    if n >= n_need:
        break
if n < n_need:
    raise SystemExit(f"yetersiz kare: {n} < {n_need}")
traj = md.join(frames) if len(frames) > 1 else frames[0]
traj.save_dcd(str(dst / "trajectory.dcd"))
run.update({"completed": True, "completed_production_steps": cut_steps,
            "preview_only": f"onizleme: suren kosunun {cut_ns:g} ns'de kirpilmis kopyasi, cikarim icin degil",
            "preview_source_completed_steps": run.get("completed_production_steps")})
(dst / "run.json").write_text(json.dumps(run, indent=1))
print(f"{src.name}: {n} kare = {cut_ns:g} ns")
