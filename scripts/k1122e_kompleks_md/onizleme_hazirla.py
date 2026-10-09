"""Önizleme: kopyalanmış yarım koşunun run.json'unu DCD'deki kare sayısına göre 'tamamlanmış' gibi işaretler.
Yalnız onizleme/ altındaki KOPYAYI değiştirir; gerçek koşuya dokunmaz. Sonuçlar çıkarım için uygun değildir."""
import json, sys
from pathlib import Path
import mdtraj as md
run = Path(sys.argv[1]); sysdir = Path(sys.argv[2])
r = json.loads((run / "run.json").read_text())
sel = json.loads((run / "trajectory_atoms.json").read_text())["atom_indices"]
top = md.load_prmtop(str(sysdir / "complex.prmtop")).subset(sel)
n = 0
for chunk in md.iterload(str(run / "trajectory.dcd"), top=top, chunk=50):
    n += len(chunk)
start, iv = r["production_start_step"], r["report_interval_steps"]
r["completed_production_steps"] = (start // iv + n) * iv - start
r["completed"] = True
r["preview_only"] = "onizleme: yarim kosunun kopyasi, cikarim icin degil"
(run / "run.json").write_text(json.dumps(r, indent=1))
print(f"{n} kare = {n * iv * r['timestep_fs'] * 1e-6:.2f} ns")
