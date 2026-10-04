"""WP1 core data freeze: pull every raw input the variant universe and paralog map need,
through the same code the kif13a-atlas MCP server exposes (so provenance is identical).

    mcp\\atlas_server\\.venv\\Scripts\\python.exe scripts\\pull_wp1_core.py
    ... pull_wp1_core.py --only clinvar,alphamissense
"""

from __future__ import annotations

import argparse
import asyncio
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core, expression, structure, variants  # noqa: E402

GENES = ["KIF13A", "KIF1A", "KIF13B", "KIF1B", "KIF1C"]


def steps() -> dict[str, list]:
    t = {g: core.target(g) for g in GENES}
    return {
        "constraint": [(variants.gnomad_constraint, g) for g in GENES],
        "gnomad": [(variants.gnomad_variants, g, "coding") for g in ("KIF13A", "KIF1A")],
        "clinvar": [(variants.clinvar_variants, g) for g in ("KIF13A", "KIF1A")],
        "alphamissense": [(variants.alphamissense, t[g]["uniprot"]) for g in GENES],
        # features are annotated on UniProt canonical numbering (KIF1A: apply isoform_offsets)
        "uniprot": [(structure.uniprot_features, t[g]["uniprot"]) for g in GENES],
        # MANE isoform models: for KIF1A that is Q12756-3; canonical Q12756 kept for AlphaMissense alignment
        "alphafold": [(structure.alphafold_model, a) for a in ("Q9H1H9", "Q12756-3", "Q12756", "Q9NQT8")],
        "pdb": [(structure.pdb_structures, t[g]["uniprot"]) for g in ("KIF13A", "KIF1A", "KIF13B")],
        "paralog": [(structure.paralog_map, "KIF1A", "KIF13A", "1-1791"), (structure.paralog_map, "KIF13B", "KIF13A", "1-1826")],
        "expression": [
            *[(expression.gtex_expression, g) for g in ("KIF13A", "KIF1A", "KIF13B")],
            *[(expression.hpa_summary, g) for g in ("KIF13A", "KIF1A", "KIF13B")],
            *[(expression.cellxgene_expression, g) for g in ("KIF13A", "KIF1A", "KIF13B")],
        ],
    }


async def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated step names: " + ",".join(steps()))
    args = ap.parse_args()
    plan = steps()
    names = args.only.split(",") if args.only else list(plan)
    failures = 0
    for name in names:
        for fn, *fargs in plan[name]:
            t0 = time.monotonic()
            label = f"{name}:{fn.__name__}({', '.join(map(str, fargs))[:40]})"
            try:
                result = await fn(*fargs)
                saved = result.get("saved")
                files = [s["path"] for s in saved] if isinstance(saved, list) else ([saved["path"]] if saved else [])
                print(f"ok    {label:<60} {time.monotonic() - t0:6.1f}s  {', '.join(files)[:110]}", flush=True)
            except Exception as exc:
                failures += 1
                print(f"FAIL  {label:<60} {type(exc).__name__}: {exc}"[:300], flush=True)
    print(json.dumps({"files_in_inventory": len(core.inventory(10_000)), "failures": failures}))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
