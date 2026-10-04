"""Extract target genes from the Zeng ... Roussos 2026 class-level cis-eQTL parquet files.

Schema (s09_inspect_parquet.py, meta-Oligo-class.parquet): Source, Level, CellType, chr (double),
start, end (end = start + 1), variant (rsID), nonassessed_allele, assessed_allele, gene (symbol),
beta, se, z, isTopHit. 101,486,848 rows in 97 row groups for Oligo. Rows are grouped by gene,
so a pyarrow filter on `gene` touches only the row groups whose min/max statistics include it.
Output: work/zeng_<Class>.tsv.gz
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, "/arf/scratch/avural/kif13a_eqtl/pylib")
import pyarrow.compute as pc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

W = "/arf/scratch/avural/kif13a_eqtl"
SP = f"{W}/raw/synapse"
TARGETS = ["KIF13A", "NUP153", "NHLRC1", "TPMT",      # locus genes
           "CR1", "ADAMTS4", "NDUFS2", "APH1B"]       # positive-control genes (AD x Oligo)


def main() -> int:
    files = sorted(f for f in os.listdir(SP) if f.endswith("-class.parquet"))
    for f in files:
        cls = f.split("-")[1]
        out = f"{W}/work/zeng_{cls}.tsv.gz"
        if os.path.exists(out):
            print(f"have {out}", flush=True)
            continue
        t0 = time.time()
        pf = pq.ParquetFile(f"{SP}/{f}")
        keep = []
        gi = pf.schema_arrow.get_field_index("gene")
        for rg in range(pf.num_row_groups):
            st = pf.metadata.row_group(rg).column(gi).statistics
            if st is not None and st.has_min_max and not any(st.min <= g <= st.max for g in TARGETS):
                continue
            t = pf.read_row_group(rg)
            m = pc.is_in(t.column("gene"), value_set=__import__("pyarrow").array(TARGETS))
            sub = t.filter(m)
            if sub.num_rows:
                keep.append(sub.to_pandas())
        import pandas as pd
        d = pd.concat(keep) if keep else pd.DataFrame()
        d.to_csv(out, sep="\t", index=False, compression="gzip")
        cnt = d.groupby("gene").size().to_dict() if len(d) else {}
        print(f"{cls}: {len(d)} rows {cnt} ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
