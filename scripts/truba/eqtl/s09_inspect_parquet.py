"""Print the schema, row-group layout and first rows of one Zeng 2026 parquet file (compute node)."""
import sys

sys.path.insert(0, "/arf/scratch/avural/kif13a_eqtl/pylib")
import pyarrow.parquet as pq  # noqa: E402

p = sys.argv[1]
pf = pq.ParquetFile(p)
md = pf.metadata
print("rows", md.num_rows, "row_groups", md.num_row_groups, "columns", md.num_columns)
print(pf.schema_arrow)
print("key-value metadata keys:", list((pf.schema_arrow.metadata or {}).keys())[:10])
for rg in (0, 1, md.num_row_groups // 2, md.num_row_groups - 1):
    g = md.row_group(rg)
    stats = []
    for c in range(g.num_columns):
        col = g.column(c)
        s = col.statistics
        if s is not None and s.has_min_max and c < 6:
            stats.append(f"{col.path_in_schema}:[{s.min}..{s.max}]")
    print("rg", rg, "rows", g.num_rows, " ".join(stats))
t = pf.read_row_group(0).slice(0, 5).to_pandas()
import pandas as pd  # noqa: E402
with pd.option_context("display.width", 250, "display.max_columns", 40):
    print(t)
