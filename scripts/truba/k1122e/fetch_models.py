"""Download exact official splice model weights with shared atlas provenance."""
from __future__ import annotations
import asyncio
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))
from atlas_mcp import core


async def main() -> None:
    sources = [
        ("spliceai", "illumina/SpliceAI", "03f42437aaf56dc5dfd822c4ccee5aec1a705079",
         [(f"spliceai/models/spliceai{i}.h5", f"spliceai{i}.h5") for i in range(1, 6)]),
        ("pangolin", "tkzeng/Pangolin", "5cf94b8db938c658391b4305cd7ce33297d44ff7",
         [(f"pangolin/models/final.{r}.{t}.3.v2", f"final.{r}.{t}.3.v2")
          for r in range(1, 4) for t in [0, 2, 4, 6]] + [("pangolin/model.py", "model.py")]),
    ]
    for name, repo, commit, files in sources:
        for source, filename in files:
            url = f"https://raw.githubusercontent.com/{repo}/{commit}/{source}"
            payload = await core.get_bytes(url)
            core.save_bytes(f"k1122e_models/{name}", filename, payload,
                            source_urls=[url], upstream_commit=commit,
                            purpose="academic research inference")
            print(f"Downloaded {name}/{filename}: {len(payload)} bytes", flush=True)


if __name__ == "__main__":
    asyncio.run(main())
