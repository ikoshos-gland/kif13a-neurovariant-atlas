"""Connect to every server in .mcp.json the way Claude Code does, list its tools and
run one real KIF13A query. Run with the atlas venv (it has the MCP client SDK):

    mcp\\atlas_server\\.venv\\Scripts\\python.exe mcp\\smoke_test.py
    ... smoke_test.py --servers gnomad,string          # subset
    ... smoke_test.py --list protein                   # tool names + argument schemas
    ... smoke_test.py --call gnomad gnomad_get_gene_constraint '{"gene":"KIF13A"}'
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
import time
from contextlib import asynccontextmanager
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.client.streamable_http import streamablehttp_client

ROOT = Path(__file__).resolve().parents[1]

# One real query per server; save=False keeps the smoke test from writing data.
TESTS: dict[str, tuple[str, dict]] = {
    "kif13a-atlas": ("gnomad_constraint", {"gene": "KIF13A", "save": False}),
    "gnomad": ("gnomad_get_gene_constraint", {"gene": "KIF13A"}),
    "biomcp": ("get", {"entity": "gene", "id": "KIF13A"}),
    "opentargets": ("search_entities", {"query_strings": ["KIF13A"]}),
    "ensembl": ("ensembl_get_homology", {"symbol": "KIF13A", "species": "homo_sapiens", "target_species": "mus_musculus", "max_results": 3}),
    "uniprot": ("uniprot_get_entry", {"accessions": ["Q9H1H9"], "fields": "accession,length,gene_names"}),
    "interpro": ("getProteinEntries", {"uniprotAccession": "Q9H1H9"}),
    "protein": ("protein_get_annotations", {"uniprot": "Q9H1H9", "limit": 5}),
    "pdbe": ("get_best_structures_observed_api_uniprot_best_structures__un", {"uniprot_accession": "Q12756"}),
    "string": ("string_resolve_proteins", {"proteins": "KIF13A", "species": "9606"}),
    "tooluniverse": ("execute_tool", {"tool_name": "gnomad_get_gene_constraints", "arguments": {"gene_symbol": "KIF13A"}}),
    "pubmed": ("search_articles", {"query": "KIF13A serotonin 1A receptor anxiety", "max_results": 3}),
    "biorxiv": ("get_categories", {}),
    "zotero": ("", {}),
}


@asynccontextmanager
async def connect(cfg: dict):
    if cfg.get("type") == "http":
        async with streamablehttp_client(cfg["url"]) as (read, write, _):
            async with ClientSession(read, write) as session:
                yield session
    else:
        params = StdioServerParameters(
            command=cfg["command"], args=cfg.get("args", []), env={**os.environ, **cfg.get("env", {})}
        )
        # SMOKE_ERRLOG=path keeps the server's stderr for debugging
        with open(os.environ.get("SMOKE_ERRLOG", os.devnull), "a", encoding="utf-8") as errlog:
            async with stdio_client(params, errlog=errlog) as (read, write):
                async with ClientSession(read, write) as session:
                    yield session


def _text(result) -> str:
    parts = []
    for item in result.content:
        parts.append(getattr(item, "text", None) or f"<{type(item).__name__}>")
    return " ".join(parts)


async def run_one(name: str, cfg: dict, mode: str, call: tuple[str, dict] | None, timeout: float) -> dict:
    t0 = time.monotonic()
    report: dict = {"server": name}
    try:
        async def body():
            async with connect(cfg) as session:
                await session.initialize()
                tools = (await session.list_tools()).tools
                report["n_tools"] = len(tools)
                if mode == "list":
                    report["tools"] = {t.name: {"desc": (t.description or "")[:160],
                                                "args": list((t.inputSchema or {}).get("properties", {}).keys()),
                                                "required": (t.inputSchema or {}).get("required", [])}
                                       for t in tools}
                    return
                report["sample_tools"] = [t.name for t in tools[:8]]
                tool, args = call or TESTS.get(name, ("", {}))
                if tool:
                    result = await session.call_tool(tool, args)
                    report["call"] = tool
                    report["is_error"] = bool(result.isError)
                    report["result_head"] = _text(result)[:700]
        await asyncio.wait_for(body(), timeout)
        report["status"] = "ERROR_RESULT" if report.get("is_error") else "OK"
    except Exception as exc:  # report and move on to the next server
        report["status"] = "FAIL"
        report["error"] = f"{type(exc).__name__}: {exc}"[:500]
    report["seconds"] = round(time.monotonic() - t0, 1)
    return report


async def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--servers", help="comma-separated subset")
    ap.add_argument("--list", metavar="SERVER", help="print tool names and argument schemas for one server")
    ap.add_argument("--call", nargs=3, metavar=("SERVER", "TOOL", "JSON_ARGS"))
    ap.add_argument("--timeout", type=float, default=300)
    ap.add_argument("--out", help="also write the JSON report to this file")
    a = ap.parse_args()

    project_cfg = ROOT / ".mcp.json"
    if project_cfg.exists():
        servers = json.loads(project_cfg.read_text(encoding="utf-8"))["mcpServers"]
    else:  # user scope, registered by setup_mcp.py; skip servers this project did not add
        user = json.loads((Path.home() / ".claude.json").read_text(encoding="utf-8")).get("mcpServers", {})
        servers = {k: v for k, v in user.items() if k in TESTS}
    if a.list:
        rep = await run_one(a.list, servers[a.list], "list", None, a.timeout)
        print(json.dumps(rep, indent=1, ensure_ascii=False))
        return 0
    if a.call:
        srv, tool, args = a.call
        rep = await run_one(srv, servers[srv], "call", (tool, json.loads(args)), a.timeout)
        print(json.dumps(rep, indent=1, ensure_ascii=False))
        return 0 if rep["status"] == "OK" else 1

    names = a.servers.split(",") if a.servers else list(servers)
    reports = []
    for name in names:
        rep = await run_one(name, servers[name], "call", None, a.timeout)
        reports.append(rep)
        head = rep.get("result_head", rep.get("error", ""))[:160].replace("\n", " ")
        print(f"{rep['status']:<12} {name:<13} tools={rep.get('n_tools', '-'):<4} {rep['seconds']:>6}s  {head}", flush=True)
    if a.out:
        Path(a.out).write_text(json.dumps(reports, indent=1, ensure_ascii=False), encoding="utf-8")
    return 0 if all(r["status"] == "OK" for r in reports) else 1


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
