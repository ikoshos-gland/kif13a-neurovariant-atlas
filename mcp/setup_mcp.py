"""Register the MCP server catalogue below with Claude Code.

Run after installing the servers (see docs/MCP_KURULUM.md) and again whenever
you add an API key to .env. Use the atlas venv python, not Microsoft Store python:

    mcp\\atlas_server\\.venv\\Scripts\\python.exe mcp\\setup_mcp.py                  # user scope (every session)
    mcp\\atlas_server\\.venv\\Scripts\\python.exe mcp\\setup_mcp.py --scope project  # only sessions opened in this folder
    mcp\\atlas_server\\.venv\\Scripts\\python.exe mcp\\setup_mcp.py --check          # report only

User scope goes through `claude mcp add-json --scope user` (stored in ~/.claude.json, keys inlined).
Project scope writes a gitignored .mcp.json. Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NODE_MODULES = ROOT / "mcp" / "node" / "node_modules" / "@cyanheads"
ATLAS_PY = ROOT / "mcp" / "atlas_server" / ".venv" / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
UV_TOOL_BIN = Path(os.environ.get("UV_TOOL_BIN_DIR", Path.home() / ".local" / "bin"))
EXE = ".exe" if os.name == "nt" else ""

# Keys each server can use. Every key is optional; see .env.example for where to get them.
KEYS = {
    "kif13a-atlas": ["NCBI_API_KEY", "NVIDIA_API_KEY"],
    "gnomad": ["NCBI_API_KEY"],
    "biomcp": ["NCBI_API_KEY", "ONCOKB_TOKEN"],
    "tooluniverse": [
        "NCBI_API_KEY", "SEMANTIC_SCHOLAR_API_KEY", "BIOGRID_ACCESS_KEY", "NVIDIA_API_KEY", "ESM_API_KEY", "OMIM_API_KEY",
    ],
    "zotero": ["ZOTERO_API_KEY", "ZOTERO_LIBRARY_ID"],
}


def load_env() -> dict[str, str]:
    values: dict[str, str] = {}
    env_file = ROOT / ".env"
    if env_file.exists():
        for line in env_file.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                v = v.strip().strip('"').strip("'")
                if v:
                    values[k.strip()] = v
    return values


def node_server(package: str, extra_env: dict[str, str] | None = None) -> tuple[dict, Path]:
    entry = NODE_MODULES / package / "dist" / "index.js"
    node = shutil.which("node") or "node"
    cfg = {
        "type": "stdio",
        "command": node,
        "args": [str(entry)],
        "env": {"MCP_TRANSPORT_TYPE": "stdio", "MCP_LOG_LEVEL": "warning", **(extra_env or {})},
    }
    return cfg, entry


def uv_tool(exe: str, args: list[str], extra_env: dict[str, str] | None = None) -> tuple[dict, Path]:
    path = UV_TOOL_BIN / f"{exe}{EXE}"
    cfg = {"type": "stdio", "command": str(path), "args": args,
           "env": {"PYTHONIOENCODING": "utf-8", **(extra_env or {})}}
    return cfg, path


def uv_tool_module(tool: str, module: str, args: list[str], extra_env: dict[str, str] | None = None) -> tuple[dict, Path]:
    """Start a uv tool through its own venv's python -m. Used when Windows Smart App Control
    blocks the generated .exe launcher (seen with zotero-mcp.exe, WinError 4551)."""
    tools_dir = Path(os.environ.get("UV_TOOL_DIR", Path(os.environ.get("APPDATA", Path.home())) / "uv" / "tools"))
    python = tools_dir / tool / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    cfg = {"type": "stdio", "command": str(python), "args": ["-m", module, *args],
           "env": {"PYTHONIOENCODING": "utf-8", **(extra_env or {})}}
    return cfg, python


def catalogue() -> dict[str, tuple[dict, Path | None]]:
    return {
        # project data layer (writes data/raw with provenance)
        "kif13a-atlas": (
            {"type": "stdio", "command": str(ATLAS_PY), "args": ["-m", "atlas_mcp.server"],
             "env": {"ATLAS_PROJECT_ROOT": str(ROOT), "PYTHONIOENCODING": "utf-8"}},
            ATLAS_PY,
        ),
        # variants and genetics
        "gnomad": node_server("gnomad-genetics-mcp-server", {"CANVAS_PROVIDER_TYPE": "duckdb"}),
        "biomcp": uv_tool("biomcp", ["serve"]),
        "opentargets": ({"type": "http", "url": "https://mcp.platform.opentargets.org/mcp"}, None),
        "ensembl": node_server("ensembl-mcp-server"),
        # protein, domains, structure
        "uniprot": node_server("uniprot-mcp-server"),
        "interpro": ({"type": "http", "url": "https://www.ebi.ac.uk/interpro/mcp"}, None),
        "protein": node_server("protein-mcp-server"),
        "pdbe": uv_tool("pdbe-mcp-server", ["--server-type", "pdbe_api_server"]),
        # interactions and broad fallback (GTEx, HPA, CADD, BioGRID, Europe PMC, OpenAlex, NIM ...)
        "string": ({"type": "http", "url": "https://mcp.string-db.org/mcp"}, None),
        "tooluniverse": uv_tool("tooluniverse", []),
        # literature
        "pubmed": ({"type": "http", "url": "https://pubmed.mcp.claude.com/mcp"}, None),
        "biorxiv": ({"type": "http", "url": "https://hcls.mcp.claude.com/biorxiv/mcp"}, None),
        # needs Zotero desktop running with "Allow other applications to communicate with Zotero"
        "zotero": uv_tool_module("zotero-mcp-server", "zotero_mcp.cli", ["serve"], {"ZOTERO_LOCAL": "true"}),
    }


def main() -> int:
    # Microsoft Store Python virtualizes %APPDATA%: it (and any child it spawns) cannot see uv tools
    # or uv's own interpreters under AppData\Roaming\uv, so every path check would fail.
    if "PythonSoftwareFoundation" in sys.executable or "WindowsApps" in sys.executable:
        print("Microsoft Store Python cannot see %APPDATA%\\uv. Run this instead:\n"
              f"  {ATLAS_PY} mcp\\setup_mcp.py {' '.join(sys.argv[1:])}")
        return 2

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="report only, do not write files")
    parser.add_argument("--scope", choices=["user", "project"], default="user")
    args = parser.parse_args()

    env = load_env()
    servers: dict[str, dict] = {}
    missing: list[str] = []
    for name, (cfg, path) in catalogue().items():
        ok = path is None or path.exists()
        keys_set = [k for k in KEYS.get(name, []) if k in env]
        print(f"  {'ok ' if ok else 'MISSING'} {name:<13} {path or cfg.get('url')}"
              + (f"   keys: {', '.join(keys_set)}" if keys_set else ""))
        if not ok:
            missing.append(name)
            continue
        if cfg["type"] == "stdio":
            cfg["env"].update({k: env[k] for k in keys_set})
        servers[name] = cfg

    if args.check:
        return 1 if missing else 0
    skipped = f"; skipped (not installed): {', '.join(missing)}" if missing else ""

    settings_path = ROOT / ".claude" / "settings.local.json"
    if args.scope == "project":
        (ROOT / ".mcp.json").write_text(json.dumps({"mcpServers": servers}, indent=2) + "\n", encoding="utf-8")
        settings_path.parent.mkdir(exist_ok=True)
        settings = json.loads(settings_path.read_text(encoding="utf-8")) if settings_path.exists() else {}
        settings["enableAllProjectMcpServers"] = True
        settings_path.write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
        print(f"\nwrote {ROOT / '.mcp.json'} with {len(servers)} servers{skipped}")
        return 0

    claude = shutil.which("claude")
    if not claude:
        print("claude CLI not found on PATH; install Claude Code or use --scope project")
        return 1
    failed = []
    for name, cfg in servers.items():
        subprocess.run([claude, "mcp", "remove", name, "--scope", "user"], capture_output=True, text=True)
        res = subprocess.run([claude, "mcp", "add-json", name, json.dumps(cfg), "--scope", "user"],
                             capture_output=True, text=True)
        if res.returncode != 0:
            failed.append(f"{name}: {(res.stderr or res.stdout).strip()[:200]}")
    # a leftover project .mcp.json would shadow the user entries in this folder
    (ROOT / ".mcp.json").unlink(missing_ok=True)
    if settings_path.exists():
        settings = json.loads(settings_path.read_text(encoding="utf-8"))
        settings.pop("enableAllProjectMcpServers", None)
        if settings:
            settings_path.write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
        else:
            settings_path.unlink()
    print(f"\nregistered {len(servers) - len(failed)} servers at user scope{skipped}")
    for line in failed:
        print("  FAILED", line)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
