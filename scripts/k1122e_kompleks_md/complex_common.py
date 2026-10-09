"""Shared helpers for the Rab10(GTP, Mg):KIF13A complex MD study.

Protocol loading, the 12-system manifest (4 arms x 3 starts), hashing and
provenance sidecars, in the same style as the monomer main arm (common.py).

    python complex_common.py --protocol complex_protocol.json            # print the manifest
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
DEFAULT_PROTOCOL = HERE / "complex_protocol.json"
ARM_ORDER = ["WT", "K1122E_AF3", "K1122E_paired", "K1122R_paired"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def write_text(path: Path, text: str) -> None:
    """LF line endings on every platform."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def load_protocol(path: Path = DEFAULT_PROTOCOL) -> dict[str, Any]:
    protocol = read_json(path)
    if protocol.get("schema_version") != 1:
        raise ValueError("Unsupported complex protocol schema")
    if list(protocol["arms"]) != ARM_ORDER:
        raise ValueError(f"Arm order must be {ARM_ORDER}")
    for arm, spec in protocol["arms"].items():
        if len(spec["starts"]) != 3 or len(set(spec["starts"])) != 3:
            raise ValueError(f"{arm}: three distinct starts required")
        for start in spec["starts"]:
            if start not in protocol["starts"]:
                raise ValueError(f"{arm}: unknown start {start}")
            mutation = spec["mutation"]
            native = protocol["starts"][start]["residue_1122"]
            if mutation and mutation["from"] != native:
                raise ValueError(f"{arm}: mutation {mutation} does not fit start {start} ({native})")
    kif = protocol["construct"]["kif13a"]
    if len(kif["sequence_wt"]) != kif["last"] - kif["first"] + 1:
        raise ValueError("KIF13A sequence length does not match its numbering")
    rab = protocol["construct"]["rab10"]
    if len(rab["sequence"]) != rab["last"] - rab["first"] + 1:
        raise ValueError("Rab10 sequence length does not match its numbering")
    stages = protocol["simulation"]["stages"]
    if stages[0]["kind"] != "minimize" or any(s["kind"] not in {"nvt", "npt"} for s in stages[1:]):
        raise ValueError("Stage list must start with one minimisation followed by MD stages")
    return protocol


def system_manifest(protocol: dict[str, Any], protocol_path: Path) -> list[dict[str, Any]]:
    """One run per (arm, start); the three starts of an arm are its replicas."""
    systems = []
    protocol_hash = sha256(protocol_path)
    for arm in ARM_ORDER:
        spec = protocol["arms"][arm]
        for replica, start in enumerate(spec["starts"], start=1):
            index = len(systems)
            systems.append({
                "index": index,
                "system_id": f"{arm}_{start}",
                "arm": arm,
                "start": start,
                "replica": replica,
                "pair_key": start if arm in {"WT", "K1122E_paired", "K1122R_paired"} else None,
                "af3_job": protocol["starts"][start]["af3_job"],
                "af3_model": protocol["starts"][start]["model"],
                "native_1122": protocol["starts"][start]["residue_1122"],
                "mutation": spec["mutation"],
                "residue_1122": spec["mutation"]["to"] if spec["mutation"] else protocol["starts"][start]["residue_1122"],
                "seed": protocol["seed_base"] + 100 * index + 1,
                "protocol_id": protocol["protocol_id"],
                "protocol_sha256": protocol_hash,
            })
    return systems


def system_by_id(protocol: dict[str, Any], protocol_path: Path, system_id: str) -> dict[str, Any]:
    for item in system_manifest(protocol, protocol_path):
        if item["system_id"] == system_id:
            return item
    raise KeyError(f"Unknown system id {system_id}")


def provenance(path: Path, inputs: list[Path], *, params: dict[str, Any],
               source: str = "Rab10:KIF13A complex MD") -> None:
    path = Path(path)
    record = {
        "source": source,
        "tool": "complex_md/scripts",
        "utc": utc_now(),
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
        "inputs": [{"path": str(Path(item).resolve()), "sha256": sha256(item)} for item in inputs if Path(item).exists()],
        "params": params,
    }
    write_json(Path(str(path) + ".prov.json"), record)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    args = parser.parse_args()
    print(json.dumps(system_manifest(load_protocol(args.protocol), args.protocol), indent=2))


if __name__ == "__main__":
    main()
