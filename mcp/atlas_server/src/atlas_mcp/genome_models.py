"""Genome language model scoring: Evo 2 40B zero-shot variant effect via NVIDIA NIM.

Score = log P(alt window) - log P(ref window) under Evo 2, with the variant centred in a GRCh38
window (8192 nt by default, as in the Evo 2 paper). Negative = alt is less likely = more disruptive.
The NIM `forward` endpoint returns the unembedding logits for every position (vocabulary = byte
values), so the sequence log-likelihood is computed locally.

Hosted calls spend NVIDIA credits: this tool is for targeted sets (max 50 variants per call).
Saturation scoring of whole genes belongs on TRUBA with local weights.
"""

from __future__ import annotations

import base64
import io
import re
import zipfile
from typing import Any

# Loaded at module import on purpose (Windows stdio deadlock if numpy is first imported inside a call).
import numpy as np

from . import core

NIM_FORWARD = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/forward"
ENSEMBL_REGION = "https://rest.ensembl.org/sequence/region/human"
VARIANT_RE = re.compile(r"^(?:chr)?([0-9XY]{1,2})[-:_ ](\d+)[-:_ ]([ACGT])[-:_>/ ]([ACGT])$", re.I)
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
MAX_VARIANTS = 50


def log_softmax(logits: np.ndarray) -> np.ndarray:
    lg = logits.reshape(-1, logits.shape[-1]).astype(np.float64)
    return lg - np.logaddexp.reduce(lg, axis=1, keepdims=True)


def token_logprobs(seq: str, logits: np.ndarray) -> np.ndarray:
    """out[j] = log P(seq[j] | seq[:j]) for j >= 1; out[0] = 0 (the first base has no prediction)."""
    lp = log_softmax(logits)
    ids = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    out = np.zeros(len(seq))
    out[1:] = lp[np.arange(len(seq) - 1), ids[1:]]
    return out


def _log_likelihood(seq: str, logits: np.ndarray) -> float:
    return float(token_logprobs(seq, logits).sum())


async def evo2_logits(seq: str) -> np.ndarray:
    """Unembedding logits (L x 512, vocabulary = byte values); row i predicts base i+1."""
    key = core.secret("NVIDIA_API_KEY")
    if not key:
        raise RuntimeError("NVIDIA_API_KEY is not set (add it to .env and rerun mcp/setup_mcp.py).")
    resp = await core.request(
        "POST", NIM_FORWARD, json_body={"sequence": seq, "output_layers": ["unembed"]},
        headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"}, timeout=600,
    )
    archive = zipfile.ZipFile(io.BytesIO(base64.b64decode(resp.json()["data"])))
    return np.load(io.BytesIO(archive.read("unembed.output.npy"))).reshape(len(seq), -1)


async def _evo2_ll(seq: str) -> float:
    return _log_likelihood(seq, await evo2_logits(seq))


async def evo2_score(
    variants: list[str], window: int = 8192, both_strands: bool = False, label: str = "evo2", save: bool = True
) -> dict[str, Any]:
    if len(variants) > MAX_VARIANTS:
        raise ValueError(f"At most {MAX_VARIANTS} variants per call (hosted credits); split the list or use TRUBA.")
    if not 128 <= window <= 8192:
        raise ValueError("window must be between 128 and 8192 nt")
    half = window // 2

    # Score = sum of per-base log-prob changes from the variant base onward. Positions before the variant
    # cannot change in a causal model, but the hosted NIM shifts them slightly when any later base changes
    # (measured 2026-09-15: p99 0.03, max 0.36, sum ~0.4 nats over 4096 upstream bases), so they are excluded.
    idx = half - 1
    ridx = window - 1 - idx
    ref_cache: dict[tuple[str, int], tuple[str, np.ndarray, np.ndarray | None]] = {}
    rows: list[dict[str, Any]] = []
    for raw in variants:
        m = VARIANT_RE.match(raw.strip())
        if not m:
            rows.append({"input": raw, "error": "use GRCh38 SNV format like 6-17764000-A-G"})
            continue
        chrom, pos_s, ref, alt = m.groups()
        pos, ref, alt = int(pos_s), ref.upper(), alt.upper()
        start = pos - half + 1  # variant sits at 0-based index half-1, window [start, start+window-1]
        key = (chrom, pos)
        if key not in ref_cache:
            url = f"{ENSEMBL_REGION}/{chrom}:{start}..{start + window - 1}:1"
            seq = (await core.request("GET", url, headers={"Accept": "text/plain"})).text.strip().upper()
            lp_fwd = token_logprobs(seq, await evo2_logits(seq))
            rc = seq[::-1].translate(COMPLEMENT)
            lp_rev = token_logprobs(rc, await evo2_logits(rc)) if both_strands else None
            ref_cache[key] = (seq, lp_fwd, lp_rev)
        seq, ref_fwd, ref_rev = ref_cache[key]
        if seq[idx] != ref:
            rows.append({"input": raw, "error": f"reference base at {chrom}:{pos} is {seq[idx]}, not {ref}"})
            continue
        alt_seq = seq[:idx] + alt + seq[idx + 1:]
        d_fwd = token_logprobs(alt_seq, await evo2_logits(alt_seq)) - ref_fwd
        row: dict[str, Any] = {
            "input": raw, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "window": window,
            "delta_ll_forward": round(float(d_fwd[idx:].sum()), 4),
            "site_term_forward": round(float(d_fwd[idx]), 4),
            "upstream_noise_forward": round(float(d_fwd[:idx].sum()), 4),
        }
        if both_strands:
            alt_rc = alt_seq[::-1].translate(COMPLEMENT)
            d_rev = token_logprobs(alt_rc, await evo2_logits(alt_rc)) - ref_rev
            row["delta_ll_reverse"] = round(float(d_rev[ridx:].sum()), 4)
            row["delta_ll_mean"] = round((row["delta_ll_forward"] + row["delta_ll_reverse"]) / 2, 4)
        rows.append(row)

    out: dict[str, Any] = {
        "model": "Evo 2 40B (NVIDIA NIM forward, unembed logits)",
        "score": "sum of log P(alt seq) - log P(ref seq) from the variant base onward; more negative = more disruptive. "
                 "Reference panel 2026-09-15: clinical/canonical splice positives -17 to -68, common intronic "
                 "variants -2.5 to +3.8; treat |score| < ~3 as indistinguishable from common variants.",
        "window": window, "both_strands": both_strands, "rows": rows,
    }
    if save and any("delta_ll_forward" in r for r in rows):
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)[:60]
        out["saved"] = core.save_tsv(
            "evo2", f"{safe}_w{window}_{core.today()}.tsv", rows, source="NVIDIA NIM arc/evo2-40b forward",
            urls=[NIM_FORWARD, ENSEMBL_REGION], tool="evo2_score",
            params={"window": window, "both_strands": both_strands},
        )
    return out
