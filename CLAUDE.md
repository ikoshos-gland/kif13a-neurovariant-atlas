# KIF13A NeuroVariant Atlas

Paralog-informed prioritisation of human KIF13A missense variants for neurological impact, calibrated on
KIF1A (KAND) clinical variants. 16-week plan, phases A-D (WP1-WP4). Conversation language: Turkish.

## Hard rules

- **Protein numbering is MANE Select.** Identifiers live in `config/targets.yaml` (verified 2026-09-15).
  - KIF13A: ENST00000259711.11 / NM_022113.6 / Q9H1H9, **1805 aa** (1770 aa is isoform Q9H1H9-2).
  - KIF1A: ENST00000498729.9 / NM_001244008.2 = **Q12756-3, 1791 aa**. UniProt canonical Q12756 is 1690 aa.
    Positions 1-393 agree; canonical 395-848 is +9 and 849-1690 is +101 in MANE.
    AlphaMissense per-protein tables and the canonical AlphaFold model use canonical numbering.
- **gnomAD variants: use the transcript query on the MANE transcript**, not the gene query. The gene query
  reports HGVSp on the most severe transcript (5 different transcripts for KIF13A missense).
- **ClinVar:** the classification field is `germline_classification`; the old `clinical_significance` is gone.
- **Pipeline data goes through the `kif13a-atlas` MCP server** (writes `data/raw` + provenance). Third-party
  servers are for exploration; do not build figures from their pasted output.
- Never cite a PDB ID, PMID or DOI that was not verified live. 8H3T is NOT KIF1A.
- **KIF1A calibration labels come from `data/interim/clinvar_curation/KIF1A_calibration_sets_*.tsv`**
  (see `docs/KIF1A_ClinVar_kurasyonu.md`). ClinVar B/LB missense has zero motor-domain variants, so the
  headline metric is within-motor AUROC against gnomAD motor controls (AC >= 3), CV grouped by residue;
  never report whole-protein P/LP vs B/LB AUROC without the "is motor domain" baseline (0.995).
  REVEL/MetaRNN/ClinPred are label-trained and stay out of the consensus.
- Heavy compute (AlphaFold3/OpenFold3, ESM saturation scoring, Evo 2 full runs, MD, CELLxGENE Census
  cell-level work) runs on TRUBA via SLURM, never on the login node.
- Never use the em dash character in any text or file.

## MCP servers

Registered at **user scope** (available in every new session, any folder) by
`mcp\atlas_server\.venv\Scripts\python.exe mcp\setup_mcp.py`
(not the Microsoft Store `python`, which cannot see `%APPDATA%\uv`).

| Need | Server | Notes |
|---|---|---|
| Project pulls with provenance | `kif13a-atlas` | gnomAD, ClinVar, VEP+plugins, AlphaMissense, CADD, Evo 2 40B (NIM, max 50 SNVs/call, spends credits), AlphaFold files, UniProt features, PDB coverage, paralog map, GTEx, HPA, CELLxGENE |
| gnomAD exploration, SQL over variant tables | `gnomad` | cyanheads, duckdb canvas |
| dbSNP, MyVariant/dbNSFP, GWAS Catalog v2 | `biomcp` | `get`/`search` tools; MyVariant gnomAD fields are v2.1.1 |
| L2G, credible sets, target evidence | `opentargets` | official hosted; do not use its constraint numbers |
| Homology, xrefs, region features | `ensembl` | single-variant VEP only, no plugins |
| UniProt entries, isoforms, ID mapping | `uniprot` | |
| Domains, InterProScan | `interpro` | official EBI hosted |
| Structure search, Foldseek, TM-align, AFDB | `protein` | |
| SIFTS, validation, interfaces | `pdbe` | official PDBe |
| Interactions, enrichment | `string` | official; `species` must be a string "9606" |
| Broad fallback (GTEx eQTL, HPA, BioGRID, Europe PMC, OpenAlex, S2, NIM) | `tooluniverse` | compact mode: `list_tools` / `grep_tools` / `execute_tool` |
| Literature | `pubmed`, `biorxiv`, `zotero` (+ user-scope `google-scholar`) | bioRxiv has no keyword search; use Europe PMC `SRC:PPR` |

Health check: `mcp\atlas_server\.venv\Scripts\python.exe mcp\smoke_test.py`

> **Interpreter status, 16 Sep 2026:** the uv-managed CPython 3.12.14 behind `mcp\atlas_server\.venv`
> is now blocked by Windows ("Permission denied" on the interpreter, "Unable to create process" via
> the venv shim), so the venv cannot launch at all. Verified working replacement: the system
> `python` (3.13.14), which already carries httpx, numpy, pyyaml, biopython, pandas, scipy and mcp,
> imports `atlas_mcp.core` and reads `.env` correctly. Run project scripts with plain `python` from
> the repo root. Every reproduction command in `HANDOFF.md` still names the dead venv, and the
> user-scope MCP registrations point at it too, so the `kif13a-atlas` server needs re-pointing
> before it will start again.

## Layout

`config/` verified IDs · `mcp/` servers, setup, smoke test · `data/raw` pulls (gitignored, provenance in
`data/raw/PROVENANCE.jsonl`) · `data/interim` derived maps · `workflow/` Nextflow/Snakemake (WP1) ·
`docs/MCP_KURULUM.md` setup guide and known traps.
