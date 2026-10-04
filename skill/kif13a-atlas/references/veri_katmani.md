# Data layer

Fourteen MCP servers are registered at user scope, so they load in any new Claude Code session, in any
folder. Regenerate the registration after adding an API key:

```powershell
mcp\atlas_server\.venv\Scripts\python.exe mcp\setup_mcp.py      # user scope (default)
mcp\atlas_server\.venv\Scripts\python.exe mcp\smoke_test.py     # every server, one real KIF13A query
claude mcp list
```

> **Broken as of 2026-09-16 and not yet repaired.** The uv-managed CPython 3.12.14 behind
> `mcp\atlas_server\.venv` is blocked by Windows ("Permission denied" on the interpreter, "Unable to
> create process" through the shim), so both commands above fail and the `kif13a-atlas` server will not
> start. The registrations still point at that dead interpreter. Project scripts run fine under the system
> `python` 3.13.14, which already carries httpx, numpy, pyyaml, biopython, pandas, scipy and mcp and
> imports `atlas_mcp.core`. Repair options, in order: rebuild the venv with a working interpreter, or
> re-run `setup_mcp.py` pointing the registrations at the system python. Note that `setup_mcp.py` refuses
> to run under Microsoft Store python by design, and the system python here **is** the Store build, so
> that guard has to be reconsidered rather than bypassed blindly.

Full setup notes and the per-source mapping live in the project's `docs/MCP_KURULUM.md`.

## The project server: kif13a-atlas

Written for this project (`mcp/atlas_server`), 16 tools. What makes it different from the third-party
servers is that every tool writes its result under `data/raw` with a `.prov.json` sidecar and a line in
`data/raw/PROVENANCE.jsonl` recording the source URL, parameters, UTC timestamp and sha256. Pipeline data
comes from here; third-party servers are for exploration.

| Tool | Use |
|---|---|
| `atlas_targets` | verified identifiers and the isoform traps from `config/targets.yaml` |
| `gnomad_constraint`, `gnomad_variants` | constraint, and MANE-numbered variants with frequencies and scores |
| `clinvar_variants` | complete ClinVar pull with review status, stars, traits, GRCh38 SPDI |
| `vep_annotate` | VEP on MANE with AlphaMissense, CADD, SpliceAI, REVEL, dbNSFP ESM1b and PrimateAI |
| `alphamissense` | per-protein table from AlphaFold DB, position filters, exact variant lookups |
| `cadd_region`, `evo2_score` | region CADD; Evo 2 40B delta log-likelihood (spends NVIDIA credits) |
| `alphafold_model`, `uniprot_features`, `pdb_structures` | structures, features, experimental coverage |
| `paralog_map` | MANE-to-MANE residue mapping with identity and local window identity |
| `gtex_expression`, `hpa_summary`, `cellxgene_expression` | bulk, protein-level and single-cell expression |
| `data_inventory` | what has already been pulled |

## Third-party servers

Variants and genetics: `gnomad` (cyanheads, duckdb canvas for SQL over variant tables), `biomcp`
(MyVariant and dbNSFP, dbSNP, GWAS Catalog), `opentargets` (official hosted; do not use its constraint
numbers, they disagree with gnomAD v4), `ensembl`.

Protein and structure: `uniprot`, `interpro` (official EBI hosted), `protein` (structure search, AlphaFold
DB, Foldseek, TM-align), `pdbe` (official).

Breadth and literature: `string` (official), `tooluniverse` (GTEx, HPA, BioGRID, Europe PMC, OpenAlex,
NVIDIA NIM; compact mode with `list_tools`, `grep_tools`, `execute_tool`), `pubmed`, `biorxiv` (no keyword
search upstream; search preprints through Europe PMC with `SRC:PPR`), `zotero` (needs Zotero desktop open),
and the pre-existing user-scope `google-scholar`.

## Sources reached without MCP, verified working (2026-09-16)

These were queried directly over HTTP during the verification round. None needs an account, a key or a
data access agreement. Record them here so nobody rediscovers them.

| Source | Route | What it gave |
|---|---|---|
| **DECIPHER patient pages** | `deciphergenomics.org/patient/<id>/genotype/table`, JavaScript rendered so a browser is needed | per-patient CNVs with size, gene content, inheritance, pathogenicity, array mean ratio, deposited variant and phenotype counts |
| **DECIPHER gene pages** | `deciphergenomics.org/gene/<SYMBOL>/overview/clinical-info` | %HI, pLI, LOEUF, sHet, pHaplo, pTriplo, pDN, pGOF, pLOF, and whether the gene appears in G2P, Morbid, GeneReviews, ClinGen and GenCC |
| **DDG2P gene list** | downloadable TSV, 2,869 records | allelic requirement and mechanism per gene-disease pair. The KIF1A mechanism split lives here |
| **gnomAD SV and exome CNV** | GraphQL, `gnomad_sv_r4` region query | structural variants with AC, AN, AF, filters and LOFTEE-style consequence. No subcohorts in v4 |
| **DGV Gold standard** | public track | control structural variants, the sanity check on any deletion claim |
| **dbVar** | E-utilities `db=dbvar` plus study pages | nstd101 is obsolete, nstd102 replaces it, both Case-Set with no control arm |
| **ClinGen dosage curation list** | `ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv` | haploinsufficiency scores, read live rather than hard-coded |
| **Europe PMC citations** | `/MED/<pmid>/citations` | who cited a case report and whether any citing paper is an independent case |
| **Ensembl overlap** | `/overlap/region/homo_sapiens/<chr>:<start>-<end>?feature=gene`, both GRCh38 and GRCh37 | exactly which genes a deletion interval contains, the check that makes a single-gene claim geometric |

## Sources used by the September and October axes, verified working

| Source | Route | Access | Used in |
|---|---|---|---|
| Kaplanis 2020, gene4denovo2, denovo-db, autism consortia | supplementary tables and downloads | open | de novo axis (9e) |
| Oxford BIG40 | `open.win.ox.ac.uk/ukbiobank/big40/` full statistics plus PheWeb | open, GRCh37 | imaging axis (9h), coloc (9j) |
| GWAS Catalog summary statistics | `ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST...` | open | UDIP-FA GCST90860937, sulcal depth GCST90095129, AD controls |
| MetaBrain cortex eQTL | `eqtlgen.molgeniscloud.org/api/v2/metabrain_cis_eqtl_v2?q=GeneSymbol==KIF13A` | open | imaging axis |
| Zeng 2026 single-nucleus brain eQTL | Synapse `syn61929919` (`cis_analysis` syn61930218, class and subclass files, coloc tables) | free Synapse account plus a personal access token with **download** scope; the user's token sits in `~/.synapseConfig` on TRUBA | coloc (9j) |
| Bryois 2022 cell-type eQTL | Zenodo record 7276971 | open | coloc (9j) |
| SingleBrain 2025 meta-analysis | Zenodo record 14908182; files are large, so `s06_singlebrain_stream.py` pulls only the KIF13A window with HTTP Range requests (a full-file job ran out of memory) | open | coloc (9j) |
| 1000 Genomes phase 3 | for LD panels (EUR and all) | open | coloc (9j) |
| CELLxGENE datasets | direct h5ad at `datasets.cellxgene.cziscience.com/<dataset_id>.h5ad`; ids in `scripts/truba/insan_yas/datasets.tsv` | open | lifespan (9i) |
| Lifespan PFC atlas extras | Synapse `syn52396927`: donor table with PMI and dreamlet results are OPEN (account needed), raw FASTQ needs a DUA | partly open | not yet |
| RGC-ME | Nature supplementary tables (open); browser `rgc-research.regeneron.com/me` needs a licence click | tables open | knockouts (9k) |
| Genes & Health | Nat Genet supplementary tables (open); the cloud bucket is 403 to anonymous access | tables open | knockouts (9k) |
| Chundru 2024 | Nat Genet supplementary tables and data file | open | knockouts (9k) |
| popEVE | supplementary tables plus `data.evemodel.org/popeve/v1.1/mutations/<NP>.csv` (KIF13A NP_071396.4, KIF1A NP_004312.2) | open | 9k |
| DDG2P | `ftp.ebi.ac.uk/.../2026_09_28/DDG2P_2026-09-28.csv.gz` | open | 9k |
| AlphaGenome Atlas / AVI | `alphagenome.google/atlas`, API | key under non-clinical terms, CC BY-NC 4.0; the user must obtain it | evaluated (9l), not used |

Identified as high value and **not yet used**: SpliceVault observed splice outcomes, the 100,000 Genomes
Project (Genomics England Research Environment, application), All of Us v8 controlled tier, and a future
UK Biobank imaging GWAS at 60,000 to 100,000 scans. No kinesin motor scan exists in MaveDB, ProteinGym or
the Human Domainome (checked 2026-10-03).

## TRUBA workspaces for this project

All under `/arf/scratch/avural/` (purged periodically; small results are mirrored locally):

| Folder | Axis | Size (2026-10-04) | Local mirror |
|---|---|---|---|
| `kif13a_allenbrain` | mouse Allen atlas (the user's original analysis) | about 103 GB cache | `C:\Users\pc\Claude\kif13a-allenbrain` (13 Sep results not mirrored) |
| `kif13a_eqtl` | cell-type eQTL coloc (9j) | 33 GB | `scripts/truba/eqtl/`, `data/interim/eqtl_hucre_tipi/` |
| `kif13a_insan_yas` | human lifespan (9i) | 38 GB | `scripts/truba/insan_yas/`, `data/interim/insan_yas/` |
| `kif13a_nakavt` | knockouts and popEVE (9k) | 0.3 GB | `scripts/truba/nakavt/`, `data/interim/nakavt/` |

Python on TRUBA: `source /arf/home/avural/miniconda3/etc/profile.d/conda.sh && conda activate allenbrain`
(anndata, pandas, scipy). Extra packages went to per-workspace folders (`Rlib/` with coloc and susieR,
`pylib/` with openpyxl and pypdf), never into home, to protect the inode quota. Internet only on the
`debug` partition. Operate through the `truba` skill.

## Keys

In the project's gitignored `.env`, inlined into the MCP registration by `setup_mcp.py`:

- `NCBI_API_KEY` set and verified (E-utilities 3 to 10 requests per second)
- `NVIDIA_API_KEY` set and verified (Evo 2 40B). Every call costs credits; ask before large runs.
- OpenAlex needs no key: keyless gives 1000 credits a day and a search costs 10.
- Synapse personal access token: created and placed by the user on TRUBA (`~/.synapseConfig`, mode 600,
  view plus download scope) on 2026-10-03. Never read or print it. An older view-only token named
  `sigma2026` was still active at that time; the user was advised to revoke it.
- Not obtained: BioGRID, ESM Forge, OMIM, DECIPHER client key, AlphaGenome API key.

## What is deliberately not MCP

Model training and heavy scoring are not data pulls: AlphaFold3 and OpenFold3 complexes, ESM saturation
scoring, full Evo 2 runs, molecular dynamics, cell-level single-cell work and colocalisation belong on
TRUBA via SLURM, never on a login node. See the workspace table above; the `truba` skill covers cluster
operation.
