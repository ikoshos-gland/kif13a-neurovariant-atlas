---
name: kif13a-atlas
description: Operate the KIF13A NeuroVariant Atlas project: paralog-calibrated missense prioritisation for an understudied gene, calibrated on its clinically characterised sibling, plus its evidence axes (structure, loss of function, de novo, regional constraint, white-matter imaging GWAS, cell-type eQTL colocalisation, human lifespan single-cell, biallelic knockouts, popEVE). Use it whenever work touches KIF13A or KIF1A variants, the folder C:\Users\pc\Claude\kif13a-neurovariant-atlas, its candidate list, calibration sets, AlphaMissense thresholds, any of those axes, its TRUBA workspaces (kif13a_eqtl, kif13a_insan_yas, kif13a_nakavt), its MCP data layer, or the NSD1/NSD3 multi-pair work. Trigger it too for "varyant onceliklendirme", "aday liste", "kalibrasyon", "paralog", "AlphaGenome", "popEVE", "kolokalizasyon", "nakavt", "ak madde", "hangi varyantlar zararli", or when the user asks to continue, explain or rebuild this project, even if they never name the gene.
---

# KIF13A NeuroVariant Atlas

A method and a live project. The method: when a gene has no clinical variant data, calibrate
variant-effect models on its closest paralog that does, then transfer the calibrated thresholds and
residue-level evidence. The project applies it to KIF13A (no known disease) using KIF1A (hundreds of
curated pathogenic variants) as the teacher.

Work in `C:\Users\pc\Claude\kif13a-neurovariant-atlas`. Read its `CLAUDE.md` before touching anything;
it holds the hard rules. This skill holds the reasoning, the state and the traps.

## Ground rules that keep results honest

- **Speak Turkish with the user.** Keep gene symbols, tool names and statistics in English. Never use
  the em dash character anywhere, in chat or in files.
- **Every number needs a source and a date.** Nothing enters a document or a message unless it came
  from a live query or a file on disk that you can point at. If a claim cannot be checked, mark it
  UNVERIFIED rather than smoothing it away.
- **Both classes of a comparison come from one source.** Scoring positives with VEP and negatives with
  gnomAD's own annotation produced a fake CADD AUROC of 0.914; one source gave 0.781. This is the
  single easiest way to fool yourself here.
- **Report the confound next to the result.** Whole-protein AUROC for this gene pair is meaningless:
  "is the variant in the motor domain" alone scores 0.995. Any headline number ships with its baseline.
- **Measure the denominator before reading a count as a signal.** "This lab has 1,293 ClinVar submissions
  and exactly one in KIF13A, so it is a deliberate candidate-gene call" was our own reasoning, and it was
  wrong: one record per gene is that lab's most common pattern, 60.6 percent of their genes. The rule
  applies equally to an absence, which means nothing until you compute the expected count. This is trap 26
  and it cost the project its most exciting-looking lead.
- **Audit your own best findings before the user gets attached to them.** When a result starts sounding
  like a discovery, hand it to an agent whose instructions are to break it, with the controls named in the
  prompt. Two such audits on 2026-09-16 dropped both leading human observations to the lowest evidence
  rung, and one of the arguments they broke was ours.
- **Protein positions are MANE Select.** The isoform traps in `references/tuzaklar.md` have already
  bitten this project twice.
- **Ask before long or costly runs.** Evo 2 calls spend the user's NVIDIA credits. Overnight background
  jobs die when the machine sleeps; say so rather than promising results by morning.
- **Never handle the user's credentials.** Synapse tokens, AlphaGenome API keys and similar are placed by
  the user themselves (for Synapse: `~/.synapseConfig` on TRUBA, mode 600). You may check that the file
  exists and print a token's *scopes*, never its value. If the user pastes a token into chat, do not use
  it; tell them to revoke it and place a new one themselves.
- **Heavy work runs on TRUBA through SLURM, never on the login node** (ban risk, see the `truba` skill).
  Only the `debug` partition has internet, so downloads run there as jobs.
- **Agents get interrupted by usage limits.** Three axis agents on 2026-10-03 were each stopped once or
  twice by rate limits; their TRUBA jobs kept running. Tell any agent to write its document early (verdict
  first) and to keep small results in the local project, so an interruption never loses the analysis.

## How to teach this project

The user is learning bioinformatics through it. Explain in small steps with concrete numbers from the
project rather than generic textbook prose, and put the conclusion first. Ask a check question only
when the user seems to want one. When they say "soru sorma", drop the Socratic questions and explain
continuously; quiz-style questions frustrate them. Keep `progress.md` in the project current.

## The method in one page

1. **Fix identity.** MANE transcript, protein length, UniProt accession for the target and the teacher.
2. **Build the variant universe.** gnomAD on the MANE transcript (not the gene query), ClinVar for both.
3. **Map residues.** Align the two MANE proteins; produce target position to teacher position with an
   identity flag. Report identity per domain, never as one global number.
4. **Curate labels on the teacher.** Positives are curated pathogenic missense with review stars,
   absent from gnomAD, with the ones whose classification leaned on in-silico evidence held out.
   Negatives must come from the same structural region as the positives, otherwise location alone
   separates the classes. Population variants from the target region work when ClinVar has none there.
5. **Measure with grouped bootstrap.** Resample residues, not variants. Report AUROC, its interval,
   sensitivity at fixed specificity, and the gene-specific threshold that specificity implies.
6. **Transfer.** Apply the calibrated threshold plus paralog evidence to the target, then add
   independent axes: structure, loss-of-function, cell type. Convergence across axes is the argument.
7. **State what would falsify the hypothesis** before anyone asks.

Details, including why each step exists, are in `references/yontem.md`.

## Where things stand

Read `references/proje_durumu.md` at the start of any session that continues the work: the "Read this
first" box, then section 12 (the scoreboard of every axis), then 6 and 6b. Short version, as of
2026-10-04:

- **The method works and has been validated out of sample.** AlphaMissense reaches 0.925 within the
  KIF1A motor domain against population controls (gene threshold 0.968 against a generic 0.564), and the
  same pipeline on NSD1/NSD3 gave 0.973. Candidate list and the 34,295 row saturation table exist.
- **Every disease hypothesis has been tested and none has been rescued.** Haploinsufficiency is largely
  out (0 de novo LoF in affected people across about 97,700 trios). Missense: the motor domain is only
  weakly constrained (o/e 0.814), the ATP pocket genuinely so (0.523, where 23 tier A candidates sit) but
  four times shallower than KIF1A (0.140); popEVE calls none of the 68 KIF1A-identical pathogenic changes
  severe in KIF13A. Recessive: no public cohort has the power to decide (even KIF1A is indistinguishable).
  Dominant negative is the only coherent mechanism left and it has never been measured.
- **The oligodendrocyte story got narrower.** Oligodendrocyte-lineage ranking replicates in human (OPC
  first in human, mature oligodendrocyte first in mouse). The mouse ageing increase does **not** replicate
  in 335 human donors (a real negative, power above 0.99). The white-matter GWAS haplotype does **not**
  colocalise with KIF13A eQTLs in oligodendrocytes or OPCs in two large single-nucleus eQTL atlases
  (PP.H3 0.95 to 0.9998); the only hint is an LD-panel-sensitive microglia signal (preliminary).
- **Six human observations exist and none carries alone.** Two were audited adversarially and fell to the
  lowest rung. Of the four loss-of-function events with known inheritance, three are maternal. Getting the
  inheritance of DECIPHER 451718 and HudsonAlpha A00166 is still the one action that could move a lead.
- **The project's real asset is the method, the dossier and the error-catching**, not a prediction. The
  honest pivot is the multi-pair method paper with KIF1A/KIF13A as the worked, mostly negative, case.

## Reference files

| File | Read it when |
|---|---|
| `references/proje_durumu.md` | starting any session; it is the state of the project and its numbers |
| `references/yontem.md` | applying or explaining the calibration method, or porting it to a new pair |
| `references/tuzaklar.md` | before writing any script that touches this data, and when something looks too good |
| `references/veri_katmani.md` | pulling data, fixing the MCP servers, adding an API key |
| `references/sonraki_adimlar.md` | deciding what to do next, or picking up an unfinished thread |

## Pipeline order

Everything is reproducible from the project scripts, in this order. Each writes provenance sidecars.

```powershell
python scripts\pull_wp1_core.py              # raw data freeze
python -u scripts\kif1a_clinvar_curation.py  # teacher labels
python -u scripts\kif1a_curation_extras.py   # motor negatives, literature
python scripts\kif1a_build_calibration_sets.py
python -u scripts\evaluate_calibration_sets.py
python scripts\build_kif13a_candidates.py    # the transfer
python scripts\build_saturation_table.py     # closes the population blind spot
python scripts\structure_axis.py             # structure axis
python -u scripts\lof_axis.py                # loss-of-function axis
python -u scripts\denovo_axis.py             # de novo axis (system python 3.13)
python -u scripts\regional_constraint.py     # regional missense constraint (ran under Anaconda)
python -u scripts\imaging_gwas_axis.py       # BIG40 white-matter imaging axis
python -u scripts\nsd_pair_calibration.py    # out-of-sample pair NSD1/NSD3
python scripts\paralog_pairs_scan.py         # multi-pair feasibility screen
```

The three axes added on 2026-10-03 ran on TRUBA, not locally. Their scripts (SLURM files included) are
mirrored under `scripts\truba\{eqtl,insan_yas,nakavt}\`, small results under
`data\interim\{eqtl_hucre_tipi,insan_yas,nakavt}\`, and the full workspaces live in
`/arf/scratch/avural/kif13a_{eqtl,insan_yas,nakavt}` (33 GB, 38 GB, 0.3 GB; scratch is purged
periodically, so anything not mirrored locally can disappear).

**Interpreter, changed 2026-09-16.** These used to run under
`mcp\atlas_server\.venv\Scripts\python.exe`. That venv's uv-managed CPython 3.12.14 is now blocked by
Windows and cannot launch at all. The system `python` 3.13.14 is the verified replacement: it carries
httpx, numpy, pyyaml, biopython, pandas, scipy and mcp, imports `atlas_mcp.core` and reads `.env`.
`HANDOFF.md` and the user-scope MCP registrations still name the dead venv and need re-pointing. See
traps 40 and 41.

## What good work looks like here

The valuable output of this project so far has not been a prediction. It has been catching the things
that would have made the predictions wrong: a protein length off by 35 residues, an isoform offset that
would have misplaced every KIF1A tail variant, a benign set that made the models look near-perfect, a
scoring pipeline that silently mixed two CADD versions, and a candidate universe that excluded exactly
the kind of variant the project is looking for. Later rounds added more of the same kind: a "16-fold"
expression gap that was a detection artefact, a Poisson p value of 0.049 that the empirical distribution
erased, a colocalisation that existed only under one LD panel, and a token that could read but not
download. When a result looks clean, spend the next few minutes trying to break it. That habit is the
deliverable.
