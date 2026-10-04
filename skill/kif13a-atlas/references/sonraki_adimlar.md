# Next steps and open decisions

State as of 2026-10-04. Read `references/proje_durumu.md` (the "Read this first" box and the section 12
scoreboard) first, and the last entries of `progress.md` in the project.

## Where the project stands

Every hypothesis the project was built on has now had its own test, and none was rescued:

- **Haploinsufficiency**: out at the severe end. Zero de novo LoF in affected people across about 97,700
  trios; less LoF-constrained than every ClinGen HI-3 neurodevelopmental gene; three of four inherited LoF
  events maternal. Mild or incompletely penetrant dose sensitivity (below 2.6-fold) survives.
- **Motor-domain missense**: weak. The ATP pocket is genuinely depleted (o/e 0.523) and holds 23 tier A
  candidates, but KIF1A's pocket is 0.140; popEVE calls none of the 68 KIF1A-identical pathogenic changes
  severe in KIF13A.
- **Recessive loss**: undecidable with public data; the same design cannot even separate KIF1A.
- **Oligodendrocyte / white matter**: expression in the lineage holds in mouse and human, but the mouse
  ageing increase is absent in human, and the white-matter GWAS haplotype is not a KIF13A eQTL in
  oligodendrocytes or OPCs (preliminary, 9j).
- **Dominant negative**: still the only coherent mechanism for tier A and still never measured.

The honest framing is now: a validated method (KIF1A/KIF13A plus NSD1/NSD3) with a thoroughly worked,
mostly negative case study and a six-observation human dossier. The pivot to the multi-pair method paper,
proposed on 2026-09-16 conditional on both axes landing negative, is due. Propose it plainly; do not
decide it for the user.

## Finish what is half done (highest value first)

1. **Write `docs/eqtl_hucre_tipi.md`** from the existing tables (`data/interim/eqtl_hucre_tipi/`,
   TRUBA log `logs/eqtl_pipe_6440189.out`). The agent was stopped by the user on 2026-10-04 after the
   pipeline finished and before the document. Verdict first ("against H_A in oligodendrocyte and OPC;
   fragile microglia hint"), then sources, validation and positive controls, the conditional table, the
   LD-panel sensitivity, power. The scratch workspace is 33 GB and will be purged eventually; mirror any
   missing scripts first.
2. **Test the microglia hint properly or retire it.** Zeng immune conditional signal rs73373199 gives
   PP.H4 0.55 to 0.62 only with the EUR LD panel and collapses with the all-1000G panel. The fix is
   in-sample LD or a method that does not need external LD (SuSiE-coloc with matched LD, or ask whether
   the Zeng authors release fine-mapped credible sets). If it survives it is the first positive link
   between the white-matter haplotype and KIF13A expression; if not, write it down as refuted.
3. **Sync `HANDOFF.md`** with sections 9i to 9l and the scoreboard. It stops at 2026-09-17.
4. **Append to `progress.md`** a 2026-10-03/04 entry covering the three TRUBA axes and the AVI evaluation.

## The actions only a human can take

- **Inheritance of the two unknown cases.** DECIPHER patient 451718 (146 kb intragenic deletion, autistic
  behaviour; Contact tab on the patient page) and ClinVar VCV002570672 / HudsonAlpha sample A00166
  (p.Lys1122Glu; Susan Hiatt, Greg Cooper, OrgID 505530). Drafts in `docs/vaka_arama.md` 9b and 9c.
- **RGC-ME browser lookup.** Accept the licence on `rgc-research.regeneron.com/me` and read the homozygote
  counts at 6-17771113-C-G and 6-17764129-C-T: are the two RGC "knockouts" the NMD-escaping last-intron
  variant c.4581+1G>C, as inferred?
- **GeneMatcher / Matchmaker Exchange.** A biallelic or dominant-negative KIF13A case can only be found by
  matchmaking. Claude drafts the submission, the user submits.
- **AlphaGenome API key**, if the splice axis is run with AVI: the user accepts the non-clinical terms and
  places the key; never paste it into chat.

## Next analyses, in order of value

1. **Dominant-negative test design.** The one mechanism left. In the lab: co-express wild-type and tier A
   mutants (ATP-pocket candidates such as K108T, C97R, G120V) and read motility or endosome distribution;
   expansion microscopy can read Rab11 endosome positioning. Computationally: Badonyi 2024 mechanism
   scores (PMID 39172982) and interface location of tier A, as a prior.
2. **Splice axis.** 26 curated KIF1A canonical splice positives, never scored. Frequency-matched intronic
   negatives; compare SpliceAI, Pangolin, AlphaGenome AVI and Evo 2 (Evo 2 separated canonical splice
   disruptions -17 to -68 from common intronic -2.5 to +3.8 in the control panel). Then K1122E, which sits
   at the third base of exon 28.
3. **AVI on the KIF1A calibration set.** Within-motor AUROC against the 0.925 of AlphaMissense; expect
   little gain because AVI's protein signal is AlphaMissense. AVI is frequency-trained, not label-trained,
   so it may enter the table.
4. **A pair the screen predicts to fail**, through the NSD pipeline. Without this false-positive control
   the "cheap metrics predict which pairs work" claim is supported only at the top of the ranking.
5. **Lifespan atlas follow-up.** PMI-corrected ageing model and the authors' own dreamlet KIF13A rows,
   both OPEN on Synapse `syn52396927` (the user's Synapse account now exists). Low priority: the negative
   is already well powered.
6. **Cell-type document.** The mouse plus human cell-type evidence (sections 9 and 9i) has never been
   written as one document; the human lifespan document now covers most of it.

## Datasets identified and still unused

- **100,000 Genomes Project gene burden** (Cipriani 2025, PMID 40011789): 34,851 rare disease cases with
  families, whole genomes, so intragenic deletions in the same patients. Genomics England Research
  Environment, application needed. The one public-adjacent source that could move the de novo or CNV
  question.
- **All of Us v8** (414,830 srWGS, 97,061 SV calls, EHR): carrier PheWAS for speech delay, autism,
  white-matter codes. Institutional controlled-tier access.
- **UK Biobank imaging at 100,000 scans** (milestone July 2025): the white-matter haplotype needs about
  46,000 imaged people; the next imaging GWAS release will cross that. No summary statistics yet.
- **Liu 2026 bioRxiv** (45,755 UKB WGS by 4,016 IDPs, rare-variant collapsing incl. white-matter tracts):
  no public browser found; KIF13A not mentioned.
- **MaveDB / ProteinGym / Human Domainome**: no kinesin motor scan exists (checked 2026-10-03). A motor
  domain DMS would be the decisive experiment, not a download.

Already queried and done, do not redo: Genebass, AstraZeneca PheWAS, All of Us public browser, Open
Targets, GWAS Catalog, BIG40 (9h), Kaplanis and the other de novo catalogues (9e), Bryois 2022, Zeng 2026,
SingleBrain 2025 (9j), CELLxGENE human brain sets and the lifespan atlas counts (9i), Genes & Health,
RGC-ME supplementary tables, Chundru 2024, popEVE (9k), DDG2P 2026-09-28.

## Smaller open threads

- **Weighted consensus.** Equal weights hurt (0.911 against 0.925). Fit weights only inside residue-grouped
  cross-validation and keep label-trained models out.
- **Interpreter and MCP repair.** The venv is dead (trap 41); `HANDOFF.md` and the user-scope MCP
  registrations still point at it.
- **Packaging.** Nextflow or Snakemake, DVC, public GitHub, Zenodo DOI: none exist, and **the repository
  still has no commits**. Offer the first commit; follow the user's commit-message rule (one short plain
  line, no co-author trailer).
- **Allen 13 Sep analyses** (cargo panel, habenula scan, kinesin correlation) exist only on TRUBA scratch;
  copy them into `kif13a-allenbrain` before a purge.

## Corrections to the original 16-week plan

Keep these visible when anyone works from the plan document:

- KIF13A is 1805 aa, so the saturation matrix is 34,295 substitutions, not 33,630.
- 8H3T is not a KIF1A structure. Use 8UTN and its siblings.
- KIF1A.org has no public variant table and DECIPHER needs a client key for deposited data, though patient
  and gene pages are readable without one. The KIF1A labels come from ClinVar plus functional literature.
- The success criterion "AUROC >= 0.85 on KIF1A" is trivially met by variant location alone. The metric
  that means something is within-motor AUROC against population controls.
- Evo 2 is not a missense model here. It belongs to the splice and noncoding part of the plan.
- The plan assumed KIF13A missense disease. The deliverable is conditional: "if a KIF13A missense disorder
  exists these are the candidates, and here is the evidence for and against that premise". The method
  paper does not depend on the premise being true.

## How to spend a session

If the user asks for progress rather than a specific task: write the missing eQTL document (item 1 above),
then sync `HANDOFF.md` and `progress.md`, then put the pivot to the method paper to the user as a decision.
If they want to keep hunting a disease link, the order is the microglia test, the splice axis, and the
dominant-negative design. Always read the last entries of `progress.md` before starting, and append to it
when the session ends.
