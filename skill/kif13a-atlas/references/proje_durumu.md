# Project state

Last updated 2026-10-04, after the three TRUBA axis agents (sections 9i to 9k) and the evaluation of
AlphaGenome Atlas (9l). Sections 1 to 9d are the 2026-09-16 state and still hold; 9e to 9h fold in what
`HANDOFF.md` recorded on 16 and 17 September. Everything below was measured in this project, not read
from a paper, unless a PMID or DOI is given. Re-verify anything older than a gnomAD or ClinVar release
before publishing it. Section 12 is the one-table scoreboard.

> **Read this first (2026-10-04).** Every disease hypothesis has now been put to a dedicated test and
> none was rescued: haploinsufficiency fell on the de novo axis, missense is weak outside the ATP pocket
> and popEVE does not call KIF13A's KIF1A-identical changes severe, the recessive question is undecidable
> with public data, the mouse ageing increase does not replicate in human, and the white-matter GWAS
> haplotype does not colocalise with KIF13A expression in oligodendrocytes or OPCs. What stands:
> oligodendrocyte-lineage expression in two species, a real but shallow ATP-pocket constraint, an untested
> dominant-negative mechanism, and a method that worked out of sample (NSD1/NSD3).

> **Read this too (2026-09-16).** On 2026-09-16 the two human observations the project leaned on were audited by
> agents instructed to break them, and both fell to the lowest evidence rung. One of the arguments that
> fell was **our own reasoning**, not someone else's data: "HudsonAlpha has 1,293 submissions and exactly
> one in KIF13A, so this is a deliberate candidate-gene call" collapsed once the submitter's baseline was
> measured (one record per gene is their most common pattern, 60.6 percent). The lesson generalises: a
> ratio is not evidence until its denominator has been measured. Section 6b carries the full audit.
>
> The current honest headline is: **KIF13A is a candidate gene with a collectable dossier of six human
> observations and no established disease, and the loss-of-function reading is weaker than it looked
> because three of the four events with known inheritance are maternal.**

## 1. Identities (verified live 2026-09-15 against Ensembl and UniProt)

| Gene | MANE transcript | RefSeq | UniProt (MANE isoform) | Length | Locus |
|---|---|---|---|---|---|
| KIF13A | ENST00000259711.11 | NM_022113.6 | Q9H1H9 | 1805 aa | chr6:17758198-17988274, minus |
| KIF1A | ENST00000498729.9 | NM_001244008.2 | Q12756-3 | 1791 aa | chr2:240713761-240824320, minus |
| KIF13B | ENST00000524189.6 | NM_015254.4 | Q9NQT8 | 1826 aa | chr8 |
| KIF1B | ENST00000676179.1 | NM_001365951.3 | O60333 | 1816 aa | chr1 |
| KIF1C | ENST00000320785.10 | NM_006612.6 | O43896 | 1103 aa | chr17 |

KIF1A canonical UniProt Q12756 is 1690 aa and is NOT the clinical numbering. Offsets from canonical to
MANE: 1-393 identical, 395-848 add 9, 849-1690 add 101.

KIF13A domains in MANE numbering: motor 5-352, neck linker 353-380, CC1 381-468, FHA 469-519, tail 520-1805.
KIF1A: motor 5-354, neck 355-365, neck coil 366-383, FHA 525-581 (MANE), PH 1676-1774.

## 2. Constraint (gnomAD v4.1, 2026-09-15)

| | LOEUF | pLI | missense Z | observed/expected missense |
|---|---|---|---|---|
| KIF13A | 0.449 | 0.9998 | 2.79 | 0.86 (1901/2209) |
| KIF1A | 0.259 | 1.0000 | 6.87 | 0.68 (1683/2486) |

Both are loss-of-function intolerant. KIF1A is under severe missense selection; KIF13A largely tolerates
missense. This asymmetry is the strongest argument against the project's original premise.

**When the user asks "KIF1A is better studied, isn't more missense expected?"** (asked 2026-10-03): no.
Constraint comes from gnomAD population sequencing, which reads every gene in the same people at the same
depth, so study effort cannot enter it. The direction is in fact the opposite: the population carries
fewer KIF1A missense than expected (1,683 of 2,486) than KIF13A (1,901 of 2,209). The control is the
synonymous o/e, 0.98 for KIF13A and 1.01 for KIF1A: if coverage or the model differed between the genes,
synonymous would deviate too. Missense Z behaves as (expected minus observed) divided by the square root
of expected, scaled by one genome-wide constant (2.345 for both genes, checked by hand). Z grows with gene
length, so o/e is the effect size: approximate 90 percent Poisson intervals 0.83 to 0.89 against 0.65 to
0.70, non-overlapping. Study effort **does** bias ClinVar (135 P/LP missense in KIF1A against 0 in
KIF13A), so "no pathogenic KIF13A variant in ClinVar" is a gap, never evidence of absence.

## 3. Variant universes

- KIF13A gnomAD v4.1 on the MANE transcript: 10,936 variants, of which 2,671 missense. Confidently
  observed missense (AC >= 1 and no AC0 filter): 2,075. AC0-filtered: 596. The gene-level query returns
  11,728 and numbers HGVSp on the most severe transcript, so it is the wrong query for this project.
- KIF13A ClinVar: 313 records, 293 single-gene, 217 VUS, 0 pathogenic or likely pathogenic missense.
  Exactly one record carries a real phenotype (see section 6).
- KIF1A ClinVar: 3,372 single-gene records; 135 P/LP missense, 90 B/LB missense, 980 VUS, 213 conflicting.
- All possible KIF13A substitutions: 1805 x 19 = 34,295 (AlphaMissense covers all of them). Of these,
  10,621 are reachable by a single nucleotide change; gnomAD has seen 2,615 distinct protein changes
  (24.6 percent of the reachable set).

## 4. Calibration sets (`data/interim/clinvar_curation/KIF1A_calibration_sets_20260915.tsv`)

| Set | n | Role |
|---|---|---|
| missense_positive_core | 111 (110 motor, 70 residues) | positives; 98 sit at residues identical in KIF13A |
| missense_negative_motor_population | 61 | negatives, gnomAD motor missense AC >= 3 |
| missense_negative_motor_population_ac2 | 32 | sensitivity |
| missense_negative_clinvar_outside_motor | 88 | ClinVar B/LB, none in the motor domain |
| splice_positive | 26 | canonical splice P/LP, untested so far |
| heldout_pp3_dependent_positive | 18 | classifications that leaned on in-silico evidence |
| sensitivity_conflict_P_vs_VUS | 28 | sensitivity only |
| reference_lof | 33 | nonsense and frameshift, not missense labels |

A255V was rescued into the core set: its only ClinVar record is a 0-star OMIM entry, but two functional
papers characterise it. 14 core positives appear in functional literature.

## 5. Measured performance (motor domain, 110 positives vs 61 gnomAD controls, residue-grouped bootstrap)

| Model | AUROC | 95 percent interval | Sensitivity at 90 percent specificity | Threshold |
|---|---|---|---|---|
| MetaRNN (label-trained, reference only) | 0.949 | 0.913-0.977 | 0.85 | 0.928 |
| AlphaMissense | 0.925 | 0.878-0.964 | 0.79 | 0.968 |
| REVEL (label-trained, reference only) | 0.920 | 0.868-0.959 | 0.71 | 0.855 |
| ESM1b | 0.851 | 0.774-0.917 | 0.45 | -14.4 |
| PrimateAI | 0.801 | 0.720-0.875 | 0.37 | 0.883 |
| CADD v1.7 | 0.781 | 0.696-0.858 | 0.55 | 31 |
| Evo 2 40B | 0.742 | 0.643-0.839 | 0.32 | -22.7 |
| Baseline: identical residue in KIF13A | 0.744 | | | |
| Baseline: "is motor domain" (whole protein) | 0.995 | | | |

Equal-weight consensus of the four non-label-trained models: 0.911, below AlphaMissense alone.
AlphaMissense plus Evo 2: 0.887. Evo 2 does not help even where AlphaMissense is uncertain (0.611 in the
0.50-0.99 band). Model agreement: AlphaMissense with MetaRNN 0.80, with REVEL 0.71, with CADD 0.49.

The practical consequence: AlphaMissense's generic cutoff 0.564 would call 982 of KIF13A's 2,671 observed
missense variants damaging. The gene-calibrated 0.968 calls 295.

## 6. What is actually known about KIF13A in humans

Six independent human observations exist. Section 6b records what happened when the two that looked
strongest were verified adversarially on 2026-09-16; read both sections together, because the headline
changed.

| Observation | Source | Inheritance | Classification |
|---|---|---|---|
| p.Lys1122Glu, tail, apraxia + seizure + moderate ID + Chiari | ClinVar VCV002570672 | **unknown** | VUS, research |
| 253 kb deletion, whole gene, speech and language delay | ClinVar VCV000442594 | **unknown** | VUS, no criteria |
| 146 kb deletion, fully intragenic, autistic behaviour | DECIPHER 451718 | **unknown** | Uncertain |
| 182 kb deletion, 73 percent of the gene, no phenotype recorded | DECIPHER 504413 | **maternal** | **likely benign** |
| c.2357dupA mosaic frameshift, hypomelanosis of Ito + delay | Lam 2021 PMID 33526817 | de novo, postzygotic | published case |
| 2 truncating variants (stop_gained, frameshift) | DDD research data | **both maternal** | both Uncertain |

**The pattern that matters: of the four loss-of-function events whose inheritance is known, three are
maternally inherited and none is constitutional de novo.** The only de novo event is Lam's, and it is
mosaic, which is not evidence for constitutional haploinsufficiency. For a dominant developmental
disorder the expectation is the opposite.

Institutional status, re-checked live 2026-09-16: KIF13A is absent from ClinGen gene-disease validity,
ClinGen dosage sensitivity (1,540 genes curated, KIF13A not among them), **DDG2P (2,869 records, 12
kinesins present, KIF13A absent)**, GenCC, Morbid, GeneReviews, PanelApp UK and Australia, and the OMIM
phenotype map. Present only in SysNDD as "Limited", resting on Lam 2021 alone.

**Paralog calibration that reframes the whole loss-of-function question.** DDG2P carries two KIF1A
entries: G2P00959 is **monoallelic_autosomal** with mechanism **"altered gene product structure"**, and
G2P00963 is **biallelic_autosomal** with **"absent gene product / loss of function"**. So even in the
teacher gene, losing one copy is not how disease happens; one copy causes disease only by producing a
structurally altered protein. Assuming KIF13A haploinsufficiency therefore assumes a mechanism that does
not operate in its closest disease-bearing sibling.

Corrections worth remembering:
- The ClinVar records for A254V and R305L are TCGA tumour splicing annotations from the Rogan lab, not
  patients.
- ClinVar holds **6** KIF13A copy-number **losses**, not 20 (the 20 counted losses plus gains). Its
  pathogenic losses span **23 to 382** genes, not 23 to 237.
- Of 313 KIF13A ClinVar records, 221 are VUS, 15 are Pathogenic and **all 15 are megabase-scale 6p CNVs**;
  there are **0 Likely pathogenic** records and **0 pathogenic intragenic variants of any class**.

## 6b. Verification round, 2026-09-16 (two adversarial audits)

Both leads were handed to agents instructed to find reasons the claim was wrong. Reports:
`docs/dogrulama_delesyon.md` (53.8 KB) and `docs/dogrulama_hudsonalpha.md` (577 lines).
**Both leads dropped to "suggestive", the lowest rung.** This is the most important state change since
the axes were measured.

### The 253 kb deletion (VCV000442594)

Everything in the record checked out technically (253,376 bp, GRCh37 chr6:17,758,616-18,011,991, ISCA
site 1, affected child, HP:0000750). The evidence value collapsed on six points.

| Finding | Value |
|---|---|
| gnomAD v4 exome CNV `436229__DEL` | 448,271 bp, **99.2 percent of KIF13A**, SC 1/464,297, SF 2.15e-06 |
| gnomAD SV `DEL_chr6_33713a79` | 164,000 bp, 71.3 percent of the gene, PASS, "lof", AF 8e-06, carrier afr XX |
| Total screened | about 527,343 individuals, at least 3 independent gene-destroying deletions, roughly 1 in 175,000 |
| Inheritance | `Origin: unknown`; age, sex, family, co-occurring variants, number tested all absent |
| Last evaluated | **2014-05-05**, before gnomAD SV existed |
| Submitter comment | a template identical across every nstd37 submission |
| Cited support | Kaminsky 2011 PMID 21844811 analysed only 14 recurrent CNV regions; **6p22.3 is not one of them**, so the citation carries zero case-control statistics for this locus |
| Regulatory content | 37 enhancers, 9 CTCF sites, 1 promoter, 1 open chromatin region |
| Nearest neighbours | NUP153 at 52.0 kb, NHLRC1 at 108.7 kb. **NUP153 is more constrained than KIF13A** (LOEUF 0.277 versus 0.449) and more haploinsufficiency-prone by DECIPHER %HI (24.4 versus 42.3) |
| DGV Gold standard | 3 overlapping records, largest loss only 2,386 bp, no whole-gene deletion (the one finding in the claim's favour) |

**Power caveat that must travel with the gnomAD finding.** gnomAD v4 has **no subcohorts**: the
`StructuralVariantDatasetId` enum offers `_controls` and `_non_neuro` only for v2.1. So the carriers of
those deletions cannot be shown to be neurologically normal, only to be outside a severe paediatric
disease ascertainment. In v2.1 controls the expected count for a 1-in-175,000 event is **0.03 to 0.06**,
so observing zero there is a power failure, not evidence.

**dbVar is not usable as a case-control source here.** nstd101 was made obsolete on 2021-11-03 and
returns nothing through E-utilities; its replacement nstd102 is `Study Type: Case-Set` with no control
arm, and all 21 KIF13A calls have empty `sample_id`, `validation` and `zygosity`.

### DECIPHER, which both adds and subtracts

- **Patient 451718** is the cleanest lesion in the whole project: 6:17,825,889-17,972,163 GRCh38,
  146,275 bp, **entirely inside the gene** (KIF13A spans 17,758,198-17,988,274, so single-gene status is
  geometric, not a database assertion), 63.6 percent of the gene removed, array mean ratio **-0.999**,
  heterozygous, Uncertain, **inheritance unknown**, phenotype "Autistic behavior". Exactly 1 variant and
  1 phenotype deposited. Verified independently on Ensembl GRCh38 and GRCh37: the only other features in
  the interval are one lncRNA and one pseudogene.
- **Patient 504413** is the counterexample: 182,672 bp, **73.0 percent of the gene**, **maternal**, called
  **likely benign** by the depositing centre, no phenotype recorded.
- DECIPHER population CNVs: the largest KIF13A deletion in the control set is **11,902 bp** (1 in 5,919).
  No gene-destroying deletion.
- DECIPHER gene metrics: %HI **42.3** (threshold 10), pLI 0.99, LOEUF 0.47, **sHet 0.049**, **pHaplo
  0.99**, pTriplo 0.60, missense Z 2.27; pDN 0.560, pGOF 0.598, pLOF 0.616 are non-discriminating.
  **Haploinsufficiency predictors disagree with each other** (pHaplo for, sHet and %HI against). Report
  the disagreement, never one side of it.

### p.Lys1122Glu (VCV002570672)

The reading that made this record look special was **wrong, and the error was ours**.

| Check | Result |
|---|---|
| Submission context | `StudyName = CSER-HudsonAlpha`, a batch of 384 variants of which 124 are VUS |
| "Exactly 1 KIF13A submission means a deliberate call" | **False.** 368 of their 607 single-gene genes (60.6 percent) have exactly 1 record. One record per gene is their **most common** pattern |
| Submitter scale | 1,293 SCVs from 1,118 probands, all method "research", 79.5 percent carrying HPO terms |
| De novo flag | 152 variants in the same batch are marked de novo. **Ours is `Origin: unknown`** |
| Affected status | `AffectedStatus = unknown`. Our `docs/fonksiyon_kaybi_ekseni.md` section 4.4 said "present"; that came from `ClinicalFeaturesAffectedStatus`, a different field meaning only that the listed HPO terms are present |
| Submitter comment | empty (43 records in the same batch carry one) |
| In their publications | KIF13A appears **0 times** across 3 full texts covering about 3,000 probands (Bowling 2017 PMID 28554332 states VUS are systematically shared via ClinVar) |
| gnomAD v4.1 | `6-17785639-T-C` not found, AC 0. Region well covered (92 variants in 40 bp, AN about 1.45M). **Not distinctive**: about 82 percent of KIF13A's 11,918 possible missense SNVs are absent from gnomAD |
| Residue conservation | **Lys in 186 of 203 orthologs**, genuinely conserved. This is the record's one real strength |
| Local context | Sits in coiled coil 1100-1138, which carries 44 gnomAD missense with AC > 0, including **p.Lys1104Glu at AC 1, the identical Lys to Glu change** |
| Scores | AlphaMissense **0.8044** confirmed, 0.164 below the calibrated 0.968 cutoff. CADD 25.2 and PolyPhen 0.996 D, but MPC 0.573, REVEL 0.642, and PrimateAI, FATHMM and DEOGEN2 call it tolerated |
| Annotation subtlety | VEP calls it `missense + splice_region` (third base of exon 28), so its phyloP 8.0 may partly reflect splice constraint rather than protein constraint |
| Chiari malformation | Base rate 3.6 percent among 14,116 children imaged for any indication (PMID 21806364); in a developmental delay cohort **6 of 6 were unrelated to aetiology** (PMID 36548185). No kinesin gene causes Chiari. One of HudsonAlpha's own 3 Chiari records carries a Likely pathogenic KMT2C frameshift |

**The single action that would move this record from "suggestive" to "moderate":** ask the submitter
(Susan Hiatt, Greg Cooper, OrgID 505530) whether trio data exists for sample A00166. A de novo call
changes its weight by a full rung. A human has to send that request.

### Lam 2021, verified live

Abstract fetched from Europe PMC 2026-09-16: **exome and genome** sequencing were both done, the variant
is de novo and postzygotic, and mosaicism was inferred from a low mutant allele ratio. Two things weigh
against it: the paper has been **cited 3 times in 5 years and none of the citing papers reports an
independent KIF13A case**, and its own supporting argument is the 6p22.3-p23 microdeletion syndrome,
which is the contiguous-gene inference the project is trying to avoid.

## 7. Candidate list (`data/processed/KIF13A_candidates_20260915.tsv`)

295 of 2,671 observed missense pass the calibrated cutoff. Tiers (confidently observed in brackets):
A 52 (26), B 73 (43), C 170 (110), D 33 (23).

- A: cutoff passed, residue identical in KIF13A, KIF1A pathogenic at that residue
- B: cutoff passed, identical residue in the motor domain, no KIF1A disease variant there
- C: cutoff passed elsewhere (mostly tail, no paralog support)
- D: below cutoff but paralog-supported; the group where the model may be wrong

Strongest confident candidates: R348G (KIF1A R350, functionally characterised), R258S (R254), P303R
(P305, cryo-EM characterised), E257K (E253), S278L (S274), R173G (R169), D306H, G283V, H222R, A254V.
Only 4 of the 52 tier A variants appear in ClinVar at all, 2 of them as VUS, none as pathogenic.

## 7b. Saturation lookup table (2026-09-16, closes the blind spot)

`data/processed/KIF13A_saturation_lookup_20260916.tsv`, all 34,295 possible substitutions, built because
89 percent of KIF1A's pathogenic missense are absent from gnomAD and a population-only candidate universe
therefore excludes what the project is hunting. Of these, 10,627 are reachable by a single nucleotide
change, 8,038 pass the calibrated cutoff (1,499 of them SNV-reachable), 31,680 have never been seen in
gnomAD. Tiers over the whole space: A 929 (234 SNV-reachable), B 1,944, C 5,165, D 211.

The number worth remembering: **68 reachable substitutions are the exact same amino acid change that is
pathogenic in KIF1A**; 46 of those have never been seen in a human. 52 are tier A, of which 35 are
unobserved. Examples: K108T in the ATP pocket (AM 0.997), C97R, C97Y, G120V, A90D, Y94D, G83R.

182 of the 234 SNV-reachable tier A rows were absent from the old population-based list, and of the 1,499
reachable substitutions over the calibrated cutoff only 176 have a confident gnomAD carrier. The
population-only universe was missing most of the evidence-backed candidates.

All nine sanity checks passed, including reproduction of the existing 2,671 candidate tiers (52/52 tier A)
and K1122E's presence with a `c.` notation matching its ClinVar record, which independently verifies the
CDS. Structure columns are filled only for positions 1-380 and carry a `structure_assessed` flag.
Document: `docs/missense_sonuc.md`; figure: `figures/KIF13A_missense_saturation_blindspot.svg`.

## 8. Structure axis (`docs/yapi_ekseni.md`, 2026-09-16)

KIF13A has no experimental structure. The AlphaFold model was superposed onto KIF1A cryo-EM structure
8UTN (KIF1A[1-393] plus tubulin, AMP-PNP, two heads bound): 1.03 A RMSD over 311 trimmed CA pairs.
Tier A distributes as 23 ATP pocket, 18 microtubule interface, 9 buried core, 1 surface, 1 neck linker,
so 79 percent sit at a functional site. Tier B is the opposite (mostly buried core), which reads as a
folding-disruption class rather than a weaker version of tier A.

Rigid-body placement inflates contact counts and produces clashes, so "contact yes or no" is reliable and
the exact distance is not. R350 contacts tubulin only in one of the two heads (2.4 A in head N, 7.4 A in
head K), so tubulin distance is the minimum over both heads.

KIF1A structures: 8UTN-8UTS wild type, 8UTT-8UTW P305L, 8UTY P364L, 9YA5 and 9YA7 R350G, 9YAB and 9YAI
R350W. 8UTX and 8H3T are NOT KIF1A.

## 9. Cell-type axis

Human CELLxGENE brain data: KIF13A peaks in the oligodendrocyte lineage (OPC first), oligodendrocyte
lineage mean score 1.49 against 0.16 for neuron types. KIF1A is the mirror image: brainstem motor neuron
first, 0.15 in oligodendrocyte lineage. KIF13B also glial.

Mouse, from the user's own TRUBA analysis at `C:\Users\pc\Claude\kif13a-allenbrain` (public repo
`ikoshos-gland/kif13a-allenbrain`; TRUBA `/arf/scratch/avural/kif13a_allenbrain`, read-only for other
work; Allen Brain Cell Atlas release 20260711, 2,341,350 reference cells plus 1,162,565 cells from 108
ageing mice):

- Mature oligodendrocyte ranks first of 265 subclasses, OPC second. **The "sixteen-fold above astrocytes"
  figure is wrong** (repo commit 4c8328f): the pooled mean of log2(CPM+1) mixes detection rate with level;
  the real split is about 1.8-fold in detection and 2.0-fold in level among detecting cells. Pericytes rank
  240th yet read higher than astrocytes where detected. Top of the ranking: yes; oligodendrocyte-specific:
  not established.
- OPC expression varies by region (high cortex and hippocampus, low cerebellum and medulla), mature
  oligodendrocyte is flat; Levene p < 0.0001, replicated across cohorts (rho +0.96).
- Ageing (53 young vs 44 aged mice, one value per mouse): mature oligodendrocyte +0.250 log2 (1.19x),
  Mann-Whitney p = 0.0002, Cliff's delta +0.435; OPC +0.131, p = 0.22. Absent from Allen's own age-DE table.
  The ageing cohort covers only 7 of 13 regions (no thalamus, so habenula is untestable for age).
  **Does not replicate in human (section 9i).**
- Motor neuron hypothesis negative: spinal MN 0.24x of oligodendrocyte (donor-paired p = 0.023), cranial
  MN 0.40x of medulla neurons; same cells rank Kif5a 1st, Kif1a 6th, Kif13a 34th of 84. The earlier "L5 ET
  last in MOp" was a pooled ranking and is wrong; donor-paired, L5 ET is indistinguishable from other
  cortical neurons and the whole cortical block is low.
- Non-glial high: medial habenula cholinergic `145 MH Tac2 Glut` is 3rd brain-wide (+2.44 log2 versus the
  same mouse's other neurons, Bonferroni p = 0.0029); lateral habenula is not. Also medulla glycinergic
  nuclei, NTS Dbh, VMH.
- **Three analyses from 13 Sep exist only on TRUBA** (not in the repo or README): a cargo and adaptor
  panel (in medial habenula Kif13a ranks 2nd of 30 and its cargo Tfrc 1st; in cranial MN and MOp L5 ET the
  adaptors Rab11a, Rab10, Rab22a, Ap1m1 rank top while Kif13a ranks 27/71 and 13/15, so the adaptors are
  there and the motor is not); a within-habenula correlation scan over 20 clusters (Tfrc rho 0.877, depth
  rho only 0.283; cholinergic identity Chat, Slc18a3, Slc5a7 at 0.56 to 0.61; Kif13b 0.195); and a kinesin
  correlation against a 2,000 random-gene null (only Kif15, Klc4, Kif5a and Kif21a track Kif13a beyond
  the 95th percentile; Kif13b 0.470, Kif1a 0.178 and Rab11a 0.366 sit in the random range). Logs:
  `logs/cargo_6353444.out`, `hbscan_6353481.out`, `kinesin_6353670.out`; tables in `results_motor/`.

Consequence for the project: the mechanism transfers from KIF1A, the phenotype does not. If KIF13A
variants cause disease, the expected picture is white matter and myelin, not the neuronal axonal
transport disease that KIF1A causes.

## 9b. Loss-of-function axis (`docs/fonksiyon_kaybi_ekseni.md`, 2026-09-16)

Of KIF13A's 204 confidently observed gnomAD LoF variants, 45 escape nonsense-mediated decay and carry 456
of the 765 LoF alleles (59.6 percent). The most common allele, p.Trp1800Ter (AC 180, the only homozygote),
removes the last 5 of 1805 residues. Genuine NMD-competent high-confidence LoF: 158, against 102 for KIF1A
and 283 for KIF13B. The first real candidate by frequency is p.Arg1395Ter at AC 29.

The fact that changes the argument: **29 to 31 NMD-competent truncating variants sit inside the motor
domain in healthy people** (R268*, R173*, R207*, A356fs, AC 2-5). Early truncation of the motor domain is
therefore tolerated heterozygously, so a motor-domain missense causing severe dominant disease would have
to act by a dominant-negative mechanism, not by dose loss. That mechanism is shown for KIF1A
(Anazawa 2022) and is unmeasured for KIF13A. It is now the only coherent rationale left for tier A.

Two corrections to earlier reasoning in this project: **KIF1A is not a haploinsufficiency gene**
(ClinGen HI = 0), so it cannot serve as the comparator for dose sensitivity; and modest missense
constraint is not evidence against a loss-of-function mechanism (ANKRD11 has missense Z -1.21 and HI 3).
Against 12 ClinGen HI-3 neurodevelopmental genes, KIF13A is less LoF-constrained than all of them
(LOEUF 0.449 versus 0.094-0.369).

ClinVar holds zero classified truncating KIF13A variants; its 6 LoF-consequence records are TCGA splicing
submissions. Of 20 ClinVar CNVs only VCV000442594 has KIF13A as its sole protein-coding gene; every
pathogenic CNV spans 23 to 237 genes. IMPC: 21 significant phenotype calls, all homozygous, and
heterozygotes were never phenotyped beyond viability, so the decisive mouse test does not exist.

Verdict: loss of function fits better than missense, but only relatively. The evidence supports "possible
dose sensitivity with incomplete penetrance", not established haploinsufficiency. Seven falsification
tests are listed in section 6.3 of the document.


## 9c. Multi-pair feasibility (`docs/olcekleme_fizibilite.md`, 2026-09-16)

445 human genes carry at least 30 curated pathogenic missense variants in ClinVar; 3,198 pairs with a
clinically empty sibling (at most 2) were scored. Only 575 of them (18 percent, across 183 rich genes)
clear all four weakness filters. Median transfer lift is 0.048, 35.5 percent of pairs have zero or negative
lift, and 54 percent have an unconstrained sibling.

**Paralog transfer is therefore not a general property; it works in roughly one pair in six.** The claim
worth publishing is the predictive one: cheap metrics identify in advance which pairs will work.

KIF1A/KIF13A ranks 152 of 3,198 (71 of 337 deduplicated) with zero weakness flags, so this case study is
in the upper quartile rather than a lucky pick. Its lift looks modest (0.128) only because the motor
domain background identity is already 62 to 70 percent, which suppresses lift while making transfer easier.
Best five pairs: NSD1/NSD3, COMP/THBS1, PRKCG/PRKCE, MTM1/MTMR4, PCDH19/PCDH10. Recommended next targets
are NSD1/NSD3 and PCDH19/PCDH10 (neurodevelopmental, tight hotspots, constrained siblings), about a week
each, dominated by manual ClinVar curation. RYR1/RYR3 scores best but 5,038 aa proteins make it impractical.

## 9d. Loss-of-function tolerance against a curated comparison panel (2026-09-16)

Same pipeline as `scripts/lof_axis.py` (gnomAD v4.1 MANE transcript, LOFTEE HC, NMD-escape removed by the
Ensembl exon structure and the 50 nt rule), run over the 12 ClinGen haploinsufficiency-3 neurodevelopmental
genes so KIF13A's numbers can be read against genes where the answer is known. Script in the session
scratchpad; rerun it from `lof_axis` module functions, it needs no new code.

| gene | aa | NMD-competent HC LoF | per 1000 aa | raw per 100k | minus top allele | **obs/exp LoF** | most common allele |
|---|---|---|---|---|---|---|---|
| TCF4 | 671 | 7 | 10.4 | 0.43 | 0.37 | **0.037** | c.41_42del AC 1 |
| STXBP1 | 594 | 5 | 8.4 | 1.98 | 0.93 | 0.038 | c.*5del AC 17 |
| ADNP | 1102 | **0** | 0.0 | n/a | n/a | 0.047 | none |
| SCN1A | 2009 | 28 | 13.9 | 85.13 | 37.23 | 0.067 | c.1086C>A AC 773 |
| SYNGAP1 | 1343 | 33 | 24.6 | 13.13 | 9.60 | 0.071 | c.2354_2357del AC 57 |
| ANKRD11 | 2663 | 59 | 22.2 | 4.58 | 4.27 | 0.087 | c.3309dup AC 5 |
| NSD1 | 2696 | 46 | 17.1 | 4.15 | 3.84 | 0.107 | c.753_754insTT AC 5 |
| CHD8 | 2581 | 69 | 26.7 | 8.92 | 7.43 | 0.111 | c.7112del AC 24 |
| DYRK1A | 754 | 17 | 22.5 | 16.73 | 6.26 | 0.113 | c.948_949insGT AC 169 |
| POGZ | 1410 | 29 | 20.6 | 2.04 | 1.86 | 0.126 | c.284-1G>T AC 3 |
| ARID1B | 2372 | 188 | 79.3 | 124.16 | 74.72 | 0.141 | c.611_612insGC AC 798 |
| SETD5 | 1442 | 59 | 40.9 | 6.63 | 5.58 | 0.287 | c.*2del AC 17 |
| **KIF13A** | **1805** | **158** | **87.5** | **19.08** | **17.29** | **0.372** | c.4183C>T AC 29 |
| KIF13B | 1826 | 283 | 155.0 | 52.85 | 38.91 | 0.572 | c.3837+1G>T AC 225 |
| KIF1A | 1791 | 102 | 57.0 | 72.68 | 21.63 | 0.202 | c.4257-2A>C AC 824 |

**Read only the obs/exp column as evidence.** Raw carrier frequency is contaminated by single common
alleles that LOFTEE keeps but selection plainly does not act on (ARID1B c.611_612insGC at AC 798 sits in a
GC-repeat region; SCN1A c.1086C>A at AC 773; KIF1A c.4257-2A>C at AC 824). ARID1B is the cautionary case:
it is an established haploinsufficiency gene yet has the second-highest raw rate in the table.

On obs/exp, KIF13A (0.372) is **less loss-of-function constrained than all 12 established genes**
(0.037 to 0.287). ADNP carries literally zero NMD-competent high-confidence LoF alleles in gnomAD.

**What this kills and what it does not.** It rules out a severe, fully penetrant dominant disorder caused
by losing one KIF13A copy. It does not rule out a mild or incompletely penetrant effect, because gnomAD
excludes severe paediatric disease but not an adult who had speech delay as a child, and speech and
language delay is exactly the phenotype in the deletion cases. Say both halves whenever this table is used.

## 9e. De novo axis (2026-09-16, `docs/denovo_ekseni.md`)

| | KIF13A | KIF1A | KIF13B | ANKRD11 | CHD8 |
|---|---|---|---|---|---|
| Kaplanis 2020 (31,058 DD trios), observed de novo LoF | **0** | 0 | 2 | 111 | 19 |
| expected LoF | 0.367 | 0.298 | 0.337 | 0.308 | 0.497 |
| observed de novo missense | 2 | **53** | 6 | 17 | 9 |
| denovoWEST p | 0.709 | 5.2e-21 | 0.286 | 1e-14 | 4.3e-16 |
| all catalogues (about 97,700 effective trios), LoF | **0** | 2 (annotation errors) | 2 | 144 | 74 |

Verdict: **evidence against** haploinsufficiency of the ANKRD11/CHD8 type. The only de novo LoF found
(p.N1649Kfs*3) is in an unaffected sibling and in the last exon. Only effects above 2.6-fold are excluded
at 95 percent, so mild dose sensitivity survives. Three motor-domain de novo missense exist (D65E, R173Q,
Q237H; one from an ALS trio); popEVE scores them -3.62, -4.13, -3.42, none severe (9k).

## 9f. Regional missense constraint (2026-09-16, `docs/bolgesel_kisitlilik.md`)

Samocha 2014 trinucleotide rates, gnomAD v4.1 exomes, expectation calibrated on the panel's synonymous.

| KIF13A region | o/e missense | 90 percent interval |
|---|---|---|
| motor 5-352 | 0.814 | 0.745-0.892 |
| CC1 381-468 | 0.738 | 0.611-0.897 |
| tail 520-1805 | 0.927 | 0.887-0.968 |
| motor: ATP pocket (25 residues) | **0.523** | 0.362-0.776 |
| motor: microtubule interface (50) | 0.746 | 0.581-0.968 |
| motor: buried core (124) | 0.678 | 0.576-0.801 |
| motor: surface (153) | **1.008** | 0.892-1.140 |

Verdict: weakly constrained. The ATP pocket is real (23 tier A candidates sit there) but KIF1A's ATP pocket
is 0.140 and its motor 0.450, so KIF13A is two to four times shallower in every class and ranks 7th of 11
autosomal kinesins. gnomAD's regional constraint algorithm found no breakpoint; MTR and CCR point to
353-519, not the motor; no KIF13A variant reaches MPC 2 (maximum 1.18), so any MPC filter removes tier A.

## 9g. Out-of-sample pair NSD1/NSD3 (2026-09-17, `docs/nsd1_nsd3_dogrulama.md`)

The screen's top pair, run end to end with the prediction row hidden. AlphaMissense AUROC **0.973**
(0.942-0.995) inside the shared PHD1-CysHis module, gene threshold 0.9732; NSD3 candidates in
`data/processed/NSD3_candidates_20260917.tsv` (tiers A 23, B 45, C 33, D 30). The prediction held, with two
caveats: the screen's lift metric did not predict the measured paralog contribution (identity baselines
0.740 vs 0.744), and NSD3's shared module is not more constrained than the rest of NSD3. **The required
next experiment for a paper is a pair the screen predicts to fail**; both tested pairs are in the top 5
percent.

## 9h. White-matter imaging axis (2026-09-17, `docs/ak_madde_goruntuleme.md`)

Verdict **weak support**. UK Biobank BIG40 (up to 33,224 people, 3,935 IDPs, GRCh37): nothing reaches 5e-8
in white-matter classes within +/-100 kb; best p 3.5e-07 (rs144963738, left anterior corona radiata L3),
only 2.2-fold better than the LCT negative-control locus. Detectable lower bound at 80 percent power is
0.126 percent of variance; the observed haplotype effect is 0.086 percent, so about 46,000 people would
be needed. Three observations keep it alive: all 8 variants on the 3.6 percent haplotype show the same
white-matter profile; **rs76091407** (GRCh38 chr6:18,020,202, 32.6 kb from the KIF13A 5' end) is genome-wide
significant for white-matter microstructure (UDIP-FA, GCST90860937, PMID 42236496, N 11,736, p 4e-09);
the haplotype is r2 0.969 with MetaBrain's third KIF13A cortex eQTL signal (rs115824472). rs75428660 is
genome-wide significant for sulcal depth and cortical thickness, not white matter. The singulum OD hint
from the AZ PheWAS (n 44) was refuted (best p 1.4e-04 at N 31,355).

## 9i. Human lifespan single-cell, axis B (2026-10-03, `docs/insan_yasam_boyu.md`)

TRUBA `/arf/scratch/avural/kif13a_insan_yas`. Raw counts from CELLxGENE, no DUA needed: the Nature 2026
lifespan atlas (PMID 42778698) is on CELLxGENE as PsychAD "Aging_Cohort" (291 donors, 1.33 M nuclei),
plus SEA-AD MTG, a white-matter set (Seeker) and two developmental sets. Pseudobulk per donor, donor as
the test unit.

- **H_B(1), ageing: did not replicate.** Mature oligodendrocyte, 232 adults: -0.018 log2 per decade
  (p 0.017); young (20-39, n 55) vs old (60+, n 82) -0.059 log2 (0.96x), Cliff -0.21. Meta over 335 donors:
  old minus young -0.057 (95 percent CI -0.096 to -0.017), I2 0. Mouse: +0.250. Power to see +0.25 was
  above 0.99 (minimum detectable 0.085), so this is a real negative. The small decline is not specific
  (empirical p 0.078 genome-wide, 0.26 among expression-matched genes, 0.10 without the 90+ donors): say
  "does not increase", not "decreases". OPC unchanged (meta p 0.41).
- **H_B(2), lineage ranking: replicated with a twist.** OPC first, mature oligodendrocyte second in all
  three human sets (mouse is the reverse). Against astrocytes: detection 1.25x (depth-matched), level 2.0x
  (mouse 1.8x and 2.0x).
- **H_B(3), development:** OPC KIF13A rises from infancy in three sources (rho 0.43, 0.42, 0.82), half of
  the change done by about age 3 to 5. KIF13B and astrocytic KIF13A rise the same way, so it is general
  glial maturation, not a KIF13A or myelination program. Overlaps the speech-delay window in time only.
- Positive controls worked (PDGFRA in OPC p 3e-8; 6 of 7 microglia ageing markers p < 0.05). Missing: PMI
  correction and the authors' own dreamlet KIF13A rows, both OPEN on Synapse (`syn52396927`) and now
  reachable since the user has a Synapse account.

## 9j. Cell-type eQTL colocalisation, axis A (2026-10-03/04, PRELIMINARY, no document)

The agent was stopped by the user before writing `docs/eqtl_hucre_tipi.md`. Results exist and were read
for this file; treat them as preliminary until a document is written. Local tables
`data/interim/eqtl_hucre_tipi/`, scripts `scripts/truba/eqtl/`, TRUBA `/arf/scratch/avural/kif13a_eqtl`
(log `logs/eqtl_pipe_6440189.out`).

- Sources: Zeng ... Roussos 2026 Nat Genet (PMID 42778721, 1,384 donors, Synapse `syn61929919`, effective
  N about 3,450 per class), SingleBrain 2025 meta-analysis, Bryois 2022 (Zenodo 7276971, effective N 200
  to 1,000). GWAS: BIG40 L3 left anterior corona radiata (UKB 25319, lead rs144963738) and UDIP-FA
  (GCST90860937, lead rs74819634). LD from 1000 Genomes EUR and all.
- Validation passed: Python coloc matches R coloc 5.2.3 (max difference 2.9e-05); positive controls
  reproduce published colocalisations exactly (Bryois S5 CR1/NDUFS2/ADAMTS4 with AD, PP.H4 0.993/0.926/
  0.768; Zeng CR1 and ADAMTS4 with AD about 0.96 and 0.94); BIG40 build test confirms GRCh37.
- **Result against H_A.** KIF13A has strong eQTLs in oligodendrocytes (Zeng lead rs113815188, p 2.7e-14)
  and OPCs (rs760849, p 3.6e-37), but they are different variants from the white-matter signal (r2 with
  the GWAS lead at most 0.05; eQTL p at the GWAS lead 0.78 in oligodendrocytes, 0.0045 marginal and about
  0.9 conditional in OPC). PP.H3 0.95 to 0.9998 and PP.H4 below 0.002 in oligodendrocyte and OPC, in Zeng
  and SingleBrain, for both GWAS. Astrocytes and neurons the same (PP.H4 at most 0.02).
- **The only hint is microglia, and it is fragile.** Zeng immune, third conditional signal rs73373199 (r2
  0.30 with the GWAS lead): PP.H4 0.55 (BIG40) and 0.62 (UDIP-FA) with the EUR panel, but 0.002 and 0.034
  with the all-1000G panel. SingleBrain immune, conditional signal is rs75428660 itself (r2 0.97 with the
  lead) with PP.H4 0.17 and an LD-mismatch flag. Not a finding until it survives in-sample LD.
- Bryois 2022 is underpowered for a 3 percent variant (power simulation, `power_sim.tsv`: even at N 1,100
  and the largest simulated effect, 0.388 SD per allele, power to reach PP.H4 >= 0.8 is only 0.32 to
  0.34), so its low PP.H4 values are uninformative.

## 9k. Biallelic knockout axis and popEVE, axis C (2026-10-03, `docs/nakavt_ekseni.md`)

Verdict **uninformative**. TRUBA `/arf/scratch/avural/kif13a_nakavt`. The design cannot even separate
KIF1A, a proven biallelic LoF disease gene (likelihood ratio 1.00 to 1.03).

| Cohort | observed KIF13A knockouts | expected | note |
|---|---|---|---|
| Genes & Health, 44,028 autozygous adults (PMID 41896352) | 0 (22 heterozygotes) | 0.17 (true LoF 0.10) | P(0) 0.84 |
| RGC-ME, 983,578 exomes (PMID 38768635) | **2 homozygotes** | 5.1 | table gives only cumulative AAF 7.73e-5, which matches c.4581+1G>C (NMD escape, in-frame if exon 38 skips); probably not true knockouts (inference, unverified); phenotypes unpublished |
| gnomAD v4.1 | 0 true-LoF homozygotes | 0.46 | |
| Chundru 2024, 29,745 DD trios (PMID 39313616) | 0 LoF/LoF | 0.0063 neutral, 0.052 under H_C | 474-fold enrichment needed for 95 percent power |

Knockout counts are over-dispersed: in RGC 59 percent of tolerant genes at the same frequency also show
zero, so a Poisson "p 0.049" is fake. At gene-set level the method works (DDG2P recessive DD genes 0.26 to
0.35-fold depleted). Human action needed: open the RGC-ME browser (licence click) and read the homozygote
counts at 6-17771113-C-G and 6-17764129-C-T.

**popEVE (Nat Genet 2025, PMID 41286104).** KIF13A is neither among the 442 genes nor the 123 new
candidates. Severe threshold popEVE <= -5.056: KIF13A tier A 1 of 52 severe (S278W), **the 68 changes that
are pathogenic in KIF1A at the identical residue: 0 severe in KIF13A**, against 86 of 111 (77.5 percent)
for KIF1A's own core positives. This is the cross-gene version of 9f: same residue, same change, severe in
KIF1A and moderate in KIF13A, because popEVE folds in each gene's human constraint. It argues against a
severe developmental disorder, not against dominant negative or milder or non-DD phenotypes. popEVE's
filtered score covers only residues 1-939 of KIF13A. popEVE's own AUROC on the KIF1A set (0.957) cannot
enter the calibration table, because its population component uses the same gnomAD data as the negatives.

## 9l. AlphaGenome Atlas and the AVI score (evaluated 2026-10-04, not run)

Announced 8 Sep 2026; preprint Cheng ... Avsec, medRxiv, doi 10.64898/2026.09.16.26363192, PMID 42779988,
PMC13596691. Precomputed predictions for all about 9 billion SNVs; AVI is a network over 18 features (10
AlphaGenome modalities as maximum across tissues, AlphaMissense, 3 termination features, PhastCons,
Cactus, 2 indel flags), trained on gnomAD frequency proxy labels (AF > 0.001 benign, < 0.001 impactful),
ClinVar used only for evaluation, PHRED-normalised genome-wide. Access by API key under non-clinical terms,
CC BY-NC 4.0; the user must accept the terms and place the key.

What it means here: it may enter the calibration table (not label-trained), but for motor-domain missense
its protein signal is AlphaMissense, so expect little gain (measure, do not assume). It carries no gene
constraint, so it will score KIF13A ATP-pocket variants like KIF1A's; disagreement with popEVE is by
design. Real uses: the untouched splice axis (26 KIF1A splice positives, K1122E at an exon boundary),
choosing the causal variant on the white-matter haplotype **only if** a colocalisation ever succeeds, and
cross-checking LoF/NMD annotation. Its tissue list has no oligodendrocyte, OPC or white-matter track.

## 10. Verified literature

Boyle 2021 HGG Adv PMID 33880452 (severity and motility classes) - Chiba 2019 PNAS PMID 31455732
(hyperactivating SPG variants) - Anazawa 2022 PNAS PMID 35917346 (dominant-negative KAND) -
Shatarupa 2026 Nat Commun PMID 41912506 (R350 salt bridge) - Kaur 2020 Hum Mutat PMID 32652677 -
Nicita 2021 J Med Genet PMID 32737135 - Lam 2021 J Hum Genet PMID 33526817 (the KIF13A mosaic case) -
Evo 2, Nature 2026 PMID 41781614.

Added 2026-10-03/04, each verified live on PubMed or Europe PMC:

| Use | Reference |
|---|---|
| cell-type brain eQTL (axis A) | Zeng ... Roussos, Nat Genet 2026, PMID 42778721, doi 10.1038/s41588-026-02733-5 |
| lifespan PFC atlas (axis B) | Nature 2026, PMID 42778698, doi 10.1038/s41586-026-10271-7 (preprint PMID 41282859) |
| Genes & Health knockouts (axis C) | Kim ... van Heel, Nat Genet 2026, PMID 41896352, doi 10.1038/s41588-026-02553-7 |
| RGC Million Exome (axis C) | Sun 2024 Nature, PMID 38768635, doi 10.1038/s41586-024-07556-0 |
| recessive DD trios (axis C) | Chundru ... Martin, Nat Genet 2024, PMID 39313616, doi 10.1038/s41588-024-01910-8 |
| popEVE | Orenbuch ... Marks, Nat Genet 2025, PMID 41286104, doi 10.1038/s41588-025-02400-1 |
| 100,000 Genomes gene burden (not used yet) | Cipriani 2025 Nature, PMID 40011789, doi 10.1038/s41586-025-08623-w |
| AlphaGenome Atlas / AVI (not used yet) | Cheng ... Avsec, medRxiv 2026, PMID 42779988, doi 10.64898/2026.09.16.26363192 |
| UDIP-FA white-matter GWAS | Zhao X, Nat Commun 2026, PMID 42236496, GCST90860937 |
| mechanism prediction (pLOF/pGOF/pDN, the DECIPHER numbers) | Badonyi 2024 PLoS One, PMID 39172982, doi 10.1371/journal.pone.0307312 |
| mechanism prevalence across disease | Badonyi 2025 Nat Commun, PMID 40998763, doi 10.1038/s41467-025-63234-3 |
| ClinGen PP3/BP4 calibration incl. AlphaMissense | Bergquist 2025 Genet Med, PMID 40084623, doi 10.1016/j.gim.2025.101402 |
| Human Domainome (500 domains, stability DMS) | Beltran 2025 Nature, PMID 39779847, doi 10.1038/s41586-024-08370-4 |

## 11. Files produced

- `config/targets.yaml` verified identifiers and the isoform trap (NSD1/NSD3 added 2026-09-17)
- `docs/MCP_KURULUM.md` data layer, `docs/KIF1A_ClinVar_kurasyonu.md` curation rules,
  `docs/kalibrasyon_olcumleri.md` measurements, `docs/yapi_ekseni.md` structure axis,
  `docs/vaka_arama.md` case-finding routes
- `docs/fonksiyon_kaybi_ekseni.md` (LoF), `docs/missense_sonuc.md` (saturation table),
  `docs/olcekleme_fizibilite.md` (multi-pair screen), `docs/dogrulama_delesyon.md` and
  `docs/dogrulama_hudsonalpha.md` (the two adversarial audits), all 2026-09-16
- `docs/denovo_ekseni.md` (9e) and `docs/bolgesel_kisitlilik.md` (9f), 2026-09-16
- `docs/nsd1_nsd3_dogrulama.md` (9g) and `docs/ak_madde_goruntuleme.md` (9h), 2026-09-17
- `docs/insan_yasam_boyu.md` (9i) and `docs/nakavt_ekseni.md` (9k), 2026-10-03
- **missing:** `docs/eqtl_hucre_tipi.md` (9j); tables exist in `data/interim/eqtl_hucre_tipi/`
- `data/processed/KIF13A_candidates_20260915.tsv`, `KIF13A_saturation_lookup_20260916.tsv`,
  `NSD3_candidates_20260917.tsv`; figures in `figures/`
- `data/raw/PROVENANCE.jsonl` every MCP-era pull; the TRUBA axes keep their own manifests
  (`data/interim/nakavt/SOURCES_MANIFEST.jsonl`, `scripts/truba/insan_yas/datasets.tsv`)
- `HANDOFF.md` at the repo root: the model-agnostic handover in Turkish. **It was last updated 2026-09-17
  and does not yet contain 9i to 9l.** Keep the two in sync; if they disagree, the repo files win.
- The git repository still has **no commits**.

## 12. Scoreboard (2026-10-04)

| Axis | Question | Verdict | Where |
|---|---|---|---|
| Calibration | do VEPs work on the teacher | yes, AlphaMissense 0.925 within motor | 5 |
| Out-of-sample pair | does the method port | yes, NSD1/NSD3 0.973 | 9g |
| Structure | where do tier A candidates sit | 79 percent at functional sites | 8 |
| Gene constraint | is KIF13A missense constrained | weakly (Z 2.79 vs 6.87) | 2 |
| Regional constraint | is the motor constrained locally | ATP pocket yes (0.523) but 4x shallower than KIF1A | 9f |
| popEVE | are KIF1A-identical changes severe in KIF13A | no, 0 of 68 | 9k |
| LoF tolerance | haploinsufficiency like HI-3 genes | no, less constrained than all 12 | 9d |
| De novo | dose sensitivity in DD trios | evidence against (0 LoF, >2.6x excluded) | 9e |
| Recessive knockouts | biallelic LoF disease | uninformative (no power, even for KIF1A) | 9k |
| Human observations | any case that carries alone | no; both leads at lowest rung | 6, 6b |
| Cell type | oligodendrocyte lineage | yes, mouse and human (OPC first in human) | 9, 9i |
| Ageing | mouse increase replicates in human | no, real negative | 9i |
| Imaging GWAS | white-matter trace at the locus | weak support, sub-threshold | 9h |
| eQTL colocalisation | WM signal acts via oligodendrocyte KIF13A | against (PP.H3 high); fragile microglia hint | 9j |
| Dominant negative | does a full-length mutant poison the motor | **never measured** | 6, 9b |
