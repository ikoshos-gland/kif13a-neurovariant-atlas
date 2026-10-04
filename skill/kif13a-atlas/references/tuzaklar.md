# Traps

Each of these was hit in this project, cost time, and would have produced a wrong or fake-clean result if
it had gone unnoticed. Read this before writing a script that touches the data, and again when a result
looks suspiciously good.

## Numbering and identity

1. **Isoform length.** The project plan said KIF13A is 1770 aa. The MANE protein is 1805 aa; 1770 is
   isoform Q9H1H9-2. A saturation matrix built on the wrong length is wrong everywhere.
2. **The teacher's clinical numbering is not its UniProt canonical.** KIF1A clinical reports use
   NM_001244008 (1791 aa, UniProt Q12756-3). The canonical Q12756 is 1690 aa. Positions agree only to 393;
   after that add 9, and after canonical 848 add 101. AlphaMissense per-protein tables and the canonical
   AlphaFold model use the canonical numbering, so a tail variant silently lands on the wrong residue.
   VEP's AlphaMissense plugin is transcript-aware and avoids this.
3. **UniProt feature annotations are on the canonical isoform.** Fetch features with the canonical
   accession and apply the offsets, otherwise domain boundaries drift.

## Variant data

4. **gnomAD's gene query is not transcript-aware.** It reports HGVSp on the most severe transcript.
   KIF13A missense rows spread across 5 transcripts that way. Use the transcript query on the MANE
   transcript: 10,936 variants, 2,671 missense, all numbered consistently.
5. **gnomAD lists only exons plus 75 flanking bases.** Deep intronic variants are simply absent, so
   "all variants in the gene" is never what you have.
6. **The AC0 filter.** 596 of the 2,671 KIF13A missense rows have allele count 0 and filter AC0: every
   carrier genotype failed quality control. They are not confidently observed variants. Eleven rows carry
   AC0 yet have one genome carrier; treat them conservatively and say which rule you used.
7. **ClinVar renamed the classification field.** It is `germline_classification`; the old
   `clinical_significance` returns nothing, which is why several third-party ClinVar MCP servers appear
   to work while returning empty classifications.
8. **ClinVar 0-star records can still be the best evidence.** A255V's only record is an OMIM entry with no
   assertion criteria, yet two functional papers characterise it. A pure star filter throws it away.
9. **Population absence is the norm for real disease variants.** 89 percent of KIF1A's curated pathogenic
   missense are absent from gnomAD. A candidate list built only from gnomAD-observed variants
   systematically excludes what the project is looking for.

## Benchmarking

10. **Location confound.** ClinVar's benign missense for KIF1A contains zero motor-domain variants while
    the positives are almost all there. "Is it in the motor domain" scores AUROC 0.995 on that comparison,
    beating every model. Negatives must come from the same region as the positives.
11. **Mixed score sources.** Positives scored by VEP and negatives by gnomAD's annotation gave CADD 0.914;
    one uniform VEP run gave 0.781. Always score both classes in the same call.
12. **Circularity.** 75 of 135 KIF1A pathogenic missense comments cite in-silico evidence. Variants whose
    classification leaned on PP3 are held out; on them AlphaMissense scores 0.961 versus 0.925 on the core
    set. REVEL, MetaRNN and ClinPred are trained on clinical labels and cannot be evaluation targets or
    consensus members here; they stay as reference rows.
13. **Pseudo-replication.** Several substitutions of one residue are not independent. Bootstrap and
    cross-validate by residue, not by variant.
14. **Population negatives have label noise.** The AC >= 3 set still contains variants that may be mildly
    damaging, and AlphaMissense was trained on population frequency, which flatters it against
    gnomAD-derived negatives. Report CADD alongside as the least circular reference.

## Tools and APIs

15. **VEP rejects unsorted VCF-style input.** Sort region lines by chromosome and position before posting,
    or the batch fails with "input file appears to be unsorted". Batches of 200 are fine.
16. **dbNSFP fields come as per-transcript lists.** ESM1b and MetaRNN values are comma-separated across all
    transcripts; index them against `Ensembl_transcriptid` to get the MANE value.
17. **AlphaFold DB serves AlphaMissense per protein.** `AF-<acc>-F1-aa-substitutions.csv` is a few hundred
    kilobytes, so the 1.2 GB whole-proteome download is unnecessary. Isoform accessions have a model but
    no AlphaMissense table.
18. **GTEx needs the versioned GENCODE id for the release** (v10 wants ENSG00000137177.20); the older id
    silently returns an empty list.
19. **CELLxGENE Census needs tiledbsoma, which has no Windows build.** The Where's My Gene API serves the
    same corpus pre-aggregated and works fine. Its `pc` field is a fraction despite the docs saying 0-100.
20. **Evo 2 on the hosted NIM is not bitwise causal.** Changing one base shifts upstream log-probabilities
    by up to 0.36 (about 0.4 nats summed over 4096 bases). Score from the variant base onward. Treat any
    absolute score below about 3 as noise.
21. **Fetch one genomic region and slice locally.** Per-variant Ensembl window fetches made Evo 2 scoring
    18 minutes per variant; one region fetch plus four concurrent NIM calls brought it to 3.3 seconds.
22. **File globs match the provenance sidecars.** `glob("*_constraint_*.json")` returns the `.prov.json`
    files too and `sorted(...)[-1]` picks the wrong one. Exclude them explicitly.

## Structure

23. **A superposed model is not a structure.** Rigid-body placement of the KIF13A AlphaFold model on the
    KIF1A complex produces sub-2 A clashes and inflates the interface residue count (50 against a native
    38). Contact yes or no is usable; the distance value is not.
24. **Both motor heads matter.** R350 contacts tubulin in one head and not the other; using a single chain
    would have hidden the contact that the literature describes.
25. **Read the chains from the file.** Do not assume which chain is the kinesin; parse the structure's own
    entity records.

## Appraising human evidence (added 2026-09-16, every one hit during the verification round)

26. **A ratio is not evidence until its denominator is measured.** This project argued that HudsonAlpha
    having 1,293 ClinVar submissions and exactly one in KIF13A meant a deliberate candidate-gene call.
    Measuring their baseline destroyed it: 368 of their 607 single-gene genes (60.6 percent) have exactly
    one record, so one record per gene is their *most common* pattern. Before reading a count as a signal,
    compute what that submitter, cohort or database does by default.
27. **ClinVar has two different "affected" fields.** `AffectedStatus` is the observation's own status;
    `ClinicalFeaturesAffectedStatus` only says the listed HPO terms are present. Reading the second as the
    first turned an `unknown` into a `present` in `docs/fonksiyon_kaybi_ekseni.md`.
28. **A ClinVar record's cited PMID may never have studied the locus.** VCV000442594 cites Kaminsky 2011,
    which analysed only 14 recurrent CNV regions and never touched 6p22.3. Open the citation before
    treating it as support, and read the record's last-evaluated date: 2014 predates gnomAD SV entirely.
29. **Submitter comments are often templates.** Every nstd37 submission carries the same comment text. A
    sentence that appears in every record of a study carries no information about any one of them.
30. **gnomAD v4 has no subcohorts.** The `StructuralVariantDatasetId` enum offers `_controls` and
    `_non_neuro` only for v2.1. You cannot state that a v4 carrier is neurologically normal, only that
    they were outside a severe paediatric disease ascertainment. And v2.1 is small: for a 1-in-175,000
    event the expected count there is 0.03 to 0.06, so a zero is a power failure, not a result. Compute the
    expected count before reporting any absence as meaningful.
31. **Intragenic and gene-spanning deletions are different evidence classes.** A deletion smaller than the
    gene and lying inside its boundaries is single-gene by geometry and needs no database assertion
    (DECIPHER 451718, 146 kb inside a 230 kb gene). A deletion larger than the gene sticks out at both ends
    and needs a regulatory-content check: the 253 kb case removes 37 enhancers, 9 CTCF sites and a
    promoter, and its nearest neighbour NUP153 is more constrained than KIF13A itself.
32. **"One variant deposited" is not "a negative exome".** DECIPHER patient 451718 shows Genotype 1, but
    most DECIPHER submissions are microarray only, so point mutations in other genes were never looked
    for. A single-gene claim from array data is a claim about copy number and nothing else.
33. **Case databases have no control arm.** DECIPHER, dbVar nstd101 and its replacement nstd102 are all
    Case-Set studies, and nstd102's 21 KIF13A calls have empty `sample_id`, `validation` and `zygosity`.
    Everyone in them has a phenotype, so a variant found in a phenotyped patient is not informative
    without a population rate to compare against.
34. **Check a phenotype's base rate before counting it as specific.** Chiari malformation appears in 3.6
    percent of children imaged for any indication, and in a developmental delay cohort 6 of 6 were judged
    unrelated to aetiology. Speech and language delay has a population prevalence of 3.2 to 25.6 percent
    and microarray referral yields only 15 to 20 percent diagnoses. An incidental finding inside a
    phenotype list inflates apparent specificity.
35. **Inheritance is the highest-value missing field.** Of KIF13A's four loss-of-function events with known
    inheritance in developmental cohorts, three are maternal and none is constitutional de novo. Maternal
    inheritance does not prove harmlessness (reduced penetrance, a mildly affected mother), but it points
    the opposite way to de novo. When a record lacks the field, state what the claim costs.
36. **Haploinsufficiency predictors disagree, so never quote one alone.** For KIF13A pHaplo 0.99 says dose
    sensitive while sHet 0.049 and DECIPHER %HI 42.3 (threshold 10) say it is not. The disagreement is the
    finding.
37. **Read the paralog's own allelic requirement before assuming a mechanism.** DDG2P gives KIF1A two
    entries: monoallelic with "altered gene product structure", biallelic with "absent gene product / loss
    of function". Losing one copy is not how KIF1A causes disease, so assuming it for KIF13A assumes
    something the teacher gene does not do.
38. **Raw loss-of-function carrier frequency is contaminated by single common alleles.** ARID1B is an
    established haploinsufficiency gene yet has the second-highest raw rate in the comparison panel,
    because of one 798-allele insertion in a GC repeat. Use observed/expected, and print the top allele
    beside any raw rate.
39. **Citation silence is usable evidence.** Lam 2021 has been cited 3 times in 5 years and not one citing
    paper reports an independent KIF13A case. For a candidate disease gene, five quiet years is data.

## Environment (Windows)

40. Microsoft Store `python` virtualises `%APPDATA%`, so it cannot see uv tools or uv's interpreters. It
    can still run the project's own scripts; see the next trap.
41. **The uv-managed CPython 3.12.14 behind `mcp\atlas_server\.venv` was blocked by Windows on
    2026-09-16** ("Permission denied" on the interpreter itself, "Unable to create process" through the
    venv shim), so the venv cannot launch at all. Verified replacement: the system `python` 3.13.14, which
    already carries httpx, numpy, pyyaml, biopython, pandas, scipy and mcp, imports `atlas_mcp.core` and
    reads `.env`. Every reproduction command in `HANDOFF.md` and the user-scope MCP registrations still
    name the dead venv and need re-pointing.
42. Windows Smart App Control blocks some generated launcher exes (seen with `zotero-mcp.exe`); run the
    module through a working python instead. Do not touch the security setting.
43. An MCP server over stdio deadlocks if numpy is first imported inside a tool call; import heavy
    libraries at module load.
44. PyPI is slow on this machine and two concurrent uv installs fight over the cache lock. Install serially.
45. Background jobs stop when the machine sleeps. A job left overnight may be exactly where you left it.

## De novo and ClinVar queries (added 2026-09-16/17, from `HANDOFF.md` 14 to 26)

46. **ClinVar Entrez has no `pathogenic[clinsig]` field.** Entrez silently rewrites it to
    `pathogenic[All Fields]` and returns an inflated count (692 for KIF1A, 371 for NSD1). The right form is
    `("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])`, which gives 135 and
    158 and matches manual curation. The trap is in the query syntax, not the gene.
47. **Kaplanis 2020's supplementary table is GRCh37 and never says so.** Read as GRCh38, all four KIF13A rows
    change consequence class. Annotate every row in both builds and test.
48. **Gene4Denovo2's `AAChange_refGene` uses NM_001105567/8 for KIF13A, 13 residues below MANE.** Use
    `AAChange_mrefGene` (NM_022113). This is a third isoform shift, distinct from the 1770 aa trap.
49. **denovo-db writes one row per transcript.** Counting rows inflates counts up to fourfold (KIF13A's
    "101 rows" are 25 variants).
50. **Gene4Denovo2's LoF label is not trustworthy alone**: a 3 bp in-frame insertion in KIF1A is labelled
    `stopgain`; Kaplanis calls it `inframe_insertion`.
51. **A Poisson p value as `1 - CDF` underflows on positive controls** and prints a fake `p = 0`. Sum the
    upper tail in log space.
52. **Zero observations is not zero effect.** Expected de novo LoF in KIF13A is 0.37 (Kaplanis) and 1.16
    (all catalogues), so only effects above 8.2x and 2.6x are excluded. State the excluded effect size.
53. **NSD3's old symbol WHSC1L1 makes Ensembl `lookup/symbol` return 400**, and NSD3's well-studied short
    isoform Q9BZ95-3 (645 aa) lacks the SET domain; its numbering cannot be compared with MANE.
54. Ensembl REST returned intermittent 500 and 503 on 2026-09-17; five retries were not enough, nine were.

## Imaging GWAS (added 2026-09-17, from `docs/ak_madde_goruntuleme.md`)

55. **BIG40 is GRCh37, GWAS Catalog REST returns GRCh38, and at this locus the two differ by only 231 bp.**
    rs75428660 falls inside KIF13A in both builds, so the wrong build looks right. Resolve every variant in
    both builds.
56. **BIG40's PheWeb copy keeps only p <= 0.01 and MAF >= 0.01.** Absence there is not absence; MAF 0.001 to
    0.01 lives only in the full summary statistics. The full files are single-member gzip (no BGZF), so no
    HTTP range access; reaching chr6:18 Mb streams about 30 percent of each file.
57. **BIG40 PheWeb returns a persistent HTTP 500 for one white-matter IDP** (TBSS FA posterior corona
    radiata L, pheno 1479). Log failures and fill them from the full statistics, or the white-matter scan
    has a hole.
58. **Ensembl's `MAF` is the 1000 Genomes global frequency** (0.0092 for rs75428660 against 0.036 in gnomAD
    NFE and 0.043 in UKB). Use population-matched frequencies before calling a mismatch.
59. **Ensembl LD `pairwise` returns an empty list both for "low r2" and "variant not in panel".** Use the
    window query and check whether the target is in the returned set.
60. **The eQTL Catalogue REST API is gone (HTTP 410).** MetaBrain's DataTables backend
    (`eqtlgen.molgeniscloud.org/api/v2/metabrain_cis_eqtl_v2?q=GeneSymbol==KIF13A`) works. In GTEx,
    `maf01: false` means "not tested", not "not an eQTL".
61. **Window width changes the answer.** +/-250 kb reaches TPMT, whose common haplotype is richer in
    white-matter associations and unlinked to KIF13A (r2 < 0.05). Use +/-100 kb for gene attribution.

## The 2026-10-03 TRUBA axes (eQTL, human lifespan, knockouts)

62. **A Synapse token can read but not download.** Synapse permissions said `canDownload: true` while every
    file request returned 403; the active token had only the `view` scope. Check scopes with
    `GET auth/v1/personalAccessToken` (print names and scopes only). Anonymous download is closed even for
    open folders (`canPublicRead` true, `canDownload` false); certification was not the cause.
63. **Conditional colocalisation with external LD is fragile.** The Zeng immune PP.H4 0.55 to 0.62 for the
    white-matter signal appears with the 1000G EUR panel and disappears with the all-1000G panel, and the
    pipeline's own LD-mismatch flag fires on several steps. Treat any PP.H4 that depends on the LD panel as
    a hypothesis until in-sample LD confirms it.
64. **Validate a home-made coloc against the reference implementation and a published result.** The Python
    port matched R coloc 5.2.3 to 2.9e-05 and reproduced Bryois S5 and Zeng AD colocalisations exactly;
    without that, a negative coloc is uninterpretable.
65. **CELLxGENE puts raw counts in different places per dataset.** In PsychAD `X` is raw and `raw` is
    absent; elsewhere `X` is normalised and raw counts are in `raw/X`. Check integer-ness before summing,
    and recompute log2(CPM+1) from raw counts.
66. **CELLxGENE's "80 year-old and over stage" means 90 and over**, and more than half of SEA-AD donors carry
    it. Treat it as censored; run the analysis with and without.
67. **PsychAD cohorts are nested.** The lifespan cohort has `is_primary_data = False`: it copies nuclei from
    HBCC, MSSM and RADC. Using HBCC as a "replication set" counts the same donors twice.
68. **The mouse per-cell metric is noisy in nucleus data.** Mean log2(CPM+1) over single-count nuclei is
    nearly log2(1e6/UMI), and "level where detected" inherits depth. Use pseudobulk for age tests and depth
    stratification for cell-type comparisons.
69. **`sbatch --export=KEYS=a,b,c` splits on commas** and passes only `a`. Use another separator.
70. **Knockout counts are not Poisson.** Founder alleles and consanguineous subcohorts over-disperse them;
    in RGC 59 percent of tolerant genes at KIF13A's frequency also show zero. Compare against the empirical
    distribution of genes in the same frequency band.
71. **RGC-ME supplementary table 6 gives cumulative homozygote AAF, not variants**; identities need the
    licensed browser. Its table 2 uses non-MANE ENST00000649096 for KIF1A, which drops the AC 824
    c.4257-2A>C and shifts q fiftyfold. Check the transcript before moving a LoF frequency between cohorts.
72. **popEVE's filtered score leaves high-gap positions empty**: KIF13A 940-1805 is only in the unfiltered
    column. Percentages use non-empty scores as the denominator (that is how its Supp Table 2 matches).
73. **Chundru's supplementary data spells a header "Exepected"** and its zip carries a `__MACOSX/._*.xlsx`
    resource fork that crashes Excel readers.
74. **Genes & Health's "publicly available" bucket returns 403 to anonymous access** (VPC Service Controls).
75. **Only TRUBA's `debug` partition has internet.** Downloads run there as jobs (4 h limit); compute
    partitions are for the work after.
76. **A gene-calibrated threshold and a cross-gene severity score answer different questions.** AlphaMissense
    0.968 ranks variants within KIF13A; popEVE asks whether a variant is severe at organism level and folds
    in each gene's human constraint. The same KIF1A-pathogenic change is severe in KIF1A and moderate in
    KIF13A by design. Never read a high within-gene score as "causes disease in this gene".
77. **2 raised to a difference of pooled log means is not a fold change.** The Allen ranking averages
    log2(CPM+1) over all cells with non-detected cells at 0, which blends detection rate and level; read as
    a fold change it gave "16x above astrocytes" where the real split is 1.8x detection and 2.0x level.
    Report the two components separately.
