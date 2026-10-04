# The method: paralog-calibrated variant prioritisation

The problem this solves: a gene with no clinical variant data cannot tell you whether your prediction
tools work on it. Any prioritisation you produce is unfalsifiable. The trick is to borrow a sibling gene
that has the clinical data, measure the tools there, and carry across only what legitimately transfers.

What transfers is the protein mechanism, because paralogs share the machine. What does not transfer is
the phenotype, because paralogs are expressed in different cells. Keeping those two apart is the whole
discipline of the method.

## Step 1: fix identity before anything else

Pick the MANE Select transcript for both genes and record the protein length, the RefSeq and Ensembl ids,
and the UniProt accession of the isoform that matches MANE. Check whether the UniProt canonical isoform
equals the MANE translation; when it does not, write the offset table down and keep it where every script
can see it. Do this first because every later number silently depends on it.

## Step 2: build the variant universe

Query gnomAD through the transcript, not the gene, so that protein positions are consistent. Keep the
quality filter status per variant, since a variant with no confidently called carrier is not the same
thing as an observed one. Pull the full ClinVar record set for both genes.

Also build the theoretical universe: every possible substitution of the target protein, flagged for
whether a single nucleotide change can produce it. Disease variants are usually absent from population
data, so the observed set alone is the wrong hunting ground.

## Step 3: map residues between the pair

Align the two MANE proteins and produce a position-to-position map with an identity flag. Report identity
by region, because a single global number hides the structure of the relationship: in the KIF1A and
KIF13A pair the motor domain is 62 percent identical, the tail 34 percent, and 386 residues of KIF1A have
no partner at all.

Check the map against known landmarks before trusting it. Conserved functional motifs are ideal:
the P-loop, switch I and switch II in kinesins, and the residues UniProt annotates as ligand-binding.

## Step 4: curate labels on the teacher gene

Positives: curated pathogenic or likely pathogenic missense with review criteria, absent or nearly absent
from population data, and with the subset whose classification leaned on in-silico predictions held out
separately. That held-out subset is a probe: if a model scores much better on it than on the core set,
the classifications were partly downstream of predictions like it.

Negatives: this is where the method usually breaks. Clinical benign variants tend to sit in different
regions from pathogenic ones, so a classifier that only knows position will beat every real model. Draw
negatives from the same region as the positives. Population variants observed in several unrelated people
work well for a severe dominant disorder, with the allele-count threshold stated and varied as a
sensitivity check.

Keep separate sets for other consequence classes (splice, loss of function) rather than mixing them in.

## Step 5: measure honestly

Score both classes in one pass through one annotation source. Compute AUROC with a bootstrap that
resamples residues rather than variants, since several substitutions of one residue are not independent
observations. Report, next to every model:

- the confidence interval
- sensitivity at a fixed specificity, and the score threshold that achieves it, which is the number you
  actually transfer
- at least one baseline that uses no model at all, such as "is this residue in the functional domain" or
  "is the residue identical in the sibling"

A model that cannot beat those baselines has told you nothing. Models trained on the same clinical labels
you are testing against belong in the table as reference rows, never as the headline or as consensus
members.

## Step 6: transfer, then add independent axes

Carry the gene-specific threshold and the residue-level paralog evidence to the target gene. Tier the
result by which evidence supports it rather than by a fitted score, so a reader can see why a variant is
where it is.

Then add axes that did not contribute to the tiering, because convergence is the argument:

- **structure**: superpose the target's predicted model onto an experimental complex of the sibling and
  ask where each candidate sits, treating contacts as qualitative
- **loss of function**: constraint metrics, surviving loss-of-function variants and their position
  relative to nonsense-mediated decay, deletions, mouse knockouts
- **cell type**: where the target and the sibling are expressed, which tells you which phenotype to expect
  and which to stop expecting

Axes this project added later, each with its own document, positive and negative controls and a power
statement: de novo enrichment in developmental-disorder trios, regional missense constraint by structural
class, a population imaging GWAS for the phenotype the cell type predicts, single-nucleus eQTL
colocalisation for that GWAS signal, a human lifespan single-cell replication of the mouse finding, and
biallelic knockout counts against an expected number. The pattern to copy: write the verdict before the
details, compute the expected count or the detectable effect for every negative, and run the same test
on the teacher gene so that "uninformative" can be told apart from "negative".

## Step 6b: check the transfer with a cross-gene severity model

A gene-calibrated threshold ranks variants inside the target gene; it says nothing about whether the
target gene can produce disease at all. Run a model that scores severity on one scale across genes and
folds in human constraint (popEVE did this here) on the variants the paralog evidence prioritised. If the
identical change is severe in the teacher and moderate in the target, the residue-level transfer is sound
but the disease-level transfer is not, and the write-up must say so.

## Step 7: say what would falsify it

Write down, before anyone asks, the observation that would kill the hypothesis. In this project that was:
a patient with a motor-domain missense variant and a phenotype would confirm it; an accumulating pattern
of deletions and truncating variants with phenotypes, and missense tolerance in population data, argues
the disease mechanism is dosage instead.

## Step 8: appraise every human observation before it enters the argument

Steps 1 to 7 produce candidates. A candidate becomes a claim only when a real person carries the variant,
and that is where this project nearly went wrong. Human observations arrive as database rows that look
authoritative and are not. Put each one through the same six questions before it is allowed into a
document, a figure, or a sentence to the user.

1. **What is the inheritance?** De novo in an affected child and maternal from an unaffected mother are
   opposite evidence. If the field is absent, say what the observation costs without it and treat getting
   it as the highest-value open action. Here three of the four loss-of-function events with known
   inheritance turned out maternal, which reversed the reading of the whole axis.
2. **Is the lesion really single-gene?** A deletion smaller than the gene and lying inside its boundaries
   is single-gene by geometry and needs no assertion from anyone. Anything larger needs its flanking
   content resolved on Ensembl in both assemblies, plus its regulatory content, plus the constraint of the
   neighbours it could affect at a distance.
3. **What else was looked for in that patient?** A microarray record says nothing about point mutations,
   and a count of deposited variants is not a negative exome.
4. **What does this submitter or database do by default?** Compute the baseline before reading any count
   as a signal, and compute the expected count before reading any absence as one. This is the mistake that
   cost this project its most exciting-looking lead.
5. **What is the phenotype's base rate?** A feature present in a few percent of all imaged children adds
   nothing to specificity, however striking it looks inside a phenotype list.
6. **Has anyone replicated it?** Fetch the citing papers. A case report that has sat five years without a
   second case is telling you something.

Then, when a finding starts sounding like a discovery, hand it to an agent whose instructions are to break
it, with the controls named in the prompt. Both audits run that way in this project changed the answer.

## Porting the method to a new pair

A pair is worth the effort when the clinically rich gene has enough curated pathogenic missense to
calibrate on (roughly 30 or more), the sibling has almost none, the shared functional domain is
identifiable and well conserved, the pathogenic variants cluster in that shared domain, and the sibling is
constrained enough that a disease is plausible at all. Check the sibling's missense constraint
specifically: a gene that tolerates missense will not have a missense disease, no matter how well the
alignment works.

`scripts/paralog_pairs_scan.py` in the project implements the screen; `docs/olcekleme_fizibilite.md`
records what it found.
