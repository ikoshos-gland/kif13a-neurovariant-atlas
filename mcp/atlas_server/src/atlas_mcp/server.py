"""kif13a-atlas MCP server: the project's data layer.

Third-party MCP servers answer exploratory questions; this server pulls the
data the pipeline depends on and writes it to data/raw with provenance.
"""

from __future__ import annotations

import logging
import sys
from typing import Any

from mcp.server.fastmcp import FastMCP

from . import core, expression, genome_models, structure, variants

logging.basicConfig(stream=sys.stderr, level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
core.load_dotenv()

INSTRUCTIONS = """Data layer of the KIF13A NeuroVariant Atlas project.
Rules: protein positions are MANE Select numbering (see atlas_targets). KIF13A MANE = UniProt Q9H1H9, 1805 aa.
KIF1A MANE = Q12756-3 (1791 aa) while UniProt canonical Q12756 is 1690 aa; positions 1-393 agree, later ones shift.
Tools with save=True write files under data/raw (or data/interim) and return paths; prefer them over pasting
large tables into the conversation. Use data_inventory to see what has already been pulled."""

mcp = FastMCP("kif13a-atlas", instructions=INSTRUCTIONS)


@mcp.tool()
async def atlas_targets(gene: str | None = None) -> dict[str, Any]:
    """Verified identifiers for the project genes (KIF13A, KIF1A, KIF13B, KIF1B, KIF1C): Ensembl gene, MANE
    transcript, RefSeq, UniProt, protein length, isoform numbering traps and reference structures."""
    return core.target(gene) if gene else core.load_targets()


@mcp.tool()
async def data_inventory(limit: int = 200) -> list[dict[str, Any]]:
    """List files already pulled into the project (latest provenance record per file)."""
    return core.inventory(limit)


@mcp.tool()
async def gnomad_constraint(gene: str, save: bool = True) -> dict[str, Any]:
    """gnomAD v4.1 gene constraint: pLI, LOEUF (oe_lof_upper), missense Z, oe_mis with bounds."""
    return await variants.gnomad_constraint(gene, save)


@mcp.tool()
async def gnomad_variants(gene: str, consequence: str = "missense", save: bool = True) -> dict[str, Any]:
    """All gnomAD v4.1 (joint exome+genome) variants on the gene's MANE transcript with MANE-numbered HGVSp,
    allele counts/frequencies, filters and gnomAD in-silico scores (CADD, REVEL, SpliceAI, Pangolin, phyloP).
    consequence: 'missense' | 'coding' | 'all'. Saves a TSV; returns counts and AF bins."""
    return await variants.gnomad_variants(gene, consequence, save)


@mcp.tool()
async def clinvar_variants(gene: str, save: bool = True) -> dict[str, Any]:
    """Complete ClinVar pull for a single gene (NCBI E-utilities, no row cap): germline classification, review
    status and stars, traits, consequence, HGVS on the title transcript with a flag for MANE, GRCh38 SPDI."""
    return await variants.clinvar_variants(gene, save)


@mcp.tool()
async def vep_annotate(variants_list: list[str], label: str = "batch", save: bool = True) -> dict[str, Any]:
    """Annotate variants with Ensembl VEP on the MANE transcript plus AlphaMissense, CADD, SpliceAI, REVEL,
    dbNSFP ESM1b/PrimateAI/MetaRNN and gnomAD AF. Accepts HGVS ('ENST00000259711.11:c.1048C>T'),
    'chr6-17764000-A-G' style or rsIDs; any number (batched by 200)."""
    return await variants.vep_annotate(variants_list, label, save)


@mcp.tool()
async def alphamissense(
    uniprot: str, positions: str | None = None, variants_list: list[str] | None = None, save: bool = True
) -> dict[str, Any]:
    """AlphaMissense scores for every possible substitution of a UniProt canonical protein (from AlphaFold DB).
    positions like '300-360,410' gives per-position mean; variants_list like ['R350W'] gives exact scores.
    Canonical numbering only: for KIF1A beyond residue 393 use vep_annotate instead."""
    return await variants.alphamissense(uniprot, positions, variants_list, save)


@mcp.tool()
async def cadd_region(chrom: str, start: int, end: int, version: str = "GRCh38-v1.7", save: bool = True) -> dict[str, Any]:
    """CADD PHRED/raw scores for every possible SNV in a GRCh38 region (max 2 kb)."""
    return await variants.cadd_region(chrom, start, end, version, save)


@mcp.tool()
async def evo2_score(
    variants_list: list[str], window: int = 8192, both_strands: bool = False, label: str = "evo2", save: bool = True
) -> dict[str, Any]:
    """Evo 2 40B zero-shot score for GRCh38 SNVs ('2-240773246-G-A'): delta log-likelihood of the alt vs ref
    genomic window centred on the variant (NVIDIA NIM; spends hosted credits, max 50 variants per call).
    More negative = more disruptive. both_strands=True also scores the reverse complement (2x calls)."""
    return await genome_models.evo2_score(variants_list, window, both_strands, label, save)


@mcp.tool()
async def alphafold_model(uniprot: str, include_pae: bool = True, save: bool = True) -> dict[str, Any]:
    """Download the latest AlphaFold DB model (mmCIF, per-residue pLDDT, optional PAE) into data/raw/alphafold and
    report confidence fractions and pLDDT>=70 segments. Works for isoform accessions such as Q12756-3."""
    return await structure.alphafold_model(uniprot, include_pae, save)


@mcp.tool()
async def uniprot_features(accession: str, save: bool = True) -> dict[str, Any]:
    """UniProt entry features grouped by type (domains, regions, coiled coils, motifs, binding sites,
    mutagenesis) with counts for modified residues and natural variants."""
    return await structure.uniprot_features(accession, save)


@mcp.tool()
async def pdb_structures(uniprot: str, download: list[str] | None = None, save: bool = True) -> dict[str, Any]:
    """Experimental structures mapped to a UniProt accession (PDBe SIFTS) with titles, method, resolution and
    residue coverage; optionally download mmCIF files for listed PDB IDs."""
    return await structure.pdb_structures(uniprot, download, save)


@mcp.tool()
async def paralog_map(source_gene: str, target_gene: str, positions: str, save: bool = True) -> dict[str, Any]:
    """Map residue positions between two project genes on MANE proteins (e.g. KIF1A 350 -> KIF13A) with a quick
    BLOSUM62 global alignment; reports identity and local 15-residue window identity. Preliminary map."""
    return await structure.paralog_map(source_gene, target_gene, positions, save)


@mcp.tool()
async def gtex_expression(gene: str, dataset: str = "gtex_v10", save: bool = True) -> dict[str, Any]:
    """GTEx median TPM per tissue (default GTEx v10) with the brain tissues listed separately."""
    return await expression.gtex_expression(gene, dataset, save)


@mcp.tool()
async def hpa_summary(gene: str, save: bool = True) -> dict[str, Any]:
    """Human Protein Atlas summary: RNA tissue/single-cell/brain specificity, subcellular location, reliability."""
    return await expression.hpa_summary(gene, save)


@mcp.tool()
async def cellxgene_expression(
    gene: str, tissue: str | None = None, min_cells: int = 100, top_n: int = 25, save: bool = True
) -> dict[str, Any]:
    """CELLxGENE Discover aggregated human single-cell expression per tissue x cell type (Where's My Gene v2),
    ranked by mean expression x fraction expressing; tissue is a substring filter such as 'brain'."""
    return await expression.cellxgene_expression(gene, tissue, min_cells, top_n, save)


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
