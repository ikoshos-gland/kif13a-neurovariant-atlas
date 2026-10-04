"""Expression sources: GTEx (bulk tissue), Human Protein Atlas, CELLxGENE aggregated single-cell.

CELLxGENE Census itself needs tiledbsoma, which has no Windows build; the
Where's My Gene v2 API used here serves the same corpus pre-aggregated per
tissue x cell type, which is enough for gene-level cell-type ranking.
Cell-level analyses (donor-aware tests) belong on TRUBA.
"""

from __future__ import annotations

from typing import Any

from . import core

GTEX = "https://gtexportal.org/api/v2"
GTEX_GENCODE = {"gtex_v10": "v39", "gtex_v8": "v26"}
HPA = "https://www.proteinatlas.org"
WMG = "https://api.cellxgene.cziscience.com/wmg/v2"


def _ensembl_gene(gene: str) -> tuple[str, str]:
    if gene.upper().startswith("ENSG"):
        return gene, gene
    t = core.target(gene)
    return t["symbol"], t["ensembl_gene"]


# --------------------------------------------------------------------------- GTEx


async def gtex_expression(gene: str, dataset: str = "gtex_v10", save: bool = True) -> dict[str, Any]:
    if dataset not in GTEX_GENCODE:
        raise ValueError(f"dataset must be one of {sorted(GTEX_GENCODE)}")
    symbol, _ = _ensembl_gene(gene)
    ref = await core.get_json(f"{GTEX}/reference/gene", params={"geneId": symbol, "gencodeVersion": GTEX_GENCODE[dataset]})
    if not ref.get("data"):
        raise ValueError(f"GTEx has no gene {symbol} in GENCODE {GTEX_GENCODE[dataset]}")
    gencode_id = ref["data"][0]["gencodeId"]
    data = await core.get_json(
        f"{GTEX}/expression/medianGeneExpression", params={"gencodeId": gencode_id, "datasetId": dataset}
    )
    rows = sorted(
        ({"tissue": d["tissueSiteDetailId"], "median_tpm": d["median"]} for d in data["data"]),
        key=lambda r: r["median_tpm"], reverse=True,
    )
    out: dict[str, Any] = {
        "gene": symbol, "gencode_id": gencode_id, "dataset": dataset, "unit": "median TPM", "n_tissues": len(rows),
        "top_tissues": rows[:10], "brain_tissues": [r for r in rows if r["tissue"].startswith("Brain")],
    }
    if save:
        out["saved"] = core.save_tsv(
            "gtex", f"{symbol}_{dataset}_median_tpm_{core.today()}.tsv", rows, source="GTEx Portal API v2",
            urls=[f"{GTEX}/expression/medianGeneExpression"], tool="gtex_expression",
            params={"gencodeId": gencode_id, "datasetId": dataset},
        )
    return out


# --------------------------------------------------------------------------- HPA

HPA_PREFIXES = ("RNA ", "Subcellular", "Reliability", "Protein class", "Disease involvement", "Secretome", "Evidence")


async def hpa_summary(gene: str, save: bool = True) -> dict[str, Any]:
    symbol, ensg = _ensembl_gene(gene)
    url = f"{HPA}/{ensg}.json"
    data = await core.get_json(url)
    compact = {
        k: v for k, v in data.items()
        if k.startswith(HPA_PREFIXES) and v not in (None, "", [], {})
    }
    out: dict[str, Any] = {"gene": data.get("Gene", symbol), "ensembl": ensg, "fields": compact}
    if save:
        out["saved"] = core.save_json(
            "hpa", f"{symbol}_hpa_{core.today()}.json", data, source="Human Protein Atlas JSON",
            urls=[url], tool="hpa_summary", params={"gene": gene},
        )
    return out


# --------------------------------------------------------------------------- CELLxGENE

_tissue_labels: dict[str, str] = {}


async def _tissues(organism: str) -> dict[str, str]:
    if not _tissue_labels:
        dims = await core.get_json(f"{WMG}/primary_filter_dimensions", timeout=300)
        for item in dims["tissue_terms"].get(organism, []):
            _tissue_labels.update(item)
    return _tissue_labels


async def cellxgene_expression(
    gene: str, tissue: str | None = None, min_cells: int = 100, top_n: int = 25, save: bool = True
) -> dict[str, Any]:
    organism = "NCBITaxon:9606"
    symbol, ensg = _ensembl_gene(gene)
    body = {
        "filter": {
            "dataset_ids": [], "development_stage_ontology_term_ids": [], "disease_ontology_term_ids": [],
            "gene_ontology_term_ids": [ensg], "organism_ontology_term_id": organism,
            "self_reported_ethnicity_ontology_term_ids": [], "sex_ontology_term_ids": [], "publication_citations": [],
        },
        "is_rollup": True,
    }
    data = await core.post_json(f"{WMG}/query", body, headers={"Content-Type": "application/json"}, timeout=300)
    tissues = await _tissues(organism)
    labels = data["term_id_labels"]["cell_types"]
    rows = []
    for tissue_id, cell_types in data["expression_summary"].get(ensg, {}).items():
        for ct_id, agg in cell_types.items():
            a = agg.get("aggregated", {})
            lab = labels.get(tissue_id, {}).get(ct_id, {}).get("aggregated", {})
            total = lab.get("total_count")
            rows.append({
                "tissue_id": tissue_id, "tissue": tissues.get(tissue_id, tissue_id),
                "cell_type_id": ct_id, "cell_type": lab.get("name", ct_id),
                "mean_expression_me": a.get("me"), "n_expressing": a.get("n"), "total_cells": total,
                "fraction_expressing_pc": a.get("pc"),
                "fraction_check_n_over_total": round(a["n"] / total, 4) if total and a.get("n") is not None else None,
            })
    pool = [r for r in rows if (r["total_cells"] or 0) >= min_cells and r["cell_type_id"] != "CL:0000000"]
    if tissue:
        pool = [r for r in pool if tissue.lower() in r["tissue"].lower()]
    ranked = sorted(pool, key=lambda r: (r["mean_expression_me"] or 0) * (r["fraction_expressing_pc"] or 0), reverse=True)
    out: dict[str, Any] = {
        "gene": symbol, "ensembl": ensg, "snapshot_id": data.get("snapshot_id"),
        "fields": "me = WMG mean expression; pc = fraction of cells of this type expressing the gene (0-1 in "
                  "practice; n/total check column included); ranking = me x pc",
        "n_tissue_celltype_pairs": len(rows), "tissue_filter": tissue, "min_cells": min_cells,
        "top": [{k: r[k] for k in ("tissue", "cell_type", "mean_expression_me", "fraction_expressing_pc", "total_cells")}
                for r in ranked[:top_n]],
    }
    if save:
        out["saved"] = core.save_tsv(
            "cellxgene", f"{symbol}_cellxgene_wmg_{data.get('snapshot_id')}.tsv", rows,
            source="CELLxGENE Discover Where's My Gene v2", urls=[f"{WMG}/query"], tool="cellxgene_expression",
            params={"gene": ensg, "organism": organism, "is_rollup": True}, upstream_snapshot=data.get("snapshot_id"),
        )
    return out
