"""External enrichment: Open Targets associations for all genes; target details, STRING and UniProt
for a candidate set (pre-scored without tractability). Results cached under data/cache/ and summarized
in data/processed/external_*.parquet / .json."""
from __future__ import annotations

import json

import pandas as pd

from .candidates import omics_gene_table
from .config import load
from .external import ot_disease_associations, ot_target, string_partners, uniprot
from .literature_entities import gene_literature_counts
from .paths import PROCESSED, RAW

N_ENRICH = 300


def pre_score(table: pd.DataFrame, ot: pd.DataFrame, lit: pd.Series) -> pd.Series:
    w = load()["convergence"]["weights"]
    s = (w["rna"] * table["rna_sig"] + w["protein"] * table["protein_sig"] + w["phospho"] * table["phospho_sig"]
         + w["concordance"] * table["concordant"] + w["cnv"] * table["cnv_freq"].fillna(0)
         + w["open_targets"] * table["gene"].map(ot.set_index("symbol")["score"]).fillna(0)
         + w["literature"] * (table["gene"].map(lit).fillna(0) / 3).clip(upper=1))
    return s


def run() -> dict:
    table = omics_gene_table()
    ot = pd.DataFrame(ot_disease_associations())
    ot = ot.sort_values("score", ascending=False).drop_duplicates("symbol")  # a few symbols map to 2 Ensembl ids
    ot.to_parquet(PROCESSED / "external_ot_associations.parquet", index=False)
    lit = gene_literature_counts()
    # include OT-only genes (e.g. mutated drivers absent from bulk DE) in the universe
    extra = sorted(set(ot["symbol"]) - set(table["gene"]))
    if extra:
        table = pd.concat([table, pd.DataFrame({"gene": extra})], ignore_index=True)
        for col in ("rna_sig", "protein_sig", "phospho_sig", "concordant"):
            table[col] = table[col].fillna(False).astype(bool)
    table["pre_score"] = pre_score(table, ot, lit)
    chosen = table.sort_values("pre_score", ascending=False)["gene"].head(N_ENRICH).tolist()

    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str)
    ens = hg.dropna().set_index("symbol")["ensembl_gene_id"].to_dict()
    ens.update(ot.set_index("symbol")["ensembl"].to_dict())
    details = {}
    for g in chosen:
        rec = {"gene": g, "ensembl": ens.get(g)}
        rec["open_targets"] = ot_target(ens[g]) if ens.get(g) else None
        rec["string"] = string_partners(g)
        rec["uniprot"] = uniprot(g)
        details[g] = rec
    (PROCESSED / "external_targets.json").write_text(json.dumps(details, indent=1))
    table.to_parquet(PROCESSED / "gene_universe.parquet", index=False)
    return {"ot_associations": len(ot), "enriched_genes": len(details),
            "with_drugs": sum(1 for d in details.values() if d["open_targets"] and d["open_targets"]["drugs"])}


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
