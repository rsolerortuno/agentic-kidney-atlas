"""External enrichment: Open Targets associations for all genes; target details, STRING and UniProt
for a candidate set (pre-scored without tractability). Results cached under data/cache/ and summarized
in data/processed/external_*.parquet / .json."""
from __future__ import annotations

import json

import pandas as pd

from .candidates import omics_gene_table
from .config import load
from .external import ot_disease_associations, ot_target, string_partners, uniprot
from .integrate import compute
from .paths import PROCESSED, RAW

N_ENRICH = 300      # top genes by convergence score without the tractability stream
N_OT_ALWAYS = 100   # top Open Targets ccRCC targets always enriched (druggability of known disease genes)


def run() -> dict:
    table = omics_gene_table()
    ot = pd.DataFrame(ot_disease_associations())
    ot = ot.sort_values("score", ascending=False).drop_duplicates("symbol")  # a few symbols map to 2 Ensembl ids
    ot.to_parquet(PROCESSED / "external_ot_associations.parquet", index=False)
    # include OT-only genes (e.g. mutated drivers absent from bulk DE) in the universe
    extra = sorted(set(ot["symbol"]) - set(table["gene"]))
    if extra:
        table = pd.concat([table, pd.DataFrame({"gene": extra})], ignore_index=True)
        for col in ("rna_sig", "protein_sig", "phospho_sig", "concordant"):
            table[col] = table[col].fillna(False).astype(bool)
    table.to_parquet(PROCESSED / "gene_universe.parquet", index=False)
    # same scoring as integrate.py, minus the tractability stream (which is what we are about to fetch)
    sw = load()["convergence"]["stream_weights"]
    scored = compute()
    scored["pre_score"] = scored["score"] - sw["tractability"] * scored["c_tractability"]
    chosen = scored.sort_values("pre_score", ascending=False)["gene"].head(N_ENRICH).tolist()
    chosen = list(dict.fromkeys(chosen + ot["symbol"].head(N_OT_ALWAYS).tolist()))

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
    return {"ot_associations": len(ot), "enriched_genes": len(details),
            "with_drugs": sum(1 for d in details.values() if d["open_targets"] and d["open_targets"]["drugs"])}


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
