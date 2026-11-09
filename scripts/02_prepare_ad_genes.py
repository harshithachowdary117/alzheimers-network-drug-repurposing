"""
02_prepare_ad_genes.py - turn the Open Targets download into the Alzheimer's gene set.

INPUT  : data/raw/opentargets_alzheimer_associations.json
OUTPUTS: data/processed/ad_genes_all.csv        every associated gene + per-datatype scores
         data/processed/ad_genes.csv            the selected gene set used downstream
         results/tables/ad_gene_set_sensitivity.csv   gene-set size at different cut-offs

Selection rule (see config.py)
  A gene is an "Alzheimer's gene" if its BEST datatype score, ignoring the
  datatypes in config.EXCLUDE_DATATYPES, is >= config.AD_MIN_SCORE.

Why exclude 'known_drug'?  That evidence comes from drugs already studied in
Alzheimer's. If it defined the disease genes, the proximity analysis would
partly rediscover those same drugs (circular reasoning).
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402


def load_associations() -> pd.DataFrame:
    """Flatten the JSON into one row per gene with one column per datatype score."""
    record = json.loads((cfg.RAW / "opentargets_alzheimer_associations.json").read_text())
    print(f"Disease: {record['disease_name']} ({record['disease_id']}), "
          f"retrieved {record['retrieved']}, {len(record['rows'])} associated targets")
    recs = []
    for row in record["rows"]:
        rec = {"ensembl_id": row["target"]["id"],
               "gene_symbol": cfg.norm_symbol(row["target"]["approvedSymbol"]),
               "overall_score": row["score"]}
        for d in row["datatypeScores"]:
            rec[f"dt_{d['id']}"] = d["score"]
        recs.append(rec)
    df = pd.DataFrame(recs)
    dt_cols = [c for c in df.columns if c.startswith("dt_")]
    df[dt_cols] = df[dt_cols].fillna(0.0)          # missing datatype = no evidence
    return df


def main():
    cfg.ensure_dirs()
    df = load_associations()

    # Best score among the datatypes we are allowed to use.
    keep_cols = [c for c in df.columns
                 if c.startswith("dt_") and c[3:] in cfg.AD_DATATYPES]
    print("Datatypes used for selection:", [c[3:] for c in keep_cols])
    df["selection_score"] = df[keep_cols].max(axis=1)
    df = df.drop_duplicates("gene_symbol").sort_values("selection_score", ascending=False)
    df.to_csv(cfg.PROCESSED / "ad_genes_all.csv", index=False)

    # How sensitive is the gene-set size to the cut-off? (guides the final choice)
    rows = []
    for cutoff in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
        rows.append({"min_selection_score": cutoff,
                     "n_genes": int((df["selection_score"] >= cutoff).sum())})
    sens = pd.DataFrame(rows)
    sens.to_csv(cfg.TABLES / "ad_gene_set_sensitivity.csv", index=False)
    print("\nGene-set size by cut-off:\n", sens.to_string(index=False))

    selected = df[df["selection_score"] >= cfg.AD_MIN_SCORE]
    selected[["gene_symbol", "ensembl_id", "overall_score", "selection_score"]].to_csv(
        cfg.PROCESSED / "ad_genes.csv", index=False)
    print(f"\nSelected {len(selected)} Alzheimer's genes (cut-off {cfg.AD_MIN_SCORE}); "
          f"top 10: {', '.join(selected.gene_symbol.head(10))}")
    print("Sanity check - are the classic AD genes present?",
          {g: g in set(selected.gene_symbol) for g in ("APP", "PSEN1", "PSEN2", "APOE", "MAPT")})


if __name__ == "__main__":
    main()
