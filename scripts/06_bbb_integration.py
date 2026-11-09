"""
06_bbb_integration.py - attach blood-brain-barrier (BBB) evidence to every drug and
compare three prioritisation strategies.

INPUTS : data/raw/B3DB_classification.tsv, data/raw/B3DB_regression.tsv
         data/processed/drugs.csv                (from script 03; has InChIKeys)
         results/tables/proximity_t<T>.csv       (from script 05; any thresholds present)
OUTPUTS: data/processed/b3db_compounds.csv       one row per B3DB compound (aggregated)
         results/tables/drug_bbb_evidence.csv    BBB evidence for each drug
         results/tables/proximity_bbb_t<T>.csv   proximity table + BBB columns + strategy ranks
         results/tables/bbb_strategy_summary.csv how the three strategies compare

HOW DRUGS ARE MATCHED TO B3DB (by chemical structure, NOT by name)
  Both ChEMBL and B3DB give a structure. We turn B3DB's InChI string into an InChIKey
  and join on its first 14 characters (the molecular-connectivity layer). That ignores
  stereochemistry, protonation and salt differences, so e.g. (S)-donepezil and donepezil
  match, while two different molecules with similar names never do.
  Drug names are used only as a visible sanity check (columns `b3db_names`).

B3DB HAS CONFLICTING RECORDS  (the same molecule can be BBB+ in one source, BBB- in another)
  Per compound we therefore count the records and set:
     BBB+          all records say BBB+
     BBB-          all records say BBB-
     conflicting   both labels occur
     no_evidence   compound not in B3DB   <-- NOT the same as "cannot cross the BBB"

THE THREE STRATEGIES (rank 1 = best)
  A  proximity only                          every drug, ranked by z-score
  B  proximity, BBB-negative drugs removed   keeps BBB+, conflicting and no_evidence
  C  proximity, BBB+ evidence required       keeps only drugs with BBB+ evidence
"""
import sys
from pathlib import Path

import pandas as pd
from rdkit import RDLogger
from rdkit.Chem import inchi

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

RDLogger.DisableLog("rdApp.*")      # silence harmless RDKit parsing warnings
TOP_N = (10, 20, 50)


def inchikey14(inchi_string):
    """First 14 characters of the InChIKey generated from an InChI string (or None)."""
    try:
        key = inchi.InchiToInchiKey(inchi_string)
        return key[:14] if key else None
    except Exception:
        return None


def build_b3db_compounds() -> pd.DataFrame:
    """Aggregate B3DB to one row per compound (InChIKey14)."""
    cls = pd.read_csv(cfg.RAW / "B3DB_classification.tsv", sep="\t")
    reg = pd.read_csv(cfg.RAW / "B3DB_regression.tsv", sep="\t")
    cls["inchikey14"] = cls["Inchi"].map(inchikey14)
    reg["inchikey14"] = reg["Inchi"].map(inchikey14)
    print(f"B3DB classification: {len(cls)} records, {cls['inchikey14'].isna().sum()} without a usable InChIKey")

    cls = cls.dropna(subset=["inchikey14"])
    cls["is_pos"] = (cls["BBB+/BBB-"] == "BBB+")
    agg = cls.groupby("inchikey14").agg(
        n_records=("is_pos", "size"),
        n_bbb_pos=("is_pos", "sum"),
        b3db_names=("compound_name", lambda s: "; ".join(sorted(set(map(str, s)))[:3])),
        b3db_groups=("group", lambda s: "".join(sorted(set(map(str, s))))),
    ).reset_index()
    agg["n_bbb_neg"] = agg["n_records"] - agg["n_bbb_pos"]
    agg["bbb_label"] = "conflicting"
    agg.loc[agg["n_bbb_neg"] == 0, "bbb_label"] = "BBB+"
    agg.loc[agg["n_bbb_pos"] == 0, "bbb_label"] = "BBB-"

    # logBB (numeric brain/blood ratio) exists only for ~1000 compounds: median over records
    logbb = (reg.dropna(subset=["inchikey14", "logBB"]).groupby("inchikey14")["logBB"]
             .median().rename("logBB_median").reset_index())
    agg = agg.merge(logbb, on="inchikey14", how="left")
    agg.to_csv(cfg.PROCESSED / "b3db_compounds.csv", index=False)
    print(f"{len(agg)} unique B3DB compounds:", agg["bbb_label"].value_counts().to_dict())
    return agg


def drug_evidence(b3db: pd.DataFrame) -> pd.DataFrame:
    """Join each drug to its B3DB evidence; unmatched drugs get 'no_evidence'."""
    drugs = pd.read_csv(cfg.PROCESSED / "drugs.csv")
    ev = drugs.merge(b3db, on="inchikey14", how="left")
    ev["bbb_label"] = ev["bbb_label"].fillna("no_evidence")
    ev.to_csv(cfg.TABLES / "drug_bbb_evidence.csv", index=False)
    print(f"\nBBB evidence for {len(ev)} drugs:\n{ev['bbb_label'].value_counts().to_string()}")
    return ev


def add_strategy_ranks(prox: pd.DataFrame) -> pd.DataFrame:
    """Rank drugs by z-score (more negative = closer) under each strategy; NaN = excluded."""
    prox = prox.copy()
    keep = {"rank_A_proximity_only": prox["bbb_label"].notna(),
            "rank_B_exclude_BBB_negative": prox["bbb_label"] != "BBB-",
            "rank_C_require_BBB_positive": prox["bbb_label"] == "BBB+"}
    for col, mask in keep.items():
        prox[col] = prox.loc[mask, "z_score"].rank(method="min")
    return prox


def top_overlap(prox, col, n):
    """How many of strategy A's top-n drugs are still in the top-n of another strategy."""
    a = set(prox.nsmallest(n, "rank_A_proximity_only")["drug_id"])
    b = set(prox.dropna(subset=[col]).nsmallest(n, col)["drug_id"])
    return len(a & b)


def main():
    cfg.ensure_dirs()
    ev = drug_evidence(build_b3db_compounds())
    summary = []
    for t in cfg.THRESHOLDS:
        f = cfg.TABLES / f"proximity_t{t}.csv"
        if not f.exists():
            print(f"(no proximity table for threshold {t} yet - skipping)")
            continue
        prox = pd.read_csv(f).merge(
            ev[["drug_id", "bbb_label", "n_records", "n_bbb_pos", "n_bbb_neg",
                "logBB_median", "b3db_names"]], on="drug_id", how="left")
        prox["bbb_label"] = prox["bbb_label"].fillna("no_evidence")
        prox = add_strategy_ranks(prox).sort_values("rank_A_proximity_only")
        prox.to_csv(cfg.TABLES / f"proximity_bbb_t{t}.csv", index=False)
        for n in TOP_N:
            summary.append({"threshold": t, "top_n": n,
                            "B_overlap_with_A": top_overlap(prox, "rank_B_exclude_BBB_negative", n),
                            "C_overlap_with_A": top_overlap(prox, "rank_C_require_BBB_positive", n),
                            "A_top_n_with_BBB_pos": int((prox.nsmallest(n, "rank_A_proximity_only")
                                                         ["bbb_label"] == "BBB+").sum()),
                            "A_top_n_with_BBB_neg": int((prox.nsmallest(n, "rank_A_proximity_only")
                                                         ["bbb_label"] == "BBB-").sum()),
                            "A_top_n_no_evidence": int((prox.nsmallest(n, "rank_A_proximity_only")
                                                        ["bbb_label"] == "no_evidence").sum())})
    if summary:
        s = pd.DataFrame(summary)
        s.to_csv(cfg.TABLES / "bbb_strategy_summary.csv", index=False)
        print("\nStrategy comparison:\n", s.to_string(index=False))


if __name__ == "__main__":
    main()
