"""
07_threshold_sweep.py - the project's additional experiment: repeat the proximity analysis
on STRING networks of different confidence (config.THRESHOLDS = 400 / 700 / 900).

INPUTS : data/processed/ppi_t<T>.csv (script 04), ad_genes.csv, drug_targets.csv, drugs.csv
OUTPUTS: results/tables/proximity_t<T>.csv          'native' variant
         results/tables/proximity_common_t<T>.csv   'common-gene' variant
         results/tables/common_ad_genes.csv

TWO VARIANTS - why both?
  native  Each network uses the Alzheimer's genes that exist in it. Higher thresholds
          have fewer of them (coverage falls from ~84% to ~44%), so a drug's rank can change
          because the DISEASE SET changed, not because the edges were filtered.
  common  All networks use the SAME Alzheimer's genes: those present in every network's
          largest connected component. Differences between thresholds now come from the
          network edges only.
  If both variants agree, the conclusion about threshold sensitivity is trustworthy.

The slow all-pairs distance matrix is computed once per network and shared by both variants.

Usage
  python scripts/07_threshold_sweep.py                     # full run (1000 iterations)
  python scripts/07_threshold_sweep.py --reuse             # skip result files that already exist
  python scripts/07_threshold_sweep.py --n-iter 100 --n-procs 2    # quick test
"""
import argparse
import importlib
import os
import sys
from pathlib import Path

import netmedpy
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

# module names starting with a digit cannot be used in an `import` statement
prox = importlib.import_module("05_calculate_proximity")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--n-iter", type=int, default=cfg.N_ITER)
    ap.add_argument("--n-procs", type=int, default=None)
    ap.add_argument("--variants", nargs="+", default=["native", "common"], choices=["native", "common"])
    ap.add_argument("--reuse", action="store_true", help="skip outputs that already exist")
    a = ap.parse_args()
    n_procs = a.n_procs or max(1, (os.cpu_count() or 2) - 1)
    cfg.ensure_dirs()

    # genes present in the largest connected component of EVERY network
    node_sets = [set(prox.load_graph(t).nodes) for t in cfg.THRESHOLDS]
    ad_all = set(pd.read_csv(cfg.PROCESSED / "ad_genes.csv")["gene_symbol"])
    common = ad_all & set.intersection(*node_sets)
    pd.Series(sorted(common), name="gene_symbol").to_csv(cfg.TABLES / "common_ad_genes.csv", index=False)
    print(f"Alzheimer's genes present in all {len(cfg.THRESHOLDS)} networks: {len(common)} of {len(ad_all)}")

    plan = {"native": (None, cfg.TABLES / "proximity_t{t}.csv"),
            "common": (common, cfg.TABLES / "proximity_common_t{t}.csv")}
    for t in cfg.THRESHOLDS:
        todo = {v: str(plan[v][1]).format(t=t) for v in a.variants}
        if a.reuse:
            todo = {v: f for v, f in todo.items() if not Path(f).exists()}
        if not todo:
            print(f"[t={t}] all requested outputs exist - skipped")
            continue
        G = prox.load_graph(t)
        dmat = netmedpy.all_pair_distances(G, distance=cfg.DISTANCE, n_processors=n_procs)
        for v, out in todo.items():
            res = prox.run_proximity(t, n_iter=a.n_iter, n_procs=n_procs,
                                     ad_subset=plan[v][0], dmat=dmat)
            res["variant"] = v
            res.to_csv(out, index=False)
            print(f"[t={t}] {v}: {len(res)} drugs -> {out}")
        del dmat                                    # free ~1 GB before the next network
    print("\nSweep finished. Next: python scripts/06_bbb_integration.py, then 08_compare_results.py")


if __name__ == "__main__":
    main()
