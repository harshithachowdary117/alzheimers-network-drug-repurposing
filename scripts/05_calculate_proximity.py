"""
05_calculate_proximity.py - network proximity of every drug to the Alzheimer's gene set.

INPUTS : data/processed/ppi_t<T>.csv, ad_genes.csv, drug_targets.csv, drugs.csv
OUTPUT : results/tables/proximity_t<T>.csv   (one row per drug)

Method (Guney et al. 2016, implemented in NetMedPy)
  * d(drug, AD) = average, over the drug's targets, of the shortest-path distance
    to the nearest Alzheimer's gene  ("raw_amspl").
  * A null distribution is built by drawing random target sets with the SAME size
    and SIMILAR degrees (degree-matched), N_ITER times.
  * z_score = (d - mean_random) / sd_random.   MORE NEGATIVE = CLOSER than chance.
  * p_single_tail = fraction of random draws at least as close as the real drug.
  * q_value = Benjamini-Hochberg false-discovery-rate correction over all drugs.

IMPORTANT: NetMedPy uses multiprocessing with the 'spawn' start method (the default
on macOS/Windows), so everything below must run under `if __name__ == "__main__"`.

Usage
  python scripts/05_calculate_proximity.py --threshold 700
  python scripts/05_calculate_proximity.py --threshold 400 --n-iter 100 --n-procs 2   # quick test
"""
import argparse
import os
import sys
from pathlib import Path

import networkx as nx
import netmedpy
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

# the only properties NetMedPy's screening() can return
PROPERTIES = ["z_score", "p_value_single_tail", "p_value_double_tail", "raw_amspl"]


def load_graph(threshold: int) -> nx.Graph:
    edges = pd.read_csv(cfg.PROCESSED / f"ppi_t{threshold}.csv")
    return nx.from_pandas_edgelist(edges, "source", "target")


def run_proximity(threshold, distance=cfg.DISTANCE, null_model=cfg.NULL_MODEL,
                  n_iter=cfg.N_ITER, min_targets=cfg.MIN_TARGETS, n_procs=None,
                  ad_subset=None, dmat=None):
    """Compute proximity for all drugs on one network and return a results DataFrame.

    ad_subset : optional set of gene names; if given, only Alzheimer's genes in this set
                are used (script 07 uses it for the 'common gene set' variant).
    dmat      : optional precomputed NetMedPy distance matrix for this network, so that
                script 07 does not recompute the slow all-pairs distances for every variant.
    """
    n_procs = n_procs or max(1, (os.cpu_count() or 2) - 1)
    G = load_graph(threshold)
    nodes = set(G.nodes)

    # --- map Alzheimer's genes and drug targets onto the network (keep only nodes that exist)
    ad_all = set(pd.read_csv(cfg.PROCESSED / "ad_genes.csv")["gene_symbol"])
    ad = ad_all & nodes
    if ad_subset is not None:
        ad = ad & set(ad_subset)
    dtab = pd.read_csv(cfg.PROCESSED / "drug_targets.csv")
    drugs = pd.read_csv(cfg.PROCESSED / "drugs.csv").set_index("drug_id")
    targets_total = dtab.groupby("drug_id")["gene_symbol"].apply(set)
    targets_in = targets_total.apply(lambda s: s & nodes)
    targets_in = targets_in[targets_in.apply(len) >= min_targets]
    print(f"[t={threshold}] network {G.number_of_nodes():,} nodes | "
          f"AD genes in network {len(ad)}/{len(ad_all)} | drugs analysed {len(targets_in)}")

    # --- all-pairs distance matrix (the slow, memory-heavy step), then screen every drug
    if dmat is None:
        dmat = netmedpy.all_pair_distances(G, distance=distance, n_processors=n_procs)
    screen = netmedpy.screening(targets_in.to_dict(), {"AD": ad}, G, dmat,
                                score="proximity", properties=PROPERTIES,
                                null_model=null_model, n_iter=n_iter, n_procs=n_procs)

    # --- assemble one tidy table (screening returns one DataFrame per property: drugs x diseases)
    res = pd.DataFrame({p: screen[p]["AD"] for p in PROPERTIES})
    res = res.rename(columns={"p_value_single_tail": "p_single", "p_value_double_tail": "p_double",
                              "raw_amspl": "distance"})
    res.index.name = "drug_id"
    # A z-score is undefined (NaN/inf) when the random draws have zero spread. Never crash:
    # flag those drugs, leave them unranked, and report how many there are.
    res["z_score"] = res["z_score"].astype(float).replace([np.inf, -np.inf], np.nan)
    res["z_valid"] = res["z_score"].notna()
    if (~res["z_valid"]).any():
        print(f"WARNING: {(~res['z_valid']).sum()} drugs have an undefined z-score and are left unranked")
    # BH false-discovery-rate correction over the drugs that have a valid p-value
    res["q_value"] = np.nan
    ok = res["p_single"].notna()
    res.loc[ok, "q_value"] = false_discovery_control(res.loc[ok, "p_single"].astype(float), method="bh")
    res["n_targets_total"] = targets_total.reindex(res.index).apply(len)
    res["n_targets_in_network"] = targets_in.reindex(res.index).apply(len)
    res["rank"] = res["z_score"].rank(method="min").astype("Int64")  # 1 = closest to AD genes
    # DIRECT HIT = at least one drug target is itself an Alzheimer's gene (distance 0 for
    # single-target drugs). For these, proximity only re-states "the target is a disease
    # gene" and its z-score is not comparable with that of multi-target drugs.
    res["n_targets_in_ad"] = [len(targets_in[d] & ad) for d in res.index]
    res["direct_hit"] = res["n_targets_in_ad"] > 0
    res["n_ad_genes_used"] = len(ad)
    res["threshold"] = threshold
    res = res.join(drugs[["drug_name", "inchikey14"]]).sort_values("rank")
    return res.reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--threshold", type=int, required=True, help="STRING score threshold")
    ap.add_argument("--distance", default=cfg.DISTANCE)
    ap.add_argument("--null-model", default=cfg.NULL_MODEL)
    ap.add_argument("--n-iter", type=int, default=cfg.N_ITER)
    ap.add_argument("--min-targets", type=int, default=cfg.MIN_TARGETS)
    ap.add_argument("--n-procs", type=int, default=None)
    ap.add_argument("--out", default=None, help="output CSV (default results/tables/proximity_t<T>.csv)")
    a = ap.parse_args()

    cfg.ensure_dirs()
    res = run_proximity(a.threshold, a.distance, a.null_model, a.n_iter, a.min_targets, a.n_procs)
    out = Path(a.out) if a.out else cfg.TABLES / f"proximity_t{a.threshold}.csv"
    res.to_csv(out, index=False)
    print(f"Saved {len(res)} drugs -> {out}")
    print(res[["rank", "drug_name", "n_targets_in_network", "distance", "z_score", "q_value"]]
          .head(15).to_string(index=False))


if __name__ == "__main__":
    main()
