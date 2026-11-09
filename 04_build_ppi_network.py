"""
04_build_ppi_network.py - build one PPI network per STRING confidence threshold
and report how well the Alzheimer's genes and drug targets map onto it.

INPUTS : data/raw/string_physical_links.txt.gz, data/raw/string_protein_info.txt.gz
         data/processed/ad_genes.csv and drug_targets.csv   (optional, for the coverage QC)
OUTPUTS: data/processed/ppi_t<threshold>.csv     edge list (source,target) of the LCC
         results/tables/network_qc.csv           size + coverage of every network

Steps for each threshold T in config.THRESHOLDS
  1. keep STRING physical edges with combined_score >= T
  2. relabel proteins (ENSP ids) with gene names
  3. keep the largest connected component (LCC) - shortest-path distances are
     only defined inside one connected component
  4. save the edges and record network + coverage statistics

This is the quality-control step: proximity is only meaningful if most Alzheimer's
genes and drug targets are actually present in the network.
"""
import sys
from pathlib import Path

import networkx as nx
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402


def load_string_edges() -> pd.DataFrame:
    """Read STRING links and translate protein ids to upper-case gene names."""
    links = pd.read_csv(cfg.RAW / "string_physical_links.txt.gz", sep=" ", compression="gzip")
    info = pd.read_csv(cfg.RAW / "string_protein_info.txt.gz", sep="\t", compression="gzip")
    id2name = dict(zip(info.iloc[:, 0], info["preferred_name"]))   # first column = protein id
    links["a"] = links["protein1"].map(id2name).map(cfg.norm_symbol, na_action="ignore")
    links["b"] = links["protein2"].map(id2name).map(cfg.norm_symbol, na_action="ignore")
    links = links.dropna(subset=["a", "b"])
    links = links[links["a"] != links["b"]]                         # drop self-loops
    return links[["a", "b", "combined_score"]]


def build_lcc(links: pd.DataFrame, threshold: int):
    """Return (full graph, largest connected component) for one score threshold."""
    sub = links[links["combined_score"] >= threshold]
    full = nx.from_pandas_edgelist(sub, "a", "b")     # duplicates (A-B and B-A) collapse
    lcc_nodes = max(nx.connected_components(full), key=len)
    return full, full.subgraph(lcc_nodes).copy()


def main():
    cfg.ensure_dirs()
    links = load_string_edges()
    print(f"STRING: {len(links):,} directed edge rows, "
          f"{pd.concat([links.a, links.b]).nunique():,} genes")

    ad_f, dt_f = cfg.PROCESSED / "ad_genes.csv", cfg.PROCESSED / "drug_targets.csv"
    ad = set(pd.read_csv(ad_f)["gene_symbol"]) if ad_f.exists() else None
    dtab = pd.read_csv(dt_f) if dt_f.exists() else None

    qc = []
    for t in cfg.THRESHOLDS:
        full, lcc = build_lcc(links, t)
        nx.to_pandas_edgelist(lcc).to_csv(cfg.PROCESSED / f"ppi_t{t}.csv", index=False)
        row = {"threshold": t, "nodes_all": full.number_of_nodes(), "edges_all": full.number_of_edges(),
               "lcc_nodes": lcc.number_of_nodes(), "lcc_edges": lcc.number_of_edges(),
               "lcc_density": round(nx.density(lcc), 5)}
        nodes = set(lcc.nodes)
        if ad is not None:
            row.update(ad_genes_total=len(ad), ad_genes_in_lcc=len(ad & nodes),
                       ad_coverage=round(len(ad & nodes) / len(ad), 3))
        if dtab is not None:
            tg = set(dtab["gene_symbol"])
            in_net = dtab[dtab["gene_symbol"].isin(nodes)]
            row.update(target_genes_total=len(tg), target_genes_in_lcc=len(tg & nodes),
                       target_coverage=round(len(tg & nodes) / len(tg), 3),
                       drugs_total=dtab["drug_id"].nunique(),
                       drugs_with_target_in_lcc=in_net["drug_id"].nunique())
        qc.append(row)
        print(f"  threshold {t}: LCC {lcc.number_of_nodes():,} nodes / {lcc.number_of_edges():,} edges")

    qc = pd.DataFrame(qc)
    qc.to_csv(cfg.TABLES / "network_qc.csv", index=False)
    print("\nNetwork QC:\n", qc.T.to_string(header=False))


if __name__ == "__main__":
    main()
