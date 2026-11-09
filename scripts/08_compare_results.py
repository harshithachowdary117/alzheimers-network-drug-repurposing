"""
08_compare_results.py - how stable are the drug rankings across STRING thresholds, and
what do the BBB data change? Produces the final tables and figures.

INPUTS : results/tables/proximity_t<T>.csv            (native variant, script 07)
         results/tables/proximity_common_t<T>.csv     (common-gene variant, script 07; optional)
         results/tables/drug_bbb_evidence.csv         (script 06)
         results/tables/network_qc.csv                (script 04)
OUTPUTS (results/tables/)
         threshold_stability.csv        Spearman correlation + top-N overlap between thresholds
         significant_counts.csv         number of drugs with q < 0.05 per threshold
         rank_stability_all_drugs.csv   rank of every drug at every threshold + BBB + class
         top_candidates_indirect.csv    non-direct-hit drugs in the top 50 at >= 2 thresholds
         direct_hit_drugs.csv           drugs whose target is itself an Alzheimer's gene
         summary.txt                    the key numbers in plain text
(results/figures/)
         fig1_rank_scatter.png, fig2_topN_overlap.png, fig3_bbb_top50.png, fig4_network_coverage.png

DESIGN NOTES
  * Only drugs analysable at EVERY threshold are compared. Higher thresholds drop drugs whose
    targets fall out of the network; comparing different drug lists would confuse "the ranking
    changed" with "the drug disappeared".
  * Two views are reported: ALL drugs, and drugs WITHOUT a direct hit (no target that is itself
    an Alzheimer's gene). Direct-hit drugs have distance 0 and very large z-scores that are not
    comparable with multi-target drugs, so they should not drive the conclusions.
  * Ranking uses the z-score (more negative = closer to the Alzheimer's genes).
"""
import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                      # no display needed
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

TOP_NS = (10, 20, 50)
BBB_COLORS = {"BBB+": "#2a9d8f", "BBB-": "#e76f51", "conflicting": "#e9c46a", "no_evidence": "#b0b0b0"}
lines = []                                  # collected for summary.txt


def say(text=""):
    print(text)
    lines.append(text)


def load_variant(variant):
    """Return {threshold: DataFrame indexed by drug_id} for one variant, or None."""
    pattern = "proximity_t{t}.csv" if variant == "native" else "proximity_common_t{t}.csv"
    out = {}
    for t in cfg.THRESHOLDS:
        f = cfg.TABLES / pattern.format(t=t)
        if not f.exists():
            return None
        df = pd.read_csv(f).set_index("drug_id")
        out[t] = df[df["z_score"].notna()]
    return out


def common_drugs(tables):
    """Drugs with a valid z-score at every threshold."""
    return set.intersection(*[set(df.index) for df in tables.values()])


def direct_hit_drugs(tables):
    """Drugs flagged as direct hits at ANY threshold."""
    ids = set()
    for df in tables.values():
        if "direct_hit" in df:
            ids |= set(df.index[df["direct_hit"].astype(bool)])
    return ids


def stability(tables, drugs, variant, view):
    """Spearman correlation and top-N overlap between every pair of thresholds."""
    z = pd.DataFrame({t: tables[t].loc[sorted(drugs), "z_score"] for t in tables})
    rows = []
    for t1, t2 in itertools.combinations(tables, 2):
        row = {"variant": variant, "view": view, "threshold_a": t1, "threshold_b": t2,
               "n_drugs": len(z), "spearman_rho": spearmanr(z[t1], z[t2])[0]}
        for n in TOP_NS:
            top1, top2 = set(z[t1].nsmallest(n).index), set(z[t2].nsmallest(n).index)
            row[f"top{n}_overlap"] = len(top1 & top2)
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    cfg.ensure_dirs()
    native = load_variant("native")
    if native is None:
        sys.exit("Missing proximity_t<T>.csv for some thresholds - run scripts 05/07 first.")
    common = load_variant("common")
    variants = {"native": native}
    if common is not None:
        variants["common"] = common

    bbb_f = cfg.TABLES / "drug_bbb_evidence.csv"
    if bbb_f.exists():
        bbb = pd.read_csv(bbb_f).set_index("drug_id")[["bbb_label", "logBB_median"]]
    else:
        print("(drug_bbb_evidence.csv not found - run script 06; BBB columns will be 'no_evidence')")
        bbb = pd.DataFrame(columns=["bbb_label", "logBB_median"])

    # ------------------------------------------------------------ threshold stability
    say("=== STABILITY ACROSS STRING THRESHOLDS ===")
    stab = []
    for variant, tables in variants.items():
        drugs = common_drugs(tables)
        direct = direct_hit_drugs(tables)
        say(f"[{variant}] drugs analysable at all thresholds: {len(drugs)}; "
            f"of these {len(drugs & direct)} are direct hits")
        stab.append(stability(tables, drugs, variant, "all_drugs"))
        stab.append(stability(tables, drugs - direct, variant, "no_direct_hit"))
    stab = pd.concat(stab, ignore_index=True)
    stab.to_csv(cfg.TABLES / "threshold_stability.csv", index=False)
    say("\n" + stab.round(3).to_string(index=False))

    # ------------------------------------------------------------ significant drugs
    sig = []
    for variant, tables in variants.items():
        for t, df in tables.items():
            sig.append({"variant": variant, "threshold": t, "drugs_analysed": len(df),
                        "q_below_0.05": int((df["q_value"] < 0.05).sum()),
                        "p_below_0.05": int((df["p_single"] < 0.05).sum()),
                        "expected_by_chance": round(0.05 * len(df), 1),
                        "direct_hits": int(df["direct_hit"].sum()) if "direct_hit" in df else None})
    sig = pd.DataFrame(sig)
    sig.to_csv(cfg.TABLES / "significant_counts.csv", index=False)
    say("\n=== DRUGS CLOSER THAN CHANCE ===\n" + sig.to_string(index=False))

    # ------------------------------------------------------------ per-drug rank table (native)
    ts = cfg.THRESHOLDS
    drugs = sorted(common_drugs(native))
    direct = direct_hit_drugs(native)
    rk = pd.DataFrame(index=drugs)
    rk["drug_name"] = native[ts[0]].loc[drugs, "drug_name"]
    rk["n_targets_in_network"] = native[ts[-1]].loc[drugs, "n_targets_in_network"]
    rk["direct_hit"] = [d in direct for d in drugs]
    for t in ts:
        rk[f"z_{t}"] = native[t].loc[drugs, "z_score"]
        rk[f"q_{t}"] = native[t].loc[drugs, "q_value"]
        rk[f"rank_{t}"] = rk[f"z_{t}"].rank(method="min").astype(int)        # among all drugs
        # rank among INDIRECT drugs only (direct hits stay NaN)
        rk[f"rank_ind_{t}"] = rk.loc[~rk["direct_hit"], f"z_{t}"].rank(method="min")
    rk["mean_rank"] = rk[[f"rank_{t}" for t in ts]].mean(axis=1)
    ind_cols = [f"rank_ind_{t}" for t in ts]
    rk["mean_rank_indirect"] = rk[ind_cols].mean(axis=1)
    rk["n_thresholds_top50_indirect"] = (rk[ind_cols] <= 50).sum(axis=1)

    # CANDIDATE CLASSES. Direct hits are kept apart: their rank only says "the target is an
    # Alzheimer's gene". Among the other drugs, 'robust' = top 50 at ALL thresholds,
    # 'recurrent' = top 50 at 2 of 3.
    rk["candidate_class"] = "none"
    rk.loc[rk["direct_hit"], "candidate_class"] = "direct_hit"
    rk.loc[~rk["direct_hit"] & (rk["n_thresholds_top50_indirect"] >= 2), "candidate_class"] = "recurrent_indirect"
    rk.loc[~rk["direct_hit"] & (rk["n_thresholds_top50_indirect"] == 3), "candidate_class"] = "robust_indirect"
    rk = rk.join(bbb, how="left")
    rk["bbb_label"] = rk["bbb_label"].fillna("no_evidence")
    rk.index.name = "drug_id"
    rk.sort_values("mean_rank").to_csv(cfg.TABLES / "rank_stability_all_drugs.csv")
    ind = (rk[rk["candidate_class"].isin(["robust_indirect", "recurrent_indirect"])]
           .sort_values(["n_thresholds_top50_indirect", "mean_rank_indirect"], ascending=[False, True]))
    ind.to_csv(cfg.TABLES / "top_candidates_indirect.csv")
    dh = rk[rk["direct_hit"]].sort_values("mean_rank")
    dh.to_csv(cfg.TABLES / "direct_hit_drugs.csv")

    say(f"\n=== CANDIDATE CLASSES (native, {len(rk)} drugs) ===\n" + rk["candidate_class"].value_counts().to_string())
    say(f"\nDirect hits by BBB evidence: {dh['bbb_label'].value_counts().to_dict()}")
    say(f"Indirect candidates (robust + recurrent) by BBB evidence: {ind['bbb_label'].value_counts().to_dict()}")
    say("\nDirect-hit drugs, best 15 (target is itself an Alzheimer's gene):\n" +
        dh[["drug_name", "n_targets_in_network", "bbb_label", "mean_rank"]].head(15).to_string())
    say("\nIndirect candidates (top 50 among non-direct drugs at >=2 thresholds):\n" +
        ind[["drug_name", "n_targets_in_network", "n_thresholds_top50_indirect", "bbb_label",
             "mean_rank_indirect"]].head(25).to_string())

    # ------------------------------------------------------------ figures
    # Fig 1: rank-vs-rank scatter between thresholds
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, (t1, t2) in zip(axes, itertools.combinations(ts, 2)):
        rho = spearmanr(rk[f"z_{t1}"], rk[f"z_{t2}"])[0]
        ax.scatter(rk[f"rank_{t1}"], rk[f"rank_{t2}"], s=6, alpha=.4, c=["#d62828" if d else "#264653" for d in rk["direct_hit"]])
        ax.set(xlabel=f"rank at STRING ≥ {t1}", ylabel=f"rank at STRING ≥ {t2}", title=f"Spearman ρ = {rho:.2f}")
    fig.suptitle("Drug ranking across STRING thresholds (red = direct hit)", y=1.02)
    fig.tight_layout(); fig.savefig(cfg.FIGURES / "fig1_rank_scatter.png", dpi=200, bbox_inches="tight"); plt.close(fig)

    # Fig 2: top-N overlap, one panel per variant
    fig, axes = plt.subplots(1, len(variants), figsize=(5.5 * len(variants), 4), sharey=True, squeeze=False)
    for ax, variant in zip(axes[0], variants):
        sub = stab[(stab.variant == variant) & (stab.view == "all_drugs")]
        width = 0.25
        for i, (_, r) in enumerate(sub.iterrows()):
            ax.bar([j + i * width for j in range(len(TOP_NS))], [r[f"top{n}_overlap"] for n in TOP_NS],
                   width, label=f"{int(r.threshold_a)} vs {int(r.threshold_b)}")
        ax.set_xticks([j + width for j in range(len(TOP_NS))]); ax.set_xticklabels([f"top {n}" for n in TOP_NS])
        ax.set(title=f"{variant} Alzheimer's gene set", ylabel="drugs shared in the top-N"); ax.legend()
    fig.tight_layout(); fig.savefig(cfg.FIGURES / "fig2_topN_overlap.png", dpi=200); plt.close(fig)

    # Fig 3: BBB evidence among the top 50 at each threshold (all drugs vs direct hits excluded)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, (prefix, title) in zip(axes, [("rank_", "all drugs"), ("rank_ind_", "direct hits excluded")]):
        bottoms = pd.Series(0, index=[str(t) for t in ts])
        for label, color in BBB_COLORS.items():
            counts = pd.Series({str(t): int((rk.loc[rk[f"{prefix}{t}"] <= 50, "bbb_label"] == label).sum()) for t in ts})
            ax.bar(counts.index, counts.values, bottom=bottoms.values, color=color, label=label)
            bottoms = bottoms + counts
        ax.set(xlabel="STRING score threshold", title=f"Top 50: {title}")
    axes[0].set_ylabel("drugs in the top 50"); axes[1].legend(fontsize=8)
    fig.suptitle("BBB evidence among top-ranked drugs", y=1.02)
    fig.tight_layout(); fig.savefig(cfg.FIGURES / "fig3_bbb_top50.png", dpi=200, bbox_inches="tight"); plt.close(fig)

    # Fig 4: how much of the data survives each threshold
    qc_f = cfg.TABLES / "network_qc.csv"
    if qc_f.exists():
        qc = pd.read_csv(qc_f)
        fig, ax = plt.subplots(figsize=(6, 4))
        x = range(len(qc))
        ax.bar([i - .27 for i in x], qc["ad_coverage"], .27, label="Alzheimer's genes in network")
        ax.bar(list(x), qc["target_coverage"], .27, label="drug-target genes in network")
        ax.bar([i + .27 for i in x], qc["drugs_with_target_in_lcc"] / qc["drugs_total"], .27, label="drugs with ≥1 target in network")
        ax.set_xticks(list(x)); ax.set_xticklabels(qc["threshold"].astype(str))
        ax.set(xlabel="STRING score threshold", ylabel="fraction retained", ylim=(0, 1.05), title="Data coverage by network")
        ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(cfg.FIGURES / "fig4_network_coverage.png", dpi=200); plt.close(fig)

    (cfg.TABLES / "summary.txt").write_text("\n".join(lines))
    print(f"\nSaved tables to {cfg.TABLES} and figures to {cfg.FIGURES}")


if __name__ == "__main__":
    main()
