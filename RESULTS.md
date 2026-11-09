# Results

All numbers below come from the pipeline outputs in `results/tables/` (STRING v12.0 physical
network, 573 Alzheimer's genes, 1,010 approved small-molecule drugs, log-binning null model,
1,000 random draws per drug). One item marked **[VERIFY]** (a literature citation) must be checked before
submission.

## Summary

1. **The Alzheimer's genes form a significant module** in the PPI network at all three STRING
   thresholds, so network proximity is a meaningful question for this gene set.
2. **Almost no drug is closer to the Alzheimer's genes than chance.** Only one drug,
   nirogacestat, has q < 0.05, and it does so at all three thresholds. The number of drugs with p < 0.05 (33–44 per network) is at or below what chance
   alone would give (36–50).
3. **The top of every ranking is dominated by "direct hits"**: drugs with a target that is itself
   an Alzheimer's gene. Their proximity is a restatement of gene-set membership, and the 63 such
   drugs hit only 20 distinct genes (ACE, DPP4, EGFR, RET, PDE3A and CSF1R account for most).
4. **Rankings of the remaining (indirect) drugs are not stable across STRING thresholds**:
   Spearman ρ = 0.26–0.52, and only 3–26 of the top 50 drugs are shared between networks.
5. **Requiring BBB evidence removes most top-ranked drugs**, and about a third of all drugs have
   no BBB record at all.

## 1. Data and network coverage

| STRING threshold | LCC nodes | LCC edges | Alzheimer's genes in network | Drug-target genes in network | Drugs with ≥1 target in network |
|---|---|---|---|---|---|
| ≥ 400 | 15,640 | 210,570 | 483 / 573 (84.3%) | 443 / 452 (98.0%) | 993 / 1,010 |
| ≥ 700 | 9,830 | 85,576 | 328 / 573 (57.2%) | 378 / 452 (83.6%) | 897 / 1,010 |
| ≥ 900 | 7,167 | 43,180 | 252 / 573 (44.0%) | 313 / 452 (69.2%) | 722 / 1,010 |

Raising the threshold removes genes and drugs, and poorly studied genes disappear first.
Only 252 Alzheimer's genes and 722 drugs exist at every threshold.

**Gene-set choice.** Using all Open Targets evidence gave 3,536 genes, and 238 of the 452
drug-target genes (53%) were themselves Alzheimer's genes. Restricting to genetic evidence
(score ≥ 0.3) gave 573 genes with 20 of 452 overlapping, leaving room for proximity to
discriminate between drugs.

**Disease module.** The Alzheimer's genes in each network form a larger connected component than
expected by chance (log-binning null, 1,000 draws):

| Threshold | Largest connected component of the genes | z-score | p |
|---|---|---|---|
| 400 | 174 | 2.47 | < 0.001 |
| 700 | 53 | 5.06 | 0.002 |
| 900 | 17 | 4.02 | 0.010 |

At 900 the module is small (17 genes), so it is the weakest network for this analysis.

## 2. Proximity results

| Variant | Threshold | Drugs analysed | q < 0.05 | p < 0.05 | Expected by chance | Direct hits |
|---|---|---|---|---|---|---|
| native | 400 | 993 | 1 | 39 | 49.7 | 63 |
| native | 700 | 897 | 1 | 44 | 44.9 | 62 |
| native | 900 | 722 | 1 | 33 | 36.1 | 56 |
| common | 400 | 993 | 1 | 43 | 49.7 | 56 |
| common | 700 | 897 | 1 | 42 | 44.9 | 56 |
| common | 900 | 722 | 1 | 32 | 36.1 | 56 |

"Native" uses the Alzheimer's genes present in each network; "common" uses the 252 genes present
in all three.

* The same single drug passes q < 0.05 at all three thresholds: **nirogacestat** (6 targets in the
  network, mean distance 0.5, p = 0 and q = 0 in each: none of the 1,000 random draws was as close).
  Nirogacestat is a gamma-secretase inhibitor, and its targets include PSEN1 and PSEN2, so this is
  the one clearly biologically expected result. It is also the only significant result that
  survives changing the STRING threshold.
* The number of nominally significant drugs is not higher than chance would produce, so there is
  no evidence that approved drugs as a group lie closer to Alzheimer's genes than random target
  sets. With discrete path lengths the empirical p-values are conservative, so this does not prove
  that no real signal exists.

## 3. Direct hits

56 of the 722 drugs present at every threshold have a target that is itself an Alzheimer's gene.
At threshold 700, 42 drugs have distance 0. Their empirical p-values lie in a narrow band
(0.024-0.097, median 0.036) and their q-values are 0.94-1.0, so none is significant. This is a
floor set by the null: any random target that is itself an Alzheimer's gene also gives distance 0,
and 328 of the 9,830 nodes in this network (3.3%) are Alzheimer's genes, similar to the median
p-value. Their z-scores are nevertheless large and negative (median -2.5, minimum -6.0), most likely
because the null distribution is discrete and dominated by distances above 0, so a normal z-score
overstates the evidence. For direct hits the p-value is the safer statistic, and their ranks are not
comparable with those of other drugs.

BBB evidence for the 56 direct hits: 5 BBB+, 19 BBB−, 4 conflicting, 28 no evidence.
Why these targets count as Alzheimer's genes (checked in `ad_genes_all.csv`; selection score =
best of `genetic_association` and `genetic_literature`):

| Gene | genetic_association | genetic_literature | Selection score |
|---|---|---|---|
| ACE | 0.888 | 0.608 | 0.888 |
| EGFR | 0.696 | 0.000 | 0.696 |
| HCN4 | 0.580 | 0.000 | 0.580 |
| DPP4 | 0.405 | 0.000 | 0.405 |

All four pass the 0.3 cut-off through Open Targets genetic-association scores. Such scores can come
from locus-to-gene assignment of genome-wide association signals, in which a signal is attributed to
a gene without the gene being shown to be causal. A drug hitting one of these genes is therefore a
direct hit only as far as that assignment is correct. FLT3 is **not** in the gene set (both scores 0); quizartinib is a direct hit through CSF1R and RET.

**Which genes explain the direct hits.** Before network filtering, 63 drugs have a target that is an
Alzheimer's gene, and these drugs hit only **20 distinct genes**. Six genes account for most of them:

| Alzheimer's gene | Drugs hitting it | Examples |
|---|---|---|
| ACE | 14 | benazepril, captopril, cilazapril, enalapril |
| DPP4 | 11 | alogliptin, anagliptin, evogliptin, gosogliptin |
| EGFR | 11 | afatinib, brigatinib, dacomitinib, erlotinib |
| RET | 8 | alectinib, pralsetinib, quizartinib, regorafenib |
| PDE3A | 5 | cilostazol, dyphylline, ibudilast, inamrinone |
| CSF1R | 4 | pazopanib, pexidartinib, quizartinib, sunitinib |

The others are single-drug or two-drug genes: HCN4 (dronedarone, ivabradine), MME (racecadotril,
sacubitril), and one each for ABCA1 (probucol), SLC6A1 (tiagabine), KCNJ3 (vernakalant), PSEN1 and
PSEN2 (nirogacestat), and NDUFAB1 and MT-ND1 (metformin).

So the 56 direct hits in the stability analysis are mostly a few drug classes (ACE inhibitors,
gliptins, EGFR inhibitors, RET inhibitors), not independent signals. Nirogacestat hits the
gamma-secretase subunits PSEN1 and PSEN2, the best-known Alzheimer's genes. MME encodes neprilysin,
an enzyme that degrades amyloid-beta, so inhibiting it (as sacubitril does) is not obviously
beneficial; this illustrates that proximity does not indicate benefit.

## 4. Additional experiment: stability across STRING thresholds

Only the 722 drugs analysable at all thresholds are compared.

**All drugs (native)**

| Pair | Spearman ρ | Top-10 shared | Top-20 shared | Top-50 shared |
|---|---|---|---|---|
| 400 vs 700 | 0.508 | 7 | 10 | 47 |
| 400 vs 900 | 0.416 | 7 | 14 | 47 |
| 700 vs 900 | 0.623 | 8 | 12 | 49 |

**Direct hits excluded (666 drugs, native)**

| Pair | Spearman ρ | Top-10 shared | Top-20 shared | Top-50 shared |
|---|---|---|---|---|
| 400 vs 700 | 0.378 | 1 | 1 | 7 |
| 400 vs 900 | 0.261 | 1 | 3 | 3 |
| 700 vs 900 | 0.524 | 2 | 6 | 26 |

The all-drug rankings look moderately stable only because the same 56 direct hits fill the top
at every threshold. Without them, the top-10 lists share one or two drugs, and the 400 and 900
networks share only 3 of their top 50.

**Is it the genes or the edges?** In the common-gene variant all networks use the same 252
Alzheimer's genes. Correlations without direct hits are similar or lower (ρ = 0.31, 0.19, 0.50), so
the instability is not explained by Alzheimer's genes dropping out at higher thresholds. It is
associated with changing the network edges.

Only 30 non-direct-hit drugs are in the top 50 at two or more thresholds, and only 3 at all three:
**pentostatin** (BBB+), **doxycycline** and **lomitapide** (both no BBB evidence).
Others recurring in two networks include omaveloxolone, nintedanib, elafibranor, alpidem (BBB+),
zaleplon (BBB+), several progestins, PI3K inhibitors (idelalisib, copanlisib, duvelisib,
leniolisib) and CDK4/6 inhibitors (abemaciclib, ribociclib, trilaciclib). Drugs of one class share
a target, so these also are not independent findings. See `top_candidates_indirect.csv`.

## 5. BBB integration

B3DB gave labels for 653 of the 1,010 drugs (BBB+ 315, BBB− 271, conflicting 67); **357 (35%) have
no record**. Absence from B3DB is treated as missing evidence, not as impermeability.

Effect on the top of the proximity ranking (strategy A = proximity only; B = BBB-negative drugs
removed; C = BBB+ evidence required):

| Threshold | Top-N | B shares with A | C shares with A | In A: BBB+ | BBB− | no evidence |
|---|---|---|---|---|---|---|
| 400 | 10 | 8 | 2 | 2 | 2 | 6 |
| 400 | 20 | 11 | 3 | 3 | 9 | 8 |
| 400 | 50 | 28 | 5 | 5 | 22 | 19 |
| 700 | 10 | 9 | 3 | 3 | 1 | 6 |
| 700 | 20 | 14 | 3 | 3 | 6 | 10 |
| 700 | 50 | 29 | 5 | 5 | 21 | 20 |
| 900 | 10 | 10 | 3 | 3 | 0 | 7 |
| 900 | 20 | 13 | 3 | 3 | 7 | 9 |
| 900 | 50 | 32 | 4 | 4 | 18 | 24 |

* Removing BBB-negative drugs (B) displaces 18–22 of the top 50, so BBB evidence changes the
  ranking substantially.
* Requiring positive BBB evidence (C) keeps only 4–5 of the top 50, mainly because 19–24 of those 50
  have no B3DB record.
* Among the 30 recurring indirect candidates, BBB evidence is: 4 BBB+, 4 conflicting, 2 BBB−,
  20 none.

## 6. What can and cannot be concluded

**Can be concluded**
* Alzheimer's genetic-evidence genes form a significant module in STRING physical networks.
* In this setting network proximity gives no overall enrichment of approved drugs near these genes.
* Rankings of non-direct-hit drugs depend strongly on the STRING threshold, so a drug should not be
  presented as a network-based candidate on the basis of a single network.
* BBB evidence is missing for about a third of approved drugs and, where present, removes
  many top-ranked drugs.

**Cannot be concluded**
* That any listed drug would treat Alzheimer's. Proximity measures closeness, not direction of
  effect: a drug near disease genes could worsen the disease (the gamma-secretase inhibitor
  semagacestat worsened cognition in a phase 3 trial; Doody et al., *NEJM* 2013 **[VERIFY citation]**).
* That a BBB+ label means therapeutic brain exposure. B3DB labels permeability, not dose-dependent
  brain levels, and efflux transport matters.
* That the instability is purely due to the network. Random-sampling noise within a threshold was
  not separately measured (see limitations).

## 7. Limitations and future work

* **Gene set.** Results depend on the choice of genetic evidence and a 0.3 cut-off. A broader gene
  set made proximity uninformative; a stricter one would reduce coverage further.
* **Drug targets.** Only ChEMBL mechanism-of-action targets were used (median one per drug).
  Single-target drugs give coarse, nearly discrete distances. Off-target effects are ignored.
* **Discrete distances and ties.** Shortest-path lengths take few values, which makes p-values
  conservative and limits power. Continuous distances (random-walk) were not tested.
* **Direct hits.** Their ranks are not informative about network proximity, and at threshold 700 their
  p-values were never below 0.024.
* **Correlated drugs.** Drugs sharing a target appear as separate entries, overstating how many
  independent candidates there are.
* **Missing data.** Coverage falls to 44% of Alzheimer's genes and 722 of 1,010 drugs at threshold 900;
  poorly studied genes drop out first, which favours well-studied biology. B3DB lacks
  records for 35% of drugs, particularly recent ones.
* **Noise floor.** The ranking was not repeated at one threshold with different random seeds, so
  the part of the instability due to sampling noise is unquantified.
* **Single interactome.** Only STRING physical interactions were used. Other networks (for example
  BioGRID, as in the NetMedPy robustness example) were not tested.
* **No experimental validation.** All findings are computational.
