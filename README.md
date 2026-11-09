# Network-Based Drug Repurposing for Alzheimer's Disease
### Integrating protein-interaction proximity with blood–brain barrier (BBB) evidence

> **Status:** pipeline complete. All results are from actual outputs; see `results/RESULTS.md`.

## Research question
Can network proximity between drug targets and Alzheimer's-associated proteins, combined with
BBB permeability evidence, improve the prioritisation of drug-repurposing candidates?

## Background
Network medicine views disease as a perturbation of a local neighbourhood of the human
interactome (Barabási et al. 2011). Guney et al. (2016) showed that drugs whose targets lie
close, in network distance, to a disease's genes are more likely to be therapeutically relevant;
closeness is judged against random target sets of matched degree. For brain diseases a candidate
must also reach the brain, so we add blood–brain barrier (BBB) evidence to the ranking.
Methodological reference: **NetMedPy** (Aldana et al., *Bioinformatics* 2025). We use its distance,
null-model and screening functions, and write our own data processing, BBB integration,
threshold sweep and evaluation code. Our dataset (Alzheimer's disease, approved drugs) and both
extensions below differ from the tutorial's Vitamin D example.

## What is new relative to the tutorial
1. **BBB integration** – drugs are re-ranked with experimental BBB labels / logBB from B3DB.
   Drugs with no B3DB record are labelled *"no evidence"*, never *"cannot cross"*.
2. **STRING confidence sweep** – the whole analysis is repeated at combined-score thresholds
   400 / 700 / 900 and ranking stability is quantified (Spearman correlation, top-N overlap).

## Data sources
| Data | Source | Used for |
|---|---|---|
| Human physical PPI | STRING v12.0 | interaction network |
| Alzheimer's-associated genes | Open Targets Platform | disease gene set |
| Drug → target | ChEMBL (approved small molecules, mechanism-of-action targets) | drug targets |
| BBB labels / logBB | B3DB (Meng et al., *Sci Data* 2021) | BBB evidence |

**References** (verify details before submission)
1. Barabási A-L, Gulbahce N, Loscalzo J. Network medicine. *Nat Rev Genet* 2011;12:56–68.
2. Guney E, et al. Network-based in silico drug efficacy screening. *Nat Commun* 2016;7:10331.
3. Aldana A, et al. NetMedPy: a Python package for large-scale network medicine screening. *Bioinformatics* 2025.
4. Szklarczyk D, et al. The STRING database in 2023. *Nucleic Acids Res* 2023;51:D638–D646.
5. Ochoa D, et al. The next-generation Open Targets Platform. *Nucleic Acids Res* 2023;51:D1353–D1359.
6. Zdrazil B, et al. The ChEMBL Database in 2023. *Nucleic Acids Res* 2024;52:D1180–D1192.
7. Meng F, et al. A curated diverse molecular database of blood-brain barrier permeability with chemical descriptors. *Sci Data* 2021;8:289.

## Installation
Use Python **3.11** (recommended: older NetMedPy releases documented conflicts above 3.11).
```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running the pipeline
```bash
python scripts/01_download_data.py          # downloads raw data into data/raw/
python scripts/02_prepare_ad_genes.py
python scripts/03_prepare_drug_targets.py
python scripts/04_build_ppi_network.py      # builds networks + coverage QC table
python scripts/05_calculate_proximity.py --threshold 700
```
> Scripts must be run from a file (not interactively) because NetMedPy uses multiprocessing
> with the `spawn` start method on macOS.

## Methodology
Full details in `docs/methodology.md`. Key choices: Open Targets `known_drug` evidence is
excluded when defining Alzheimer's genes (circularity); drug targets are ChEMBL mechanism-of-action
targets; the primary Alzheimer's gene set is genetic evidence only (Open Targets genetic_association + genetic_literature, score >= 0.3; chosen because literature/text-mining evidence inflated the set to >3000 genes and made over half of all drug targets disease genes); distances are shortest paths; null model is degree-matched, 1000 draws.

## Results
Full tables, numbers and caveats are in [`results/RESULTS.md`](results/RESULTS.md). Main findings:

* The 573 Alzheimer's genetic-evidence genes form a **significant connected module** in STRING
  physical networks at thresholds 400, 700 and 900.
* **No overall enrichment**: only nirogacestat (a gamma-secretase inhibitor hitting PSEN1/PSEN2) has
  q < 0.05, at all three thresholds, and the number of drugs with p < 0.05 is at or below chance
  expectation.
* The top of the ranking is filled by **direct hits** (a target that is itself an Alzheimer's
  gene). These are reported separately.
* **Extra experiment 1 (STRING threshold sweep):** rankings of non-direct-hit drugs are unstable
  (Spearman ρ 0.26–0.52; 3–26 of the top 50 shared between networks). The same holds when all
  networks use an identical gene set, so it is not just gene dropout.
* **Extra experiment 2 (BBB evidence):** 35% of drugs have no B3DB record; excluding BBB-negative
  drugs removes 18–22 of the top 50, and requiring BBB+ evidence keeps only 4–5.

Figures: `results/figures/` (rank scatter, top-N overlap, BBB evidence, data coverage).

## Limitations
* Results depend on the gene-set definition (genetic evidence, score >= 0.3); a broader set made
  proximity uninformative.
* ChEMBL mechanism-of-action targets only (median one target per drug), so distances are coarse and
  off-target effects are ignored.
* Direct hits are not informative about proximity, and drugs sharing a target are not independent.
* Coverage and BBB data are incomplete (44% of Alzheimer's genes at threshold 900; 35% of drugs lack
  B3DB records, and missing is not the same as impermeable).
* Proximity shows closeness, not direction of effect; BBB+ is a permeability label, not proof of
  therapeutic brain exposure.
* The noise floor of the random sampling and other interactomes were not tested; no experimental
  validation was done.

See `results/RESULTS.md` section 7 for details and future work.
