# Network-Based Drug Repurposing for Alzheimer's Disease
### Integrating protein-interaction proximity with blood–brain barrier (BBB) evidence

> **Status:** work in progress. Sections marked `TODO` are completed from real outputs only.

## Research question
Can network proximity between drug targets and Alzheimer's-associated proteins, combined with
BBB permeability evidence, improve the prioritisation of drug-repurposing candidates?

## Background
TODO (network medicine, proximity measure of Guney et al. 2016, why the BBB matters for CNS drugs).
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

TODO: full citations.

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
# 06_bbb_integration.py, 07_threshold_sweep.py, 08_compare_results.py  (TODO)
```
> Scripts must be run from a file (not interactively) because NetMedPy uses multiprocessing
> with the `spawn` start method on macOS.

## Methodology
TODO – see `docs/methodology.md`. Key choices so far: Open Targets `known_drug` evidence is
excluded when defining Alzheimer's genes (circularity); drug targets are ChEMBL mechanism-of-action
targets; distances are shortest paths; null model is degree-matched, 1000 draws.

## Results
TODO – see `results/RESULTS.md` (filled from actual outputs only).

## Limitations
TODO.
