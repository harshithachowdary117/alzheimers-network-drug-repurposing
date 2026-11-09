"""
config.py - single place for paths, URLs and analysis parameters.

Every other script does `import config as cfg`, so a parameter changed here
is changed consistently across the whole pipeline (important for reproducibility).
"""
from pathlib import Path

# ---------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"


def ensure_dirs():
    """Create the output folders if they do not exist yet."""
    for p in (RAW, PROCESSED, TABLES, FIGURES):
        p.mkdir(parents=True, exist_ok=True)


def norm_symbol(s):
    """Gene symbols are compared in UPPER case so 'Apoe' and 'APOE' match."""
    return str(s).strip().upper()


# ---------------------------------------------------------------- data sources
# STRING v12.0 *physical* interaction network (the same file the NetMedPy
# tutorial uses). Scores run 0-1000 ("combined_score").
STRING_LINKS_URL = ("https://stringdb-downloads.org/download/"
                    "protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz")
STRING_INFO_URL = ("https://stringdb-downloads.org/download/"
                   "protein.info.v12.0/9606.protein.info.v12.0.txt.gz")
STRING_ALIASES_URL = ("https://stringdb-downloads.org/download/"
                      "protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz")

# B3DB: blood-brain barrier labels (BBB+/BBB-) and logBB values (CC0 licence)
B3DB_BASE = "https://raw.githubusercontent.com/theochem/B3DB/main/B3DB"
B3DB_CLASS_URL = f"{B3DB_BASE}/B3DB_classification.tsv"
B3DB_REG_URL = f"{B3DB_BASE}/B3DB_regression.tsv"

# Open Targets Platform GraphQL API (Alzheimer's disease-target associations)
OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
# Open Targets has used both IDs for Alzheimer's disease; the script tries them in order
AD_DISEASE_IDS = ["MONDO_0004975", "EFO_0000249"]

# ChEMBL REST API (approved drugs + mechanism-of-action targets)
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

# ---------------------------------------------------------------- analysis parameters
# --- Alzheimer's gene set (script 02)
# Open Targets "known_drug" evidence is built from drugs already tested in the
# disease. Using it to define the disease genes and then asking which drugs hit
# those genes is circular, so that datatype is excluded from the gene-selection score.
EXCLUDE_DATATYPES = {"known_drug", "clinical"}
AD_MIN_SCORE = 0.3      # keep genes whose best non-excluded datatype score >= this

# --- PPI thresholds on STRING combined_score (0-1000). STRING's own labels:
# 400 = medium confidence, 700 = high, 900 = highest.
THRESHOLDS = [400, 700, 900]

# --- proximity settings (script 05)
DISTANCE = "shortest_path"   # cheapest metric; random-walk matrices are much heavier
NULL_MODEL = "log_binning"  # NetMedPy options: degree_match, log_binning, uniform
N_ITER = 1000                # random draws per drug
MIN_TARGETS = 1              # minimum drug targets present in the network

# --- PRIMARY Alzheimer's gene set: genetic evidence only (decided after inspecting
# per-datatype gene counts). 'literature' (text-mining) and 'affected_pathway' inflate the
# set to >3000 genes, and over half of all drug targets then become disease genes.
AD_DATATYPES = {"genetic_association", "genetic_literature"}
