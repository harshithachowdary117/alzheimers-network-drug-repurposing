"""
03_prepare_drug_targets.py - build the approved-drug -> protein-target table.

INPUTS : data/raw/chembl_approved_molecules.json
         data/raw/chembl_mechanisms.json
         data/raw/chembl_targets.json
         data/raw/string_protein_aliases.txt.gz   (UniProt accession -> STRING gene name)
OUTPUTS: data/processed/drugs.csv           one row per approved small-molecule drug
         data/processed/drug_targets.csv    one row per (drug, target gene)

Design decisions (document these in the README)
  * Only approved (max_phase 4) SMALL MOLECULES - BBB data (B3DB) covers small molecules.
  * Targets = ChEMBL mechanism-of-action targets, i.e. the curated PRIMARY targets,
    not every weak bioassay hit. This is conservative; many drugs end up with 1-3 targets.
  * Human SINGLE PROTEIN / PROTEIN COMPLEX targets only.
  * Targets are converted to gene names through the UniProt accession -> STRING alias file,
    so the gene names match the PPI network exactly. If an accession is missing there we
    fall back to ChEMBL's GENE_SYMBOL synonym and record that in `mapping_method`.
  * Compounds are identified by InChIKey (stored) so script 06 can join to B3DB by
    chemical structure and not by drug name.
"""
import gzip
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

KEEP_TARGET_TYPES = {"SINGLE PROTEIN", "PROTEIN COMPLEX"}


def load_uniprot_to_gene() -> dict:
    """UniProt accession -> STRING preferred gene name, via the STRING alias file."""
    aliases = pd.read_csv(cfg.RAW / "string_protein_aliases.txt.gz", sep="\t",
                          compression="gzip", dtype=str)
    aliases.columns = ["protein_id", "alias", "source"]
    info = pd.read_csv(cfg.RAW / "string_protein_info.txt.gz", sep="\t", dtype=str)
    id2name = dict(zip(info.iloc[:, 0], info["preferred_name"]))
    uni = aliases[aliases["source"].str.contains("UniProt_AC", na=False)].copy()
    uni["gene"] = uni["protein_id"].map(id2name)
    uni = uni.dropna(subset=["gene"]).drop_duplicates("alias")
    return {a: cfg.norm_symbol(g) for a, g in zip(uni["alias"], uni["gene"])}


def target_genes(target: dict, uni2gene: dict):
    """Return [(gene, method), ...] for one ChEMBL target record."""
    if target.get("organism") != "Homo sapiens" or target.get("target_type") not in KEEP_TARGET_TYPES:
        return []
    out = []
    for comp in target.get("target_components") or []:
        acc = comp.get("accession")
        if acc in uni2gene:
            out.append((uni2gene[acc], "uniprot_alias"))
            continue
        # fallback: ChEMBL's own gene-symbol synonyms (may include old aliases)
        syms = [s["component_synonym"] for s in comp.get("target_component_synonyms") or []
                if s.get("syn_type") == "GENE_SYMBOL"]
        out.extend((cfg.norm_symbol(s), "chembl_symbol") for s in syms[:1])
    return out


def main():
    cfg.ensure_dirs()
    mols = pd.DataFrame(json.loads((cfg.RAW / "chembl_approved_molecules.json").read_text()))
    mechs = json.loads((cfg.RAW / "chembl_mechanisms.json").read_text())
    targets = {t["target_chembl_id"]: t for t in
               json.loads((cfg.RAW / "chembl_targets.json").read_text())}
    uni2gene = load_uniprot_to_gene()
    print(f"{len(uni2gene)} UniProt accessions mapped to STRING gene names")

    # ---- drugs: approved small molecules; group salts/forms under their parent compound
    mols = mols[mols["molecule_type"] == "Small molecule"].copy()
    mols["drug_id"] = mols["parent_chembl_id"].fillna(mols["molecule_chembl_id"])
    mols["inchikey14"] = mols["standard_inchi_key"].str[:14]    # 14 = connectivity layer
    print(f"{len(mols)} approved small-molecule records, {mols.drug_id.nunique()} parent compounds")
    mol2drug = dict(zip(mols["molecule_chembl_id"], mols["drug_id"]))

    # ---- mechanism records -> (drug, gene) pairs
    pairs, dropped = [], {"not_approved_small_molecule": 0, "no_human_protein_gene": 0}
    for m in mechs:
        drug = mol2drug.get(m.get("molecule_chembl_id"))
        if drug is None:
            dropped["not_approved_small_molecule"] += 1
            continue
        genes = target_genes(targets.get(m.get("target_chembl_id"), {}), uni2gene)
        if not genes:
            dropped["no_human_protein_gene"] += 1
        for gene, method in genes:
            pairs.append({"drug_id": drug, "gene_symbol": gene,
                          "target_chembl_id": m.get("target_chembl_id"),
                          "action_type": m.get("action_type"),
                          "mechanism_of_action": m.get("mechanism_of_action"),
                          "mapping_method": method})
    dt = pd.DataFrame(pairs).drop_duplicates(["drug_id", "gene_symbol"])
    print("Mechanism records dropped:", dropped)

    # ---- drug table restricted to drugs with at least one usable target
    drugs = (mols[mols["drug_id"].isin(dt["drug_id"])]
             .sort_values("pref_name", na_position="last")
             .drop_duplicates("drug_id")
             [["drug_id", "pref_name", "standard_inchi_key", "inchikey14"]]
             .rename(columns={"pref_name": "drug_name"}))
    dt = dt.merge(drugs[["drug_id", "drug_name"]], on="drug_id")
    drugs.to_csv(cfg.PROCESSED / "drugs.csv", index=False)
    dt.to_csv(cfg.PROCESSED / "drug_targets.csv", index=False)

    per_drug = dt.groupby("drug_id").size()
    print(f"\n{len(drugs)} drugs with >=1 target; {dt.gene_symbol.nunique()} distinct target genes")
    print("Targets per drug:\n", per_drug.describe().round(2).to_string())
    print("Mapping methods:\n", dt["mapping_method"].value_counts().to_string())
    print("Drugs missing an InChIKey:", int(drugs["standard_inchi_key"].isna().sum()))


if __name__ == "__main__":
    main()
