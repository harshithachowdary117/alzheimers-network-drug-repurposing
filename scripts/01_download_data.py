"""
01_download_data.py - download every raw dataset the project needs.

INPUTS : none (needs an internet connection)
OUTPUTS: files in data/raw/ plus data/raw/download_manifest.json
         (URL, size, SHA-256 and timestamp for each file, for reproducibility)

Sources
  STRING v12.0   human physical PPI links, protein info, protein aliases
  B3DB           blood-brain-barrier classification + logBB tables (GitHub)
  Open Targets   Alzheimer's disease -> target associations (GraphQL API)
  ChEMBL         approved small-molecule drugs, their mechanism-of-action
                 targets and the target annotations (REST API)

Re-running is safe: files that already exist are skipped (use --force to redo).

Usage
  python scripts/01_download_data.py
  python scripts/01_download_data.py --only string b3db
"""
import argparse
import datetime as dt
import hashlib
import json
import sys
import time
from pathlib import Path

import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

CHUNK_TARGET_IDS = 50     # ChEMBL target ids per request (keeps URLs short)
PAGE = 1000               # ChEMBL page size (API maximum)


# ------------------------------------------------------------------ helpers
def make_session():
    """A requests session with a polite User-Agent."""
    s = requests.Session()
    s.headers.update({"User-Agent": "AD-network-repurposing-student-project/1.0"})
    return s


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download_file(session, url: str, dest: Path, force: bool = False):
    """Stream a (possibly large) file to disk with a progress bar.

    The file is first written to '<name>.part' and renamed only when complete,
    so an interrupted download never leaves a half-written file behind.
    """
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"  [skip] {dest.name} already exists")
        return
    tmp = dest.with_name(dest.name + ".part")
    with session.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(tmp, "wb") as f, tqdm(total=total, unit="B", unit_scale=True,
                                        desc=dest.name) as bar:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
                bar.update(len(chunk))
    tmp.replace(dest)


def request_json(session, method, url, retries=4, **kwargs):
    """GET/POST returning parsed JSON, with exponential back-off on failures."""
    for attempt in range(retries):
        try:
            r = session.request(method, url, timeout=120, **kwargs)
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, ValueError) as err:
            wait = 2 ** attempt
            print(f"    request failed ({err}); retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"Giving up on {url}")


# ------------------------------------------------------------------ STRING
def get_string(session, force):
    print("\n== STRING v12.0 ==")
    download_file(session, cfg.STRING_LINKS_URL, cfg.RAW / "string_physical_links.txt.gz", force)
    download_file(session, cfg.STRING_INFO_URL, cfg.RAW / "string_protein_info.txt.gz", force)
    download_file(session, cfg.STRING_ALIASES_URL, cfg.RAW / "string_protein_aliases.txt.gz", force)


# ------------------------------------------------------------------ B3DB
def get_b3db(session, force):
    print("\n== B3DB ==")
    download_file(session, cfg.B3DB_CLASS_URL, cfg.RAW / "B3DB_classification.tsv", force)
    download_file(session, cfg.B3DB_REG_URL, cfg.RAW / "B3DB_regression.tsv", force)


# ------------------------------------------------------------------ Open Targets
OT_QUERY = """
query AssociatedTargets($id: String!, $index: Int!, $size: Int!) {
  disease(efoId: $id) {
    id
    name
    associatedTargets(page: {index: $index, size: $size}) {
      count
      rows {
        target { id approvedSymbol }
        score
        datatypeScores { id score }
      }
    }
  }
}
"""


def get_open_targets(session, force):
    """Download all targets associated with Alzheimer's disease, with per-datatype scores.

    Datatype scores (genetic_association, literature, known_drug, ...) are kept
    because script 02 needs them to exclude the circular 'known_drug' evidence.
    """
    print("\n== Open Targets (Alzheimer's disease) ==")
    out = cfg.RAW / "opentargets_alzheimer_associations.json"
    if out.exists() and not force:
        print(f"  [skip] {out.name} already exists")
        return
    size = 500
    for disease_id in cfg.AD_DISEASE_IDS:
        rows, index, name, count = [], 0, None, None
        while True:
            payload = {"query": OT_QUERY,
                       "variables": {"id": disease_id, "index": index, "size": size}}
            resp = request_json(session, "POST", cfg.OT_GRAPHQL, json=payload)
            disease = (resp.get("data") or {}).get("disease")
            if disease is None:      # this ID is unknown to the current release
                break
            name = disease["name"]
            block = disease["associatedTargets"]
            count = block["count"]
            rows.extend(block["rows"])
            print(f"  {disease_id}: {len(rows)}/{count} targets", end="\r")
            if len(rows) >= count or not block["rows"]:
                break
            index += 1
        if rows:
            print()
            record = {"disease_id": disease_id, "disease_name": name, "count": count,
                      "retrieved": dt.date.today().isoformat(), "rows": rows}
            out.write_text(json.dumps(record))
            print(f"  saved {len(rows)} associations for '{name}' ({disease_id})")
            return
    raise RuntimeError("Open Targets returned no data for any Alzheimer's disease ID in "
                       "config.AD_DISEASE_IDS - check the API / IDs.")


# ------------------------------------------------------------------ ChEMBL
def chembl_paginate(session, endpoint, list_key, params):
    """Yield every record of a ChEMBL endpoint, following offset pagination."""
    url = f"{cfg.CHEMBL_API}/{endpoint}.json"
    offset, params = 0, dict(params, limit=PAGE)
    while True:
        data = request_json(session, "GET", url, params=dict(params, offset=offset))
        records = data[list_key]
        yield from records
        if not data["page_meta"].get("next"):
            break
        offset += PAGE


def get_chembl(session, force):
    print("\n== ChEMBL ==")
    mol_f = cfg.RAW / "chembl_approved_molecules.json"
    mech_f = cfg.RAW / "chembl_mechanisms.json"
    targ_f = cfg.RAW / "chembl_targets.json"
    if all(p.exists() for p in (mol_f, mech_f, targ_f)) and not force:
        print("  [skip] ChEMBL files already exist")
        return

    # 1) approved drugs (max_phase = 4). Keep only the fields we need.
    print("  fetching approved molecules ...")
    molecules = []
    for m in chembl_paginate(session, "molecule", "molecules", {"max_phase": 4}):
        molecules.append({
            "molecule_chembl_id": m.get("molecule_chembl_id"),
            "pref_name": m.get("pref_name"),
            "molecule_type": m.get("molecule_type"),
            "max_phase": m.get("max_phase"),
            "standard_inchi_key": (m.get("molecule_structures") or {}).get("standard_inchi_key"),
            "parent_chembl_id": (m.get("molecule_hierarchy") or {}).get("parent_chembl_id"),
        })
    mol_f.write_text(json.dumps(molecules))
    print(f"  saved {len(molecules)} approved molecules")

    # 2) mechanism-of-action records (drug -> target). Filtered to approved drugs in script 03.
    print("  fetching mechanisms of action ...")
    mechanisms = list(chembl_paginate(session, "mechanism", "mechanisms", {}))
    mech_f.write_text(json.dumps(mechanisms))
    print(f"  saved {len(mechanisms)} mechanism records")

    # 3) annotations (organism, type, UniProt accession, gene symbols) of the targets used
    target_ids = sorted({m["target_chembl_id"] for m in mechanisms if m.get("target_chembl_id")})
    print(f"  fetching {len(target_ids)} target annotations ...")
    targets = []
    for i in tqdm(range(0, len(target_ids), CHUNK_TARGET_IDS), desc="targets"):
        chunk = target_ids[i:i + CHUNK_TARGET_IDS]
        data = request_json(session, "GET", f"{cfg.CHEMBL_API}/target.json",
                            params={"target_chembl_id__in": ",".join(chunk), "limit": PAGE})
        targets.extend(data["targets"])
    targ_f.write_text(json.dumps(targets))
    print(f"  saved {len(targets)} targets")


# ------------------------------------------------------------------ manifest + main
def write_manifest():
    """Record size and checksum of every downloaded file (reproducibility record)."""
    files = sorted(p for p in cfg.RAW.iterdir()
                   if p.is_file() and p.name not in (".gitkeep", "download_manifest.json"))
    manifest = {"created": dt.datetime.now().isoformat(timespec="seconds"), "files": {}}
    print("\n== Verification ==")
    for p in files:
        size = p.stat().st_size
        manifest["files"][p.name] = {"bytes": size, "sha256": sha256(p)}
        print(f"  {p.name:<45} {size / 1e6:9.2f} MB")
    (cfg.RAW / "download_manifest.json").write_text(json.dumps(manifest, indent=2))


SOURCES = {"string": get_string, "b3db": get_b3db,
           "opentargets": get_open_targets, "chembl": get_chembl}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--only", nargs="+", choices=SOURCES, help="download only these sources")
    ap.add_argument("--force", action="store_true", help="re-download even if files exist")
    args = ap.parse_args()

    cfg.ensure_dirs()
    session = make_session()
    for name in (args.only or SOURCES):
        SOURCES[name](session, args.force)
    write_manifest()
    print("\nDone.")


if __name__ == "__main__":
    main()
