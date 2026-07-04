#!/usr/bin/env python3
"""
retrieve_expression.py — fetch rat IN-VIVO LIVER expression (GPL1355) for the cohort.

Option A (clean liver, single platform):
  - DrugMatrix Affymetrix liver = GSE57815  (this file: fetch_drugmatrix)
  - TG-GATEs liver in vivo       = E-MTAB-799 / E-MTAB-800  (fetch_tggates, task #5)

Outputs (data/expression/):
  drugmatrix_liver_expr.parquet   probes x samples, float32 (normalised intensities, GEO VALUE)
  drugmatrix_liver_manifest.csv   one row per sample: compound/dose/time/vehicle/tissue/
                                  is_control + connectivity + in_labelled flag
  drugmatrix_liver_coverage.txt   human-readable coverage summary

Both sources are Affymetrix Rat 230 2.0 -> same probe space -> ready for
per-compound logFC (treated - matched vehicle) -> ComBat -> PCA downstream.
"""
import os, re, gzip
import pandas as pd

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT   = os.path.join(ROOT, "data", "expression"); os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(ROOT, "data", "_raw");        os.makedirs(CACHE, exist_ok=True)
GSE57815_GZ = os.path.join(CACHE, "GSE57815_series_matrix.txt.gz")

def norm(s):
    return re.sub(r"\s+", " ", str(s).strip().upper())

def build_name_index():
    """normalised compound name -> InChIKey connectivity, from dm_keyed (+ cohort fallback)."""
    coh = pd.read_csv(os.path.join(ROOT, "master_cohort.csv"), dtype=str)
    dmk = pd.read_csv(os.path.join(ROOT, "data", "dm_keyed.csv"), dtype=str)
    name2conn = {norm(r["name"]): r["connectivity"] for _, r in dmk.iterrows()}
    for _, r in coh.iterrows():
        name2conn.setdefault(norm(r["compound_name"]), r["connectivity"])
    labelled = set(coh.loc[coh.in_tox21.eq("Y"), "connectivity"])
    cohort   = set(coh.connectivity)
    return name2conn, labelled, cohort

# --------- series-matrix parsing (no GEOparse; we know the exact format) --------
def _parse_series_matrix_header(path):
    """Return per-sample manifest rows from the !Sample_ header block, and the
    line index of !series_matrix_table_begin (0-based)."""
    sample_lines = {}   # tag -> list of value-lists (chars appears many times)
    gsm_ids, titles = [], []
    table_begin = None
    with gzip.open(path, "rt", encoding="latin-1") as f:
        for i, line in enumerate(f):
            if line.startswith("!series_matrix_table_begin"):
                table_begin = i; break
            if not line.startswith("!Sample_"): continue
            p = line.rstrip("\n").split("\t")
            tag = p[0]; vals = [x.strip().strip('"') for x in p[1:]]
            if tag == "!Sample_geo_accession": gsm_ids = vals
            elif tag == "!Sample_title":       titles = vals
            elif tag == "!Sample_characteristics_ch1":
                sample_lines.setdefault("char", []).append(vals)
    return gsm_ids, titles, sample_lines.get("char", []), table_begin

# ============================ DrugMatrix (GSE57815) ============================
def fetch_drugmatrix():
    name2conn, labelled, cohort = build_name_index()

    print("Parsing GSE57815 series-matrix header…", flush=True)
    gsm_ids, titles, char_lines, table_begin = _parse_series_matrix_header(GSE57815_GZ)
    n = len(gsm_ids)

    def char_dict(i):
        d = {}
        for line in char_lines:
            v = line[i] if i < len(line) else ""
            if ":" in v:
                k, val = v.split(":", 1); d[k.strip().lower()] = val.strip()
        return d

    rows = []
    for i in range(n):
        ch = char_dict(i)
        compound = ch.get("compound", "")
        dose     = ch.get("dose", "")
        is_ctrl  = (compound == "") or dose.startswith("0 ")
        conn     = name2conn.get(norm(compound), "") if compound else ""
        rows.append({
            "sample": gsm_ids[i],
            "title":  titles[i] if i < len(titles) else "",
            "compound": compound, "dose": dose,
            "time": ch.get("time", ""), "vehicle": ch.get("vehicle", ""),
            "tissue": ch.get("tissue", "Liver"),
            "is_control": is_ctrl,
            "connectivity": conn,
            "in_cohort":   conn in cohort,
            "in_labelled": conn in labelled,
            "dataset": "drugmatrix_liver",
        })
    man = pd.DataFrame(rows)

    # keep: controls (all — needed for logFC) + treated samples of cohort compounds
    keep = man["is_control"] | man["in_cohort"]
    man  = man[keep].reset_index(drop=True)

    # --- expression table: rows after table_begin; header row = ID_REF + GSMs ---
    print("Parsing expression matrix (2218 cols)…", flush=True)
    expr = pd.read_csv(GSE57815_GZ, sep="\t", skiprows=table_begin + 1,
                       quotechar='"', na_values=["", "null"], engine="c")
    expr = expr[~expr.iloc[:, 0].astype(str).str.startswith("!")]  # drop table_end marker
    expr = expr.set_index(expr.columns[0])
    expr.index.name = "probe"
    expr = expr.apply(pd.to_numeric, errors="coerce").astype("float32")
    expr = expr[man["sample"].tolist()]                            # subset to kept samples

    expr.to_parquet(os.path.join(OUT, "drugmatrix_liver_expr.parquet"))
    man.to_csv(os.path.join(OUT, "drugmatrix_liver_manifest.csv"), index=False)

    # --- coverage summary ---
    n_labelled_covered = man.loc[man.in_labelled, "connectivity"].nunique()
    summary = (
        f"GSE57815 DrugMatrix Affymetrix liver\n"
        f"  total samples parsed : {len(rows)}\n"
        f"  kept (cohort+ctrl)   : {len(man)}\n"
        f"  vehicle controls     : {int(man.is_control.sum())}\n"
        f"  probes               : {expr.shape[0]}\n"
        f"  distinct treated cpd : {man.loc[~man.is_control,'compound'].nunique()}\n"
        f"  -> labelled compounds covered : {n_labelled_covered} / 613\n"
    )
    with open(os.path.join(OUT, "drugmatrix_liver_coverage.txt"), "w") as f:
        f.write(summary)
    print(summary, flush=True)
    print(f"wrote {OUT}/drugmatrix_liver_expr.parquet  ({expr.shape[0]} x {expr.shape[1]})", flush=True)

if __name__ == "__main__":
    fetch_drugmatrix()
