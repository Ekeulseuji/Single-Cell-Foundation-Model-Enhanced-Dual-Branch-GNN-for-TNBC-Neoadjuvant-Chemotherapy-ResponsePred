# ============================================================
# GSE22358 preprocessing: SOFT -> CH2DL_MEAN -> gene symbol
# Output: genes x samples
# ============================================================
import gzip
import re
import os
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

SOFT_PATH = "/Users/ekeulseuji/Downloads/GSE22358_family.soft.gz"
OUT_DIR = "/Users/ekeulseuji/Downloads"

print("=" * 70)
print("GSE22358 preprocessing")
print("=" * 70)

# ---------- 1. Read SOFT ----------
print(f"\n[1] Read {SOFT_PATH}")
with gzip.open(SOFT_PATH, "rt", encoding="latin1") as f:
    lines = f.readlines()
print(f"  Total lines: {len(lines)}")

gsm_starts = []
for i, line in enumerate(lines):
    if line.startswith("^SAMPLE"):
        m = re.match(r"\^SAMPLE\s*=\s*(GSM\d+)", line)
        if m:
            gsm_starts.append((i, m.group(1)))

print(f"  GSM count: {len(gsm_starts)}")
print(f"  First 5: {[g[1] for g in gsm_starts[:5]]}")

# ---------- 2. Parse each GSM table ----------
FULL_COLS = [
    "ID_REF", "VALUE", "SPOT", "CH1_MEAN", "CH1_SD",
    "CH1_BKD_MEDIAN", "CH1_BKD_SD", "CH2_MEAN", "CH2_SD",
    "CH2_BKD_MEDIAN", "CH2_BKD_SD", "TOT_BPIX", "TOT_SPIX",
    "CH2BN_MEDIAN", "CH2IN_MEAN", "CH1DL_MEAN", "CH2DL_MEAN",
    "LOG_RAT2N_MEAN", "CORR", "FLAG", "UNF_VALUE"
]

def parse_gsm_table(lines, start_idx):
    table_begin = None
    table_end = None
    for i in range(start_idx + 1, len(lines)):
        line = lines[i]
        if "!sample_table_begin" in line:
            table_begin = i + 1
        elif "!sample_table_end" in line:
            table_end = i
            break
        elif line.startswith("^SAMPLE"):
            break

    if table_begin is None or table_end is None:
        return None

    data = []
    for i in range(table_begin + 1, table_end):
        parts = lines[i].rstrip("\n").split("\t")
        if len(parts) < len(FULL_COLS):
            continue
        data.append(parts[:len(FULL_COLS)])

    if not data:
        return None
    return pd.DataFrame(data, columns=FULL_COLS)


print(f"\n[2] Parse {len(gsm_starts)} GSMs")
sample_ids = []
expr_matrices = []
probe_ids_ref = None

for idx, (start_i, gsm_id) in enumerate(tqdm(gsm_starts, desc="  Parsing")):
    df = parse_gsm_table(lines, start_i)

    if df is None:
        print(f"    WARNING {gsm_id}: parse failed")
        continue

    if "ID_REF" not in df.columns or "CH2DL_MEAN" not in df.columns:
        print(f"    WARNING {gsm_id}: missing columns, got {list(df.columns)}")
        continue

    ids = df["ID_REF"].astype(str).values
    cy5 = pd.to_numeric(df["CH2DL_MEAN"], errors="coerce").values

    if probe_ids_ref is None:
        probe_ids_ref = ids

    if len(ids) != len(probe_ids_ref):
        print(f"    WARNING {gsm_id}: length mismatch")
        continue

    sample_ids.append(gsm_id)
    expr_matrices.append(cy5)

print(f"\nParsed {len(sample_ids)} samples")
print(f"  Probes: {len(probe_ids_ref) if probe_ids_ref is not None else 0}")

if len(sample_ids) == 0:
    print("\nNo samples parsed. Dump first GSM diagnostics:")
    first_start = gsm_starts[0][0]
    for i in range(first_start, min(first_start + 100, len(lines))):
        print(f"    {i:>6} | {lines[i][:120].rstrip()}")
    raise SystemExit("Diagnostics printed above.")

# ---------- 3. Build matrix ----------
expr_matrix = np.stack(expr_matrices, axis=1)
expr_df = pd.DataFrame(expr_matrix, index=probe_ids_ref, columns=sample_ids)
print(f"\n[3] Matrix: {expr_df.shape}")
print(f"  Range: [{expr_df.values.min():.2f}, {expr_df.values.max():.2f}]")
print(f"  Mean: {expr_df.values.mean():.2f}")

# ---------- 4. log2 ----------
print(f"\n[4] log2(x+1)")
expr_log = np.log2(expr_df.clip(lower=0) + 1)
print(f"  After log2: [{expr_log.values.min():.3f}, {expr_log.values.max():.3f}]")
print(f"  Mean: {expr_log.values.mean():.3f}")
print(f"  Per-sample sum: [{expr_log.sum(0).min():.1f}, {expr_log.sum(0).max():.1f}]")

# ---------- 5. Map to gene symbol ----------
print(f"\n[5] Map to gene symbol")
import GEOparse

gpl = GEOparse.get_GEO(geo="GPL5325", destdir="./geo_cache", silent=True)
tbl = gpl.table.copy()
for c in tbl.columns:
    tbl[c] = tbl[c].astype(str).str.strip().str.strip('"')

print(f"  GPL columns: {list(tbl.columns)}")
print(f"  ID_REF first 5:      {list(probe_ids_ref[:5])}")
print(f"  Probe Name first 5:  {list(tbl['Probe Name'].head(5)) if 'Probe Name' in tbl.columns else 'N/A'}")
print(f"  ID first 5:          {list(tbl['ID'].head(5)) if 'ID' in tbl.columns else 'N/A'}")

first_id = str(probe_ids_ref[0]).strip().upper()
if first_id.isdigit() and "ID" in tbl.columns:
    id_col = "ID"
elif first_id.startswith("A_"):
    id_col = "Probe Name"
elif first_id[:2] in ("NM", "NR", "AA", "AI"):
    id_col = "GB_ACC"
else:
    id_col = "Probe Name"

sym_col = "Blast Gene Symbol"
print(f"  Using: {id_col} -> {sym_col}")

probe_to_sym = {}
for _, row in tbl.iterrows():
    key = str(row[id_col]).strip().upper()
    sym = str(row[sym_col]).strip().upper()
    if key and key not in ("NAN", "") and sym and sym not in ("NAN", ""):
        for sep in ["///", "//", ";", ",", "|"]:
            if sep in sym:
                sym = sym.split(sep)[0].strip()
        if sym:
            probe_to_sym[key] = sym

print(f"  {id_col} -> SYM: {len(probe_to_sym)} mappings")

gene_syms = [probe_to_sym.get(str(p).strip().upper(), None) for p in probe_ids_ref]
valid_mask = [s is not None for s in gene_syms]
print(f"  Matched: {sum(valid_mask)}/{len(probe_ids_ref)}")

expr_genes = expr_log.values[valid_mask, :]
gene_ids_valid = [s for s in gene_syms if s is not None]

expr_gene_df = pd.DataFrame(expr_genes, index=gene_ids_valid, columns=sample_ids)
expr_gene_df = expr_gene_df.groupby(level=0).mean()
print(f"  Gene matrix (genes x samples): {expr_gene_df.shape}")

# ---------- 6. Align with labels ----------
print(f"\n[6] Align with labels")
lab = pd.read_csv(f"{OUT_DIR}/GSE22358_TNBC_labels.csv")
print(f"  label columns: {list(lab.columns)}")
print(f"  label head:\n{lab.head(3)}")
common = sorted(set(expr_gene_df.columns) & set(lab["sample_id"]))
print(f"  Common: {len(common)}/{len(lab)}")

if len(common) < len(lab) * 0.8:
    print(f"  WARNING: intersection too small")
    print(f"     SOFT first 5:   {list(expr_gene_df.columns[:5])}")
    print(f"     Labels first 5: {list(lab['sample_id'].head(5))}")

# ---------- 7. Save ----------
if len(common) > 0:
    expr_final = expr_gene_df[common]
    labels_final = lab.set_index("sample_id").loc[common].reset_index()

    expr_final.to_csv(f"{OUT_DIR}/GSE22358_TNBC_expression_FIXED.csv")
    labels_final.to_csv(f"{OUT_DIR}/GSE22358_TNBC_labels_FIXED.csv", index=False)

    s_col = expr_final.sum(axis=0)
    print(f"\n[7] Final:")
    print(f"  Shape: {expr_final.shape}  (genes x samples)")
    print(f"  per-sample sum: [{s_col.min():.1f}, {s_col.max():.1f}]")
    print(f"  negative: {int((expr_final.values < 0).sum())}")
    print(f"  NaN: {int(expr_final.isna().sum().sum())}")
    print(f"  Mean: {expr_final.values.mean():.3f}")
    print(f"  pCR={(labels_final.pcr_rd=='pCR').sum()}, RD={(labels_final.pcr_rd=='RD').sum()}")

    ref = pd.read_csv(f"{OUT_DIR}/GSE25066_TNBC_expression.csv", index_col=0)
    ratio = expr_final.values.mean() / ref.values.mean()
    print(f"  mean ratio vs GSE25066: {ratio:.2f}")

    print(f"\nSaved: {OUT_DIR}/GSE22358_TNBC_expression_FIXED.csv")
