# ============================================================
# GSE32603 preprocessing: GPR -> Cy5 -> GB_ACC -> gene symbol
# Output: genes x samples
# ============================================================
import os, gc
import numpy as np
import pandas as pd
import GEOparse
from tqdm.auto import tqdm

GPR_DIR = "/Users/ekeulseuji/Downloads/GSE32603_RAW"
OUT_DIR = "/Users/ekeulseuji/Downloads"
LABEL_FILE = "/Users/ekeulseuji/Downloads/GSE32603_TNBC_labels.csv"

print("=" * 70)
print("GSE32603 preprocessing")
print("=" * 70)

# ---------- 1. Load GPL14668 annotation ----------
print("\n[1] Load GPL14668")
gpl = GEOparse.get_GEO(geo="GPL14668", destdir="./geo_cache", silent=True)
tbl = gpl.table.copy()
for c in tbl.columns:
    tbl[c] = tbl[c].astype(str).str.strip().str.strip('"')

def clean_symbol(s):
    s = str(s).upper().strip()
    if s in ("", "NAN", "---", "N/A", "NONE", "NULL"):
        return None
    for sep in ["///", "//", ";", ",", "|"]:
        if sep in s:
            s = s.split(sep)[0].strip()
    return s if s else None

tbl["SYM"] = tbl["GENE SYMBOL"].apply(clean_symbol)

gb_acc_to_sym = {}
for _, row in tbl.iterrows():
    acc = str(row["GB_ACC"]).strip().upper()
    sym = row["SYM"]
    if acc and sym and acc not in gb_acc_to_sym:
        gb_acc_to_sym[acc] = sym

print(f"  GB_ACC -> SYM: {len(gb_acc_to_sym)} mappings")

# ---------- 2. Read GPR files ----------
def read_gpr(fpath):
    with open(fpath, "r", encoding="latin1") as f:
        for i, line in enumerate(f):
            if line.lstrip().startswith('"Block"'):
                hdr = i
                break
    df = pd.read_csv(fpath, skiprows=hdr, sep="\t",
                     encoding="latin1", quoting=3,
                     on_bad_lines="skip", low_memory=False)
    df.columns = [c.strip().strip('"') for c in df.columns]
    return df

gpr_files = sorted([f for f in os.listdir(GPR_DIR) if f.lower().endswith(".gpr")])
print(f"\n[2] Found {len(gpr_files)} GPR files")

ref_df = read_gpr(os.path.join(GPR_DIR, gpr_files[0]))
probe_ids = ref_df["ID"].astype(str).str.strip().str.strip('"').values
n_probes = len(probe_ids)
print(f"  Probes: {n_probes}")

# ---------- 3. Extract Cy5 single-channel intensity ----------
print(f"\n[3] Extract Cy5 (635nm) net intensity")
sample_names = []
expr_matrix = np.zeros((n_probes, len(gpr_files)), dtype=np.float32)

for i, fname in enumerate(tqdm(gpr_files, desc="  Reading")):
    df = read_gpr(os.path.join(GPR_DIR, fname))

    fg = pd.to_numeric(df["F635 Mean"], errors="coerce").values
    bg = pd.to_numeric(df["B635 Mean"], errors="coerce").values
    intensity = fg - bg

    intensity = np.nan_to_num(intensity, nan=0.0, posinf=0.0, neginf=0.0)
    intensity[intensity < 0] = 0

    file_ids = df["ID"].astype(str).str.strip().str.strip('"').values
    if len(intensity) != n_probes or not (file_ids == probe_ids).all():
        id_to_int = dict(zip(file_ids, intensity))
        intensity = np.array([id_to_int.get(p, 0.0) for p in probe_ids], dtype=np.float32)

    expr_matrix[:, i] = intensity
    sample_names.append(fname.replace(".gpr", ""))

print(f"  Raw matrix: {expr_matrix.shape}")
print(f"  Range: [{expr_matrix.min():.2f}, {expr_matrix.max():.2f}]")
print(f"  Per-sample sum: [{expr_matrix.sum(0).min():.0f}, {expr_matrix.sum(0).max():.0f}]")

# ---------- 4. log2 transform ----------
print(f"\n[4] log2(x+1) transform")
if expr_matrix.max() > 100:
    expr_matrix = np.log2(expr_matrix + 1)
    print(f"  After log2: [{expr_matrix.min():.3f}, {expr_matrix.max():.3f}]")
    print(f"  Mean: {expr_matrix.mean():.3f}")

# ---------- 5. Map to gene symbol ----------
print(f"\n[5] GB_ACC -> gene symbol")
probe_syms = [gb_acc_to_sym.get(p.upper(), None) for p in probe_ids]
valid_mask = [s is not None for s in probe_syms]
n_valid = sum(valid_mask)
print(f"  Matched probes: {n_valid} / {n_probes} ({n_valid/n_probes*100:.1f}%)")

expr_genes = expr_matrix[valid_mask, :]
gene_ids_valid = [s for s in probe_syms if s is not None]

expr_df = pd.DataFrame(expr_genes, index=gene_ids_valid, columns=sample_names)
expr_df = expr_df.groupby(level=0).mean()

print(f"  Gene matrix (genes x samples): {expr_df.shape}")
print(f"  Per-sample sum: [{expr_df.sum(axis=0).min():.1f}, {expr_df.sum(axis=0).max():.1f}]")
print(f"  Mean: {expr_df.values.mean():.3f}")
print(f"  Negative: {int((expr_df.values < 0).sum())}")

# ---------- 6. Align with labels ----------
print(f"\n[6] Align with labels")
lab = pd.read_csv(LABEL_FILE)
common_samples = sorted(set(expr_df.columns) & set(lab["sample_id"]))
print(f"  Common samples: {len(common_samples)} / {len(lab)}")

if len(common_samples) < len(lab) * 0.8:
    print(f"  WARNING: intersection too small, check GPR filenames vs labels")
    print(f"     GPR first 5:    {sample_names[:5]}")
    print(f"     Labels first 5: {list(lab['sample_id'].head(5))}")

if len(common_samples) > 0:
    expr_final = expr_df[common_samples]
    labels_final = lab.set_index("sample_id").loc[common_samples].reset_index()

    expr_path = os.path.join(OUT_DIR, "GSE32603_TNBC_expression_FIXED.csv")
    lab_path = os.path.join(OUT_DIR, "GSE32603_TNBC_labels_FIXED.csv")

    expr_final.to_csv(expr_path)
    labels_final.to_csv(lab_path, index=False)

    print(f"\n  Saved: {expr_path}  shape={expr_final.shape}")
    print(f"  Saved: {lab_path}  shape={labels_final.shape}")

    s_col = expr_final.sum(axis=0)
    print(f"\n  Final checks:")
    print(f"    per-sample sum: [{s_col.min():.1f}, {s_col.max():.1f}]")
    print(f"    negative: {int((expr_final.values < 0).sum())}")
    print(f"    NaN: {int(expr_final.isna().sum().sum())}")
    print(f"    pCR={(labels_final.pcr_rd=='pCR').sum()}, RD={(labels_final.pcr_rd=='RD').sum()}")

    try:
        ref = pd.read_csv("/Users/ekeulseuji/Downloads/GSE25066_TNBC_expression.csv", index_col=0)
        print(f"\n  Compare with GSE25066:")
        print(f"    GSE25066: mean={ref.values.mean():.3f}, "
              f"per-sample sum=[{ref.sum(0).min():.1f}, {ref.sum(0).max():.1f}]")
        print(f"    GSE32603: mean={expr_final.values.mean():.3f}, "
              f"per-sample sum=[{s_col.min():.1f}, {s_col.max():.1f}]")
        ratio = expr_final.values.mean() / ref.values.mean()
        print(f"    mean ratio: {ratio:.2f}  (ideal 0.5~2.0)")
    except Exception as e:
        print(f"    comparison failed: {e}")
