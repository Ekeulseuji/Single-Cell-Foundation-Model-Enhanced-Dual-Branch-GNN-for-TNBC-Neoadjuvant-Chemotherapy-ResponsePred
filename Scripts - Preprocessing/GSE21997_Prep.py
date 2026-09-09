import pandas as pd
import numpy as np
import os
import gzip

download_dir = '/Users/ekeulseuji/Downloads'

series_path = os.path.join(download_dir, 'GSE21997-GPL1390_series_matrix.txt.gz')

if not os.path.exists(series_path):
    for platform in ['GPL5325', 'GPL7504']:
        alt_path = os.path.join(download_dir, f'GSE21997-{platform}_series_matrix.txt.gz')
        if os.path.exists(alt_path):
            series_path = alt_path
            print(f"Using platform: {platform}")
            break

print(f"Using file: {series_path}")

print("\nParsing Series Matrix file...")

with gzip.open(series_path, 'rt') as f:
    lines = f.readlines()

table_start = None
table_end = None
header_lines = []

for i, line in enumerate(lines):
    if line.startswith('!series_matrix_table_begin'):
        table_start = i + 1
    elif line.startswith('!series_matrix_table_end'):
        table_end = i
        break
    elif line.startswith('!') and not line.startswith('!series_matrix_table'):
        header_lines.append(line.strip())

if table_start is None or table_end is None:
    raise ValueError("Data table markers not found")

data_lines = lines[table_start:table_end]

header_raw = data_lines[0].strip().split('\t')
header = [h.strip('"') for h in header_raw]

if header[0] != 'ID_REF':
    raise ValueError(f"First column is not ID_REF, got: {header[0]}")

sample_cols = header[1:]
print(f"Found {len(sample_cols)} samples")

data_rows = []
probe_ids = []

for line in data_lines[1:]:
    parts = line.strip().split('\t')
    parts_clean = [p.strip('"') for p in parts]
    if len(parts_clean) >= len(header):
        probe_ids.append(parts_clean[0])
        row_values = []
        for val in parts_clean[1:]:
            if val in ['null', 'NA', '']:
                row_values.append(np.nan)
            else:
                try:
                    row_values.append(float(val))
                except ValueError:
                    row_values.append(np.nan)
        data_rows.append(row_values)

print(f"Found {len(data_rows)} probes")

expr_matrix = np.array(data_rows)
expr = pd.DataFrame(expr_matrix, index=probe_ids, columns=sample_cols).T
print(f"Expression matrix dimensions: {expr.shape}")

expr_values = expr.values.flatten()
expr_values = expr_values[~np.isnan(expr_values)]
if len(expr_values) > 0:
    print(f"Expression range: {expr_values.min():.3f} ~ {expr_values.max():.3f}")

print("\nExtracting clinical information...")

char_lines = [line for line in header_lines if 'characteristics_ch2' in line]
print(f"Found {len(char_lines)} characteristics_ch2 lines")

char_dict = {}

for line in char_lines:
    parts = line.split('\t')
    if len(parts) > 1:
        first_val = parts[1].strip('"')
        if ': ' in first_val:
            feature_name = first_val.split(': ')[0]
        else:
            print(f"Warning: Cannot extract feature name from first value: {first_val}")
            continue
    else:
        print(f"Warning: Line has no values: {line[:50]}")
        continue
    
    values = []
    for i in range(1, len(parts)):
        val = parts[i].strip('"')
        if ': ' in val:
            val = val.split(': ')[1]
        values.append(val)
    
    if len(values) >= len(sample_cols):
        char_dict[feature_name] = values[:len(sample_cols)]
    else:
        char_dict[feature_name] = values + [np.nan] * (len(sample_cols) - len(values))

print(f"Successfully parsed {len(char_dict)} clinical features")

clinical_data = []
for i, gsm in enumerate(sample_cols):
    sample_info = {'sample_id': gsm}
    for key, values in char_dict.items():
        sample_info[key] = values[i] if i < len(values) else np.nan
    clinical_data.append(sample_info)

clinical = pd.DataFrame(clinical_data)
if 'sample_id' in clinical.columns:
    clinical.set_index('sample_id', inplace=True)

print(f"Clinical data dimensions: {clinical.shape}")
print("\nClinical data columns:")
print(list(clinical.columns))

print("\n========== Extracting ER/PR/HER2 Status ==========")

er_col = 'er (0=negative; 1=positive; 9 = n/a)'
pgr_col = 'pgr  (0=negative; 1=positive; 9 = n/a)'
her2_col = 'her2 (0=negative; 1=positive; 9 = n/a)'

if er_col in clinical.columns:
    clinical['er_status'] = pd.to_numeric(clinical[er_col], errors='coerce')
    print(f"  Extracted ER status, non-NA count: {clinical['er_status'].notna().sum()}")
else:
    print(f"  Column not found: {er_col}")
    clinical['er_status'] = np.nan

if pgr_col in clinical.columns:
    clinical['pgr_status'] = pd.to_numeric(clinical[pgr_col], errors='coerce')
    print(f"  Extracted PR status, non-NA count: {clinical['pgr_status'].notna().sum()}")
else:
    print(f"  Column not found: {pgr_col}")
    clinical['pgr_status'] = np.nan

if her2_col in clinical.columns:
    clinical['her2_status'] = pd.to_numeric(clinical[her2_col], errors='coerce')
    print(f"  Extracted HER2 status, non-NA count: {clinical['her2_status'].notna().sum()}")
else:
    print(f"  Column not found: {her2_col}")
    clinical['her2_status'] = np.nan

print("\nFirst 5 samples ER/PR/HER2 status:")
print(clinical[['er_status', 'pgr_status', 'her2_status']].head())

other_mapping = {
    'study': 'study',
    'age': 'age',
    'grade': 'histologic grade (1=grade i (low);  2= grade ii (intermediate); 3= grade iii (high); 9 = n/a)',
    'histology': 'histology (1=necrosis; 2=ductal carcinoma; 3=lobular; 4=mixed ductal/lobular carcinoma; 5=other  6=no invasive tumor present; 9 = n/a )',
    'tumor_size_cm': 'tumor size (pre-chemotherapy), cm (-1 = n/a)',
    'tumor_size_category': 'tumor size (-1 = n/a;  1 = <= 2cm; 2 = 2-5 cm; 3 = > 5cm)',
    'chemotherapy': 'neoadjuvant chemotherapy',
    'positive_nodes': 'number of positive nodes (post chemotherapy)',
    'pcr_raw': 'pcr',
    'rcb_class': 'residual cancer burden index class',
    'pam50_subtype': 'pam50 + claudin-low'
}

for new_col, old_col in other_mapping.items():
    if old_col in clinical.columns:
        if new_col in ['age', 'tumor_size_cm', 'tumor_size_category', 'positive_nodes']:
            clinical[new_col] = pd.to_numeric(clinical[old_col], errors='coerce')
        else:
            clinical[new_col] = clinical[old_col]
        print(f"  Mapped: '{old_col}' -> {new_col}")
    else:
        print(f"  Column not found: {old_col}")
        clinical[new_col] = np.nan

if 'pcr_raw' in clinical.columns:
    clinical['pcr'] = clinical['pcr_raw'].apply(lambda x: 1 if str(x).upper().strip() == 'YES' else (0 if str(x).upper().strip() == 'NO' else np.nan))
    print(f"  Converted pcr: Yes -> 1, No -> 0")
else:
    clinical['pcr'] = np.nan

print("\n========== Marking TNBC Samples ==========")

print("\nER status distribution:")
print(clinical['er_status'].value_counts(dropna=False))
print("\nPR status distribution:")
print(clinical['pgr_status'].value_counts(dropna=False))
print("\nHER2 status distribution:")
print(clinical['her2_status'].value_counts(dropna=False))

clinical['is_tnbc'] = 0
tnbc_mask = ((clinical['er_status'] == 0) &
             (clinical['pgr_status'] == 0) &
             (clinical['her2_status'] == 0))
clinical.loc[tnbc_mask, 'is_tnbc'] = 1

print(f"\nTotal samples: {len(clinical)}")
print(f"TNBC samples: {clinical['is_tnbc'].sum()}")
print(f"Non-TNBC samples: {len(clinical) - clinical['is_tnbc'].sum()}")

if 'pcr' in clinical.columns:
    print("\nAll samples pCR distribution:")
    print(clinical['pcr'].value_counts(dropna=False))

tnbc_samples = clinical[clinical['is_tnbc'] == 1].index.tolist()
tnbc_samples_in_expr = [s for s in tnbc_samples if s in expr.index]

print(f"\nTNBC samples in expression matrix: {len(tnbc_samples_in_expr)}")

if len(tnbc_samples_in_expr) == 0:
    print("Warning: No TNBC samples found, using all samples...")
    tnbc_samples_in_expr = expr.index.tolist()
    clinical['is_tnbc'] = 1
else:
    print(f"Found {len(tnbc_samples_in_expr)} TNBC samples")

tnbc_expr = expr.loc[tnbc_samples_in_expr].copy()
tnbc_expr_gene = tnbc_expr.T
print(f"TNBC expression matrix dimensions: {tnbc_expr_gene.shape}")

tnbc_clinical = clinical.loc[tnbc_samples_in_expr].copy()

def map_pcr_label(pcr_value):
    if pd.isna(pcr_value):
        return np.nan
    if pcr_value == 1:
        return 'pCR'
    elif pcr_value == 0:
        return 'RD'
    else:
        return np.nan

if 'pcr' in tnbc_clinical.columns:
    tnbc_clinical['pcr_rd'] = tnbc_clinical['pcr'].apply(map_pcr_label)
else:
    tnbc_clinical['pcr_rd'] = np.nan

print(f"\nTNBC pCR/RD label distribution:")
print(tnbc_clinical['pcr_rd'].value_counts(dropna=False))

print("\nTNBC clinical preview:")
preview_cols = ['pam50_subtype', 'pcr', 'pcr_rd', 'er_status', 'pgr_status', 'her2_status', 'age', 'grade']
preview_cols = [c for c in preview_cols if c in tnbc_clinical.columns]
if preview_cols:
    print(tnbc_clinical[preview_cols].head(10))

out_clinical = os.path.join(download_dir, 'GSE21997_all_clinical_data.csv')
clinical.to_csv(out_clinical, encoding='utf-8')
print(f"\nSaved: {out_clinical}")

out_expr = os.path.join(download_dir, 'GSE21997_TNBC_expression.csv')
tnbc_expr_gene.to_csv(out_expr, encoding='utf-8')
print(f"Saved: {out_expr}")

labels_cols = ['pcr_rd']
if labels_cols:
    labels_df = tnbc_clinical[labels_cols].copy()
    out_labels = os.path.join(download_dir, 'GSE21997_TNBC_labels.csv')
    labels_df.to_csv(out_labels, encoding='utf-8')
    print(f"Saved: {out_labels}")

print("\n========== Preprocessing Complete ==========")
print(f"   Total samples: {len(clinical)}")
print(f"   TNBC samples: {len(tnbc_clinical)}")
print(f"   Features: {tnbc_expr_gene.shape[0]}")
print(f"   pCR samples: {tnbc_clinical['pcr_rd'].value_counts().get('pCR', 0)}")
print(f"   RD samples: {tnbc_clinical['pcr_rd'].value_counts().get('RD', 0)}")

print("\nFinal TNBC samples ER/PR/HER2 confirmation:")
if len(tnbc_clinical) > 0:
    print(tnbc_clinical[['er_status', 'pgr_status', 'her2_status']].head(10))
