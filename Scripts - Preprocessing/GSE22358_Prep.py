import pandas as pd
import numpy as np
import os
import gzip

download_dir = '/Users/ekeulseuji/Downloads'
series_path = os.path.join(download_dir, 'GSE22358_series_matrix.txt.gz')

print("Parsing Series Matrix file...")

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

print("\nClinical data first 5 rows:")
print(clinical.head())

field_mapping = {
    'er (0 = negative; 1 = positive)': 'er_status',
    'pgr (0 = negative; 1 = positive)': 'pgr_status',
    'her2 (0 = negative; 1 = positive)': 'her2_status',
    'intrinsic subtypes by pam50': 'pam50_subtype',
    'p53 status by amplichip': 'p53_amplichip',
    'p53 measured by ihc': 'p53_ihc',
    'response': 'response',
    'tumor size (cm)': 'tumor_size_cm',
    'grade': 'grade',
    'neoadjuvant chemotherapy': 'chemotherapy',
    'slidename': 'slide_name',
    'sample': 'sample_name',
    'study': 'study'
}

for old_key, new_key in field_mapping.items():
    if old_key in clinical.columns:
        clinical[new_key] = clinical[old_key]
        print(f"  Mapped: {old_key} -> {new_key}")
    else:
        matched_cols = [col for col in clinical.columns if old_key.lower() in col.lower()]
        if matched_cols:
            clinical[new_key] = clinical[matched_cols[0]]
            print(f"  Fuzzy matched: {matched_cols[0]} -> {new_key}")
        else:
            clinical[new_key] = np.nan

numeric_fields = ['er_status', 'pgr_status', 'her2_status', 'tumor_size_cm', 'grade']
for field in numeric_fields:
    if field in clinical.columns:
        clinical[field] = pd.to_numeric(clinical[field], errors='coerce')

clinical['is_tnbc'] = ((clinical['er_status'] == 0) &
                       (clinical['pgr_status'] == 0) &
                       (clinical['her2_status'] == 0)).astype(int)

print(f"\nTotal samples: {len(clinical)}")
print(f"TNBC samples: {clinical['is_tnbc'].sum()}")

if 'response' in clinical.columns:
    print("\nAll samples Response distribution:")
    print(clinical['response'].value_counts(dropna=False))

tnbc_samples = clinical[clinical['is_tnbc'] == 1].index.tolist()
tnbc_samples_in_expr = [s for s in tnbc_samples if s in expr.index]

print(f"\nTNBC samples in expression matrix: {len(tnbc_samples_in_expr)}")

if len(tnbc_samples_in_expr) == 0:
    print("Warning: No TNBC samples found, using all samples...")
    tnbc_samples_in_expr = expr.index.tolist()
    clinical['is_tnbc'] = 1

tnbc_expr = expr.loc[tnbc_samples_in_expr].copy()
tnbc_expr_gene = tnbc_expr.T
print(f"TNBC expression matrix (before mapping): {tnbc_expr_gene.shape}")

print("\nMapping probe IDs to Gene Symbols...")

def map_probes_to_symbols(expr_df, gpl_id='GPL5325'):
    try:
        import GEOparse
    except ImportError:
        raise ImportError("Please install GEOparse: pip install GEOparse")
    
    print(f"  Downloading {gpl_id} annotation from GEO...")
    gpl = GEOparse.get_GEO(geo=gpl_id, destdir="geo_cache", silent=True)
    
    id_col = "ID" if "ID" in gpl.table.columns else gpl.table.columns[0]
    
    symbol_col = None
    for col in gpl.table.columns:
        if "gene symbol" in col.lower() or "symbol" in col.lower():
            symbol_col = col
            break
    if symbol_col is None:
        for col in gpl.table.columns:
            if "gene" in col.lower() and "symbol" in col.lower():
                symbol_col = col
                break
    if symbol_col is None:
        raise ValueError(f"Gene Symbol column not found. Available columns: {list(gpl.table.columns)}")
    
    mapping = {}
    for _, row in gpl.table[[id_col, symbol_col]].iterrows():
        probe = str(row[id_col]).strip()
        symbol_raw = str(row[symbol_col]).strip().upper()
        if symbol_raw in ['', 'NA', 'N/A', '---', 'NONE', 'NULL']:
            continue
        if '///' in symbol_raw:
            symbol_raw = symbol_raw.split('///')[0].strip()
        elif '//' in symbol_raw:
            symbol_raw = symbol_raw.split('//')[0].strip()
        elif ';' in symbol_raw:
            symbol_raw = symbol_raw.split(';')[0].strip()
        symbol = symbol_raw.strip('"').strip()
        if symbol and symbol not in ['', '---']:
            mapping[probe] = symbol
    
    print(f"  Built {len(mapping)} probe-gene mappings")
    
    expr_df = expr_df.copy()
    expr_df.index = expr_df.index.astype(str).str.strip()
    new_index = [mapping.get(x, None) for x in expr_df.index]
    valid_mask = [x is not None for x in new_index]
    expr_df = expr_df.loc[valid_mask]
    expr_df.index = [x for x in new_index if x is not None]
    
    expr_df = expr_df.groupby(level=0).mean()
    
    mapped_count = len(expr_df)
    total_probes = expr_df.shape[0]
    print(f"  Successfully mapped {mapped_count}/{total_probes} probes ({mapped_count/total_probes*100:.1f}%)")
    
    if mapped_count / len(expr_df.index.unique()) < 0.5:
        print("  Warning: Mapping rate below 50%. Please check GPL version or probe ID compatibility.")
        print("  Consider manually downloading GPL5325 annotation file and verifying ID column.")
    
    return expr_df

tnbc_expr_gene = map_probes_to_symbols(tnbc_expr_gene, gpl_id='GPL5325')
print(f"Genes after mapping: {tnbc_expr_gene.shape[0]}")

tnbc_clinical = clinical.loc[tnbc_samples_in_expr].copy()

def map_pcr_label(response):
    if pd.isna(response):
        return np.nan
    resp = str(response).upper().strip()
    if any(kw in resp for kw in ['COMPLETE', 'PCR', 'NEAR']):
        return 'pCR'
    else:
        return 'RD'

if 'response' in tnbc_clinical.columns:
    tnbc_clinical['pcr_label'] = tnbc_clinical['response'].apply(map_pcr_label)
else:
    tnbc_clinical['pcr_label'] = np.nan

print(f"\nTNBC pCR/RD label distribution:")
print(tnbc_clinical['pcr_label'].value_counts(dropna=False))

print("\nTNBC clinical preview:")
preview_cols = ['pam50_subtype', 'response', 'pcr_label', 'er_status', 'pgr_status', 'her2_status']
preview_cols = [c for c in preview_cols if c in tnbc_clinical.columns]
if preview_cols:
    print(tnbc_clinical[preview_cols].head(10))

out_clinical = os.path.join(download_dir, 'GSE22358_all_clinical_data.csv')
clinical.to_csv(out_clinical, encoding='utf-8')
print(f"\nSaved: {out_clinical}")

out_expr = os.path.join(download_dir, 'GSE22358_TNBC_expression.csv')
tnbc_expr_gene.to_csv(out_expr, encoding='utf-8')
print(f"Saved: {out_expr}")

labels_cols = ['pcr_label']
labels_cols = [c for c in labels_cols if c in tnbc_clinical.columns]
if labels_cols:
    labels_df = tnbc_clinical[labels_cols].copy()
    out_labels = os.path.join(download_dir, 'GSE22358_TNBC_labels.csv')
    labels_df.to_csv(out_labels, encoding='utf-8')
    print(f"Saved: {out_labels}")

print("\n========== Preprocessing Complete ==========")
print(f"   Total samples: {len(clinical)}")
print(f"   TNBC samples: {len(tnbc_clinical)}")
print(f"   Features: {tnbc_expr_gene.shape[0]}")
