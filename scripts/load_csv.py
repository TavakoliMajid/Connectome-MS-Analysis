import os
import pandas as pd
import numpy as np

# =======================
# CONFIG
# =======================
base_dir = "C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj"  # root folder containing INsIDER_Cxxx / INsIDER_Pxxx
output_dir = os.path.join(base_dir, "Processed_Connectomes")
os.makedirs(output_dir, exist_ok=True)

# Whether to keep full matrix or only upper triangle (to avoid duplicate undirected edges)
use_upper_triangle_only = True
include_diagonal = False  # set True if you want to keep self-connections


# =======================
# HELPERS
# =======================
def load_connectome_csv(path: str) -> pd.DataFrame:
    """Load a CSV connectome (no headers) into a DataFrame."""
    return pd.read_csv(path, header=None)

def infer_group_from_subject(subj: str) -> str | None:
    """
    Infer group from subject ID:
      INsIDER_C### -> 'control'
      INsIDER_P### -> 'patient'
    Returns None if not recognizable.
    """
    try:
        # Expect something like 'INsIDER_C019' or 'INsIDER_P39'
        after_underscore = subj.split('_', 1)[1]  # 'C019' or 'P39'
        first = after_underscore[0].upper()
        if first == 'C':
            return 'control'
        if first == 'P':
            return 'patient'
    except Exception:
        pass
    return None

def flatten_matrix(df: pd.DataFrame,
                   upper_only: bool = True,
                   keep_diag: bool = False) -> (np.ndarray, np.ndarray, np.ndarray):
    """
    Flatten a square matrix to a 1D vector.
    Returns (values, row_idx, col_idx).
    If upper_only is True, returns upper-triangle (optionally including diag).
    """
    mat = df.values
    n = mat.shape[0]
    if upper_only:
        k = 0 if keep_diag else 1
        r, c = np.triu_indices(n, k=k)
        vals = mat[r, c]
        return vals, r, c
    else:
        r, c = np.indices(mat.shape)
        return mat.ravel(), r.ravel(), c.ravel()


# =======================
# SCAN & LOAD
# =======================
act_connectomes = {}
trekker_connectomes = {}
subjects_meta = []  # to save a small metadata table (subject, group, has_ACT, has_TREKKER)

for subj in sorted(os.listdir(base_dir)):
    subj_path = os.path.join(base_dir, subj)
    if not (os.path.isdir(subj_path) and subj.startswith("INsIDER")):
        continue

    group = infer_group_from_subject(subj)

    act_path = os.path.join(subj_path, "ACT", "connectome_ACT.csv")
    trekker_path = os.path.join(subj_path, "TREKKER", "connectome_TREKKER.csv")

    has_act = False
    has_trekker = False

    # Load ACT
    if os.path.exists(act_path):
        df_act = load_connectome_csv(act_path)
        act_connectomes[subj] = df_act
        has_act = True
        # Save a clean copy
        df_act.to_csv(os.path.join(output_dir, f"{subj}_ACT.csv"), index=False, header=False)

    # Load TREKKER
    if os.path.exists(trekker_path):
        df_trk = load_connectome_csv(trekker_path)
        trekker_connectomes[subj] = df_trk
        has_trekker = True
        # Save a clean copy
        df_trk.to_csv(os.path.join(output_dir, f"{subj}_TREKKER.csv"), index=False, header=False)

    subjects_meta.append({
        "subject": subj,
        "group": group,
        "has_ACT": has_act,
        "has_TREKKER": has_trekker
    })

# Save metadata table for quick overview
meta_df = pd.DataFrame(subjects_meta).sort_values("subject")
meta_path = os.path.join(output_dir, "subjects_metadata.csv")
meta_df.to_csv(meta_path, index=False)

print(f"Found ACT: {len(act_connectomes)} subjects | TREKKER: {len(trekker_connectomes)} subjects")
print(f"Saved per-subject CSVs to: {output_dir}")
print(f"Subjects metadata saved to: {meta_path}")


# =======================
# BUILD MERGED DATASETS
# =======================
long_rows = []  # long format
wide_rows = []  # wide format

def append_subject(subj: str, method: str, df: pd.DataFrame, group: str | None):
    vals, r_idx, c_idx = flatten_matrix(df, upper_only=use_upper_triangle_only, keep_diag=include_diagonal)

    # Long format rows
    for v, r, c in zip(vals, r_idx, c_idx):
        long_rows.append({
            "subject": subj,
            "method": method,
            "edge_row": int(r),
            "edge_col": int(c),
            "weight": float(v),
            "group": group
        })

    # Wide format row (one row per subject-method; columns are e_<r>_<c>)
    wide_entry = {
        "subject": subj,
        "method": method,
        "group": group
    }
    for v, r, c in zip(vals, r_idx, c_idx):
        wide_entry[f"e_{int(r)}_{int(c)}"] = float(v)
    wide_rows.append(wide_entry)

# Append ACT
for subj, df in act_connectomes.items():
    append_subject(subj, "ACT", df, infer_group_from_subject(subj))

# Append TREKKER
for subj, df in trekker_connectomes.items():
    append_subject(subj, "TREKKER", df, infer_group_from_subject(subj))

# Create DataFrames
df_long = pd.DataFrame(long_rows)
df_wide = pd.DataFrame(wide_rows)

# Save merged outputs
merged_long_path = os.path.join(output_dir, "connectomes_merged_long.csv")
merged_wide_path = os.path.join(output_dir, "connectomes_merged_wide.csv")
df_long.to_csv(merged_long_path, index=False)
df_wide.to_csv(merged_wide_path, index=False)

print(f"Long format saved to: {merged_long_path}")
print(f"Wide format saved to: {merged_wide_path}")

# =======================
# OPTIONAL: basic sanity prints
# =======================
if not df_wide.empty:
    n_edges = df_wide.shape[1] - 3  # subject, method, group
    print(f"Merged wide shape: {df_wide.shape} (edges={n_edges})")
if not df_long.empty:
    print(f"Merged long shape: {df_long.shape}")
    print(df_long.head())
