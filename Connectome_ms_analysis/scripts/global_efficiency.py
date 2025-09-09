import os, re
import numpy as np
import pandas as pd
import bct  

# --------- PATHS ---------
BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes"  
OUTDIR = os.path.join(BASE, "metrics")
os.makedirs(OUTDIR, exist_ok=True)

# Match filenames like: INsIDER_C019_ACT.csv, INsIDER_P39_TREKKER.csv
FILE_RE = re.compile(r"^(INsIDER_[^_]+)_(ACT|TREKKER)\.csv$", re.IGNORECASE)

def infer_group(subj: str):
    try:
        tag = subj.split("_", 1)[1][0].upper()
        if tag == "C":
            return "control"
        if tag == "P":
            return "patient"
    except Exception:
        pass
    return None

def load_and_clean_matrix(path: str) -> np.ndarray:
    df = pd.read_csv(path, header=None)
    W = df.values.astype(np.float32)
    W = 0.5 * (W + W.T)      # symmetrize
    np.fill_diagonal(W, 0.0) # no self-loops
    W[W < 0] = 0.0           # no negative weights
    return W

rows = []
files = sorted([f for f in os.listdir(BASE) if f.lower().endswith(".csv")])

print(f"Found {len(files)} CSV files in {BASE}")
for idx, fname in enumerate(files, start=1):
    m = FILE_RE.match(fname)
    if not m:
        continue

    subj, method = m.group(1), m.group(2).upper()
    group = infer_group(subj)
    fpath = os.path.join(BASE, fname)

    print(f"[{idx}/{len(files)}] Processing {fname} (Subject={subj}, Method={method}, Group={group})...")

    W = load_and_clean_matrix(fpath)

    # --- Weighted Global Efficiency ---
    try:
        Eglob = float(bct.efficiency_wei(W))
        print(f"   -> Global efficiency = {Eglob:.4f}")
    except Exception as e:
        Eglob = np.nan
        print(f"   [WARN] efficiency failed for {fname}: {e}")

    n = W.shape[0]
    density = (np.count_nonzero(np.triu(W, 1) > 0) * 2) / (n * (n - 1)) if n > 1 else 0.0

    rows.append({
        "subject": subj,
        "method": method,
        "group": group,
        "global_efficiency": Eglob,
        "n_nodes": int(n),
        "density": float(density)
    })

out = pd.DataFrame(rows).sort_values(["method", "group", "subject"])
out_path = os.path.join(OUTDIR, "global_efficiency_bct.csv")
out.to_csv(out_path, index=False)

print("="*60)
print(f"Finished! Results saved to: {out_path}")
print(out.head())
