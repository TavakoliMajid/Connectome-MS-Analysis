# save as: compute_clustering_coefficient_bct.py
import os, re, time
import numpy as np
import pandas as pd
import bct  # pip install bctpy

BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes"  # <-- set your path
OUTDIR = os.path.join(BASE, "metrics")
os.makedirs(OUTDIR, exist_ok=True)

PROPORTIONAL_THRESHOLD = None   # e.g., 0.2 or None
NORMALIZE_WEIGHTS = True

FILE_RE = re.compile(r"^(INsIDER_[^_]+)_(ACT|TREKKER)\.csv$", re.IGNORECASE)

def infer_group(subj: str):
    try:
        ch = subj.split("_",1)[1][0].upper()
        return "control" if ch == "C" else ("patient" if ch == "P" else None)
    except Exception:
        return None

def load_W(path: str) -> np.ndarray:
    W = pd.read_csv(path, header=None).values.astype(np.float32)
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 0.0)
    W[W < 0] = 0.0
    return W

def maybe_threshold(W: np.ndarray, p):
    return bct.threshold_proportional(W, p) if p is not None else W

def maybe_normalize(W: np.ndarray, enable=True):
    if not enable: return W
    m = np.max(W)
    return W / m if m > 0 else W

def density(W: np.ndarray) -> float:
    n = W.shape[0]
    if n < 2: return 0.0
    e = np.count_nonzero(np.triu(W, 1) > 0)
    return (2.0 * e) / (n * (n - 1))

files = sorted([f for f in os.listdir(BASE) if FILE_RE.match(f)])
print(f"Found {len(files)} connectome CSVs.")

rows = []
t_all = time.time()
for i, fname in enumerate(files, 1):
    subj, method = FILE_RE.match(fname).groups()
    method = method.upper()
    group = infer_group(subj)
    fpath = os.path.join(BASE, fname)

    t0 = time.time()
    print(f"[{i}/{len(files)}] {fname} | {subj} | {method} | {group}")

    W = load_W(fpath)
    if PROPORTIONAL_THRESHOLD is not None:
        W = maybe_threshold(W, PROPORTIONAL_THRESHOLD)
    W = maybe_normalize(W, NORMALIZE_WEIGHTS)

    n = W.shape[0]; dens = density(W)
    try:
        C = bct.clustering_coef_wu(W)   # node-wise weighted clustering
        C_mean = float(np.nanmean(C))
        C_std  = float(np.nanstd(C))
        print(f"   -> meanC={C_mean:.6f} (±{C_std:.6f}) | n={n} | dens={dens:.4f} | {time.time()-t0:.2f}s")
    except Exception as e:
        C_mean = np.nan; C_std = np.nan
        print(f"   [WARN] Clustering failed: {e}")

    rows.append({
        "subject": subj, "method": method, "group": group,
        "mean_clustering": C_mean, "std_clustering": C_std,
        "n_nodes": int(n), "density": float(dens), "file": fname
    })

out = pd.DataFrame(rows).sort_values(["method","group","subject"])
out_path = os.path.join(OUTDIR, "clustering_coefficient_bct.csv")
out.to_csv(out_path, index=False)
print("="*64)
print(f"Saved: {out_path} | Total time: {time.time()-t_all:.2f}s")
print(out.head())
