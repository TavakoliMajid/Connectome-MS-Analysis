# save as: compute_characteristic_path_length_bct.py
import os, re, time
import numpy as np
import pandas as pd
import bct  # pip install bctpy

# ---------- CONFIG ----------
BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes"  # <-- set your path
OUTDIR = os.path.join(BASE, "metrics")
os.makedirs(OUTDIR, exist_ok=True)

# OPTIONAL SPEED-UP:
# Keep only the strongest p fraction of edges (same p for all subjects/methods).
# Example: 0.2 keeps top 20% strongest edges. Set to None to disable.
PROPORTIONAL_THRESHOLD = None  # e.g., 0.2 or None

# Files like: INsIDER_C019_ACT.csv, INsIDER_P39_TREKKER.csv
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

def load_and_clean_W(path: str) -> np.ndarray:
    """Load CSV as weights W (streamline counts), symmetrize, zero diag, clamp negatives to 0."""
    df = pd.read_csv(path, header=None)
    W = df.values.astype(np.float32)
    W = 0.5 * (W + W.T)      # symmetrize
    np.fill_diagonal(W, 0.0) # no self-loops
    W[W < 0] = 0.0
    return W

def maybe_threshold(W: np.ndarray, p: float | None) -> np.ndarray:
    """Optionally apply proportional thresholding to keep top p fraction of weights."""
    if p is None:
        return W
    return bct.threshold_proportional(W, p)

def characteristic_path_length_weighted(W: np.ndarray) -> float:
    """
    Weighted characteristic path length using bctpy.
      1) Convert weights to lengths L = 1/W (zero weights -> inf)
      2) D = distance_wei(L)
      3) charpath(D, include_diagonal=False)
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        L = np.where(W > 0, 1.0 / W, np.inf)

    D, _ = bct.distance_wei(L)

    # Compatible across bctpy versions (no 'nan' or 'disconnected' kw)
    lam, _, _, _, _ = bct.charpath(D, include_diagonal=False)
    return float(lam)

def matrix_density(W: np.ndarray) -> float:
    n = W.shape[0]
    if n < 2:
        return 0.0
    e = np.count_nonzero(np.triu(W, 1) > 0)  # undirected edges
    return (2.0 * e) / (n * (n - 1))

# ---------- MAIN ----------
files = sorted([f for f in os.listdir(BASE) if f.lower().endswith(".csv")])
targets = [f for f in files if FILE_RE.match(f)]
total = len(targets)

print(f"Scanning {BASE}")
print(f"Found {total} connectome CSVs matching pattern.")

rows = []
start_all = time.time()

for idx, fname in enumerate(targets, start=1):
    m = FILE_RE.match(fname)
    subj, method = m.group(1), m.group(2).upper()
    group = infer_group(subj)
    fpath = os.path.join(BASE, fname)

    t0 = time.time()
    print(f"[{idx}/{total}] {fname} | Subject={subj} | Method={method} | Group={group}")

    W = load_and_clean_W(fpath)
    if PROPORTIONAL_THRESHOLD is not None:
        W = maybe_threshold(W, PROPORTIONAL_THRESHOLD)

    n = W.shape[0]
    dens = matrix_density(W)

    try:
        cpl = characteristic_path_length_weighted(W)
        print(f"   -> CPL = {cpl:.6f} | n={n} | density={dens:.4f} | time={time.time()-t0:.2f}s")
    except Exception as e:
        cpl = np.nan
        print(f"   [WARN] CPL failed for {fname}: {e}")

    rows.append({
        "subject": subj,
        "method": method,
        "group": group,
        "char_path_length": cpl,
        "n_nodes": int(n),
        "density": float(dens),
        "file": fname
    })

out = pd.DataFrame(rows).sort_values(["method", "group", "subject"])
out_path = os.path.join(OUTDIR, "characteristic_path_length_bct.csv")
out.to_csv(out_path, index=False)

elapsed = time.time() - start_all
print("=" * 64)
print(f"Finished {len(rows)}/{total} files in {elapsed:.2f}s")
print(f"Saved: {out_path}")
print(out.head())
