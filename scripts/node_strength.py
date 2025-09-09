# save as: compute_node_strength.py
import os, re, time
import numpy as np
import pandas as pd

BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes"  # <-- set your path
OUTDIR = os.path.join(BASE, "metrics")
PER_NODE_DIR = os.path.join(OUTDIR, "node_metrics")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(PER_NODE_DIR, exist_ok=True)

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
    W = 0.5*(W + W.T)
    np.fill_diagonal(W, 0.0)
    W[W < 0] = 0.0
    return W

def threshold_proportional(W: np.ndarray, p):
    # simple proportional threshold (manual to avoid requiring bct here)
    if p is None: return W
    U = W.copy()
    tri = np.triu_indices_from(U, 1)
    vals = U[tri]
    k = int(np.floor(p * vals.size))
    if k <= 0: return np.zeros_like(W)
    thr = np.partition(vals, -k)[-k]
    mask = U >= thr
    U = np.where(mask, U, 0.0)
    U = np.maximum(U, U.T)
    np.fill_diagonal(U, 0.0)
    return U

def maybe_normalize(W: np.ndarray, enable=True):
    if not enable: return W
    m = np.max(W)
    return W / m if m > 0 else W

def density(W: np.ndarray) -> float:
    n = W.shape[0]
    if n < 2: return 0.0
    e = np.count_nonzero(np.triu(W,1) > 0)
    return (2.0 * e) / (n*(n-1))

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
    W = threshold_proportional(W, PROPORTIONAL_THRESHOLD)
    W = maybe_normalize(W, NORMALIZE_WEIGHTS)

    n = W.shape[0]; dens = density(W)
    # Node strength = sum of weights incident to each node
    strength = np.sum(W, axis=1).astype(np.float32)

    # Save per-node strengths (optional but handy)
    per_node_path = os.path.join(PER_NODE_DIR, f"{subj}_{method}_strength.csv")
    pd.DataFrame({"node": np.arange(n), "strength": strength}).to_csv(per_node_path, index=False)

    mean_s = float(np.mean(strength))
    med_s  = float(np.median(strength))
    std_s  = float(np.std(strength))

    print(f"   -> mean_strength={mean_s:.6f} | median={med_s:.6f} | std={std_s:.6f} | n={n} | dens={dens:.4f} | {time.time()-t0:.2f}s")

    rows.append({
        "subject": subj, "method": method, "group": group,
        "mean_strength": mean_s, "median_strength": med_s, "std_strength": std_s,
        "n_nodes": int(n), "density": float(dens), "file": fname,
        "per_node_file": per_node_path
    })

out = pd.DataFrame(rows).sort_values(["method","group","subject"])
out_path = os.path.join(OUTDIR, "node_strength.csv")
out.to_csv(out_path, index=False)
print("="*64)
print(f"Saved summary: {out_path}")
print(f"Per-node files in: {PER_NODE_DIR}")
print(out.head())
