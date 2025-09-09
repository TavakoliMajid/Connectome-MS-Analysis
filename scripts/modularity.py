# save as: compute_modularity_robust.py
import os, re, time
import numpy as np
import pandas as pd
import networkx as nx
import bct
from community import community_louvain  # package name: python-louvain

# =========================
# CONFIG
# =========================
BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes"  # <-- set your path
OUTDIR = os.path.join(BASE, "metrics")
os.makedirs(OUTDIR, exist_ok=True)

# Standardize sparsity across subjects/methods
TARGET_DENSITY = 0.08     # desired undirected density (0.05–0.10 typical). Will be capped to what exists.
KNN_K = 3                 # per-node minimal neighbors among existing non-zero edges (0 disables kNN backbone)

# Connectivity handling
CONNECTIVITY_MODE = "giant"  # "giant" -> analyze largest connected component only; "all" -> keep all components

# Modularity settings
GAMMA = 1.0                 # resolution parameter for modularity
FORCE_LOUVAIN_ONLY = True   # True -> skip BCT (quiet & robust on big matrices)

# File pattern (your per-subject copies)
FILE_RE = re.compile(r"^(INsIDER_[^_]+)_(ACT|TREKKER)\.csv$", re.IGNORECASE)

# =========================
# HELPERS
# =========================
def infer_group(subj: str):
    try:
        ch = subj.split("_", 1)[1][0].upper()
        return "control" if ch == "C" else ("patient" if ch == "P" else None)
    except Exception:
        return None

def load_W_csv(path: str) -> np.ndarray:
    """Strict numeric load -> float64 contiguous, symmetrize, zero diagonal, clamp negatives, rm NaN/Inf."""
    W = pd.read_csv(path, header=None, dtype=np.float64).values
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 0.0)
    W[W < 0] = 0.0
    W = np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)
    return np.ascontiguousarray(W, dtype=np.float64)

def maybe_normalize(W: np.ndarray, enable=True) -> np.ndarray:
    if not enable:
        return W
    m = float(np.max(W))
    return (W / m) if m > 0 else W

def density(W: np.ndarray) -> float:
    n = W.shape[0]
    if n < 2: return 0.0
    e = np.count_nonzero(np.triu(W, 1) > 0)
    return (2.0 * e) / (n * (n - 1))

def threshold_to_target_density(W: np.ndarray, target: float, knn_k: int = 0) -> tuple[np.ndarray, float, float]:
    """
    Keep edges so final undirected density ≈ target using ONLY existing non-zero weights.
    Returns (W_thr, target_used, available_density).
    """
    n = W.shape[0]
    if n < 2:
        return np.zeros_like(W), 0.0, 0.0

    r, c = np.triu_indices(n, 1)
    vals = W[r, c]
    pos = vals > 0
    pos_vals = vals[pos]
    pos_r, pos_c = r[pos], c[pos]
    m = pos_vals.size  # #available positive undirected edges
    total_possible = n * (n - 1) // 2
    avail_density = (m * 2.0) / (n * (n - 1)) if total_possible > 0 else 0.0

    if m == 0:
        return np.zeros_like(W), 0.0, avail_density

    # You cannot exceed what's available
    target_used = min(target, avail_density)
    if target_used <= 0:
        return np.zeros_like(W), target_used, avail_density

    k_edges = int(np.floor(target_used * total_possible))
    k_edges = max(1, min(k_edges, m))  # clamp to [1, m]

    # keep top k_edges positive weights
    kth = np.partition(pos_vals, -k_edges)[-k_edges]
    keep_mask = pos_vals >= kth

    M = np.zeros((n, n), dtype=bool)
    M[pos_r[keep_mask], pos_c[keep_mask]] = True
    M |= M.T
    np.fill_diagonal(M, False)

    # Optional k-NN backbone among existing non-zero neighbors
    if knn_k and knn_k > 0:
        for i in range(n):
            wi = W[i, :].copy()
            wi[i] = 0.0
            nz = np.where(wi > 0)[0]
            if nz.size == 0:
                continue
            k = min(knn_k, nz.size)
            topk = nz[np.argpartition(wi[nz], -k)[-k:]]
            M[i, topk] = True
            M[topk, i] = True

    U = np.where(M, W, 0.0)
    U = np.maximum(U, U.T)
    np.fill_diagonal(U, 0.0)
    return U, target_used, avail_density

def largest_component_submatrix(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (W_gcc, idx) where idx are node indices kept for the giant component."""
    n = W.shape[0]
    G = nx.from_numpy_array(W > 0.0)  # topology from >0 weights
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return W, np.arange(n)
    comps = list(nx.connected_components(G))
    if len(comps) <= 1:
        return W, np.arange(n)
    gcc = max(comps, key=len)
    idx = np.array(sorted(gcc), dtype=int)
    W_gcc = W[np.ix_(idx, idx)]
    return W_gcc, idx

def modularity_bct(W: np.ndarray, gamma: float) -> tuple[float,int,str]:
    Ci, Q = bct.modularity_und(W, gamma=gamma)
    return float(Q), int(len(np.unique(Ci))), "BCT"

def modularity_louvain(W: np.ndarray, gamma: float) -> tuple[float,int,str]:
    """Sparse weighted graph + Louvain at resolution=gamma."""
    n = W.shape[0]
    rows, cols = np.triu_indices(n, 1)
    weights = W[rows, cols]
    nz = weights > 0
    G = nx.Graph()
    G.add_nodes_from(range(n))
    if np.any(nz):
        edges = [(int(r), int(c), float(w)) for r, c, w in zip(rows[nz], cols[nz], weights[nz])]
        G.add_weighted_edges_from(edges, weight="weight")
    if G.number_of_edges() == 0:
        return float("nan"), 0, "LOUVAIN"
    part = community_louvain.best_partition(G, weight="weight", resolution=gamma)
    Q = community_louvain.modularity(part, G, weight="weight")
    n_comms = len(set(part.values()))
    return float(Q), int(n_comms), "LOUVAIN"

# =========================
# MAIN
# =========================
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

    # Load & normalize weights
    W = load_W_csv(fpath)
    W = maybe_normalize(W, True)

    # Sparsify to target density (capped by availability); optional kNN backbone
    W, target_used, avail = threshold_to_target_density(W, TARGET_DENSITY, knn_k=KNN_K)
    if target_used < TARGET_DENSITY:
        print(f"   [INFO] Target density {TARGET_DENSITY:.4f} capped to available {avail:.4f} (used {target_used:.4f}).")

    # Optionally restrict to largest connected component
    if CONNECTIVITY_MODE.lower() == "giant":
        W_gcc, kept_idx = largest_component_submatrix(W)
        if W_gcc.shape[0] < W.shape[0]:
            print(f"   [INFO] Using giant component: kept {W_gcc.shape[0]}/{W.shape[0]} nodes.")
        W = W_gcc

    n = W.shape[0]
    dens = density(W)

    # Modularity
    if FORCE_LOUVAIN_ONLY:
        Q, n_comms, which = modularity_louvain(W, GAMMA)
    else:
        try:
            Q, n_comms, which = modularity_bct(W, GAMMA)
        except Exception as e:
            print(f"   [INFO] BCT modularity failed ({e}); falling back to Louvain…")
            Q, n_comms, which = modularity_louvain(W, GAMMA)

    print(f"   -> [{which}] Q={Q:.6f} | communities={n_comms} | n={n} | dens={dens:.4f} | {time.time()-t0:.2f}s")

    rows.append({
        "subject": subj, "method": method, "group": group,
        "modularity_Q": Q, "n_communities": n_comms,
        "n_nodes": int(n), "density": float(dens), "file": fname,
        "algorithm": which, "used_density": float(target_used), "available_density": float(avail)
    })

out = pd.DataFrame(rows).sort_values(["method","group","subject"])
out_path = os.path.join(OUTDIR, "modularity_robust.csv")
out.to_csv(out_path, index=False)

print("="*64)
print(f"Saved: {out_path} | Total time: {time.time()-t_all:.2f}s")
print(out.head())
