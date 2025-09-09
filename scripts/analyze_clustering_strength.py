# save as: analyze_clustering_strength.py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =========================
# CONFIG
# =========================
BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes\\metrics"  # <-- set your path
OUTFIG = os.path.join(BASE, "figures")
OUTSUM = os.path.join(BASE, "summaries")
os.makedirs(OUTFIG, exist_ok=True)
os.makedirs(OUTSUM, exist_ok=True)

CLUSTERING_FILE = "clustering_coefficient_bct.csv"
STRENGTH_FILE   = "node_strength.csv"

TITLES = {
    "mean_clustering": "Mean Clustering Coefficient",
    "std_clustering": "Clustering Coefficient (Std)",
    "mean_strength": "Mean Node Strength",
    "median_strength": "Median Node Strength",
    "std_strength": "Node Strength (Std)",
}

# =========================
# 1) Load CSVs
# =========================
clu_path = os.path.join(BASE, CLUSTERING_FILE)
str_path = os.path.join(BASE, STRENGTH_FILE)

if not os.path.exists(clu_path):
    raise FileNotFoundError(f"Missing file: {clu_path}")
if not os.path.exists(str_path):
    raise FileNotFoundError(f"Missing file: {str_path}")

clu = pd.read_csv(clu_path)
strn = pd.read_csv(str_path)

# Expected columns (we’ll keep what exists)
keep_clu = ["subject","method","group","mean_clustering","std_clustering","n_nodes","density","file"]
keep_str = ["subject","method","group","mean_strength","median_strength","std_strength","n_nodes","density","file"]
clu = clu[[c for c in keep_clu if c in clu.columns]].copy()
strn = strn[[c for c in keep_str if c in strn.columns]].copy()

# =========================
# 2) Merge → master table
# =========================
df = pd.merge(clu, strn, on=["subject","method","group"], how="outer", suffixes=("_clu","_str"))
df = df.replace([np.inf, -np.inf], np.nan)

master_path = os.path.join(BASE, "metrics_clustering_strength_master.csv")
df.to_csv(master_path, index=False)
print("Merged table preview:")
print(df.head())
print(f"\nSaved master table → {master_path}\n")

# =========================
# 3) Plots
# =========================
sns.set_context("talk")
sns.set_style("whitegrid")

def _safe_boxplot(df_in, ycol: str, title: str, fname: str):
    vals = df_in[ycol].dropna() if ycol in df_in.columns else pd.Series(dtype=float)
    if vals.empty:
        print(f"Skipping boxplot for {ycol}: no valid data.")
        return
    plt.figure(figsize=(8,5))
    ax = sns.boxplot(data=df_in, x="group", y=ycol, hue="method", palette="Set2", showfliers=False)
    sns.stripplot(data=df_in, x="group", y=ycol, dodge=True, alpha=0.45, color="black", size=3)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(title)
    ax.legend(title="Method")
    plt.tight_layout()
    outpath = os.path.join(OUTFIG, fname)
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved {outpath}")

def _safe_scatter_act_trekker(df_in, valcol: str, title: str, fname: str):
    if valcol not in df_in.columns:
        print(f"Skipping scatter for {valcol}: column missing.")
        return
    piv = df_in.pivot_table(index=["subject","group"], columns="method", values=valcol, aggfunc="mean").reset_index()
    if not {"ACT","TREKKER"}.issubset(set(piv.columns)):
        print(f"Skipping scatter for {valcol}: ACT/TREKKER columns missing.")
        return
    piv = piv.dropna(subset=["ACT","TREKKER"], how="any")
    if piv.empty:
        print(f"Skipping scatter for {valcol}: no paired rows after dropna.")
        return
    plt.figure(figsize=(6,6))
    ax = sns.scatterplot(data=piv, x="ACT", y="TREKKER", hue="group", s=60)
    vmin = np.nanmin([piv["ACT"].min(), piv["TREKKER"].min()])
    vmax = np.nanmax([piv["ACT"].max(), piv["TREKKER"].max()])
    pad = 0.05 * (vmax - vmin) if np.isfinite(vmax - vmin) else 1.0
    lims = [vmin - pad, vmax + pad]
    ax.plot(lims, lims, "k--", alpha=0.6)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    ax.set_xlabel("ACT"); ax.set_ylabel("TREKKER")
    plt.tight_layout()
    outpath = os.path.join(OUTFIG, fname)
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved {outpath}")

# Boxplots
_safe_boxplot(df, "mean_clustering", "Mean Clustering Coefficient", "boxplot_mean_clustering.png")
_safe_boxplot(df, "mean_strength",   "Mean Node Strength",        "boxplot_mean_strength.png")
# Optional:
_safe_boxplot(df, "median_strength", "Median Node Strength",      "boxplot_median_strength.png")

# ACT vs TREKKER scatters (optional but useful)
_safe_scatter_act_trekker(df, "mean_clustering", "ACT vs TREKKER (Mean Clustering)", "scatter_mean_clustering.png")
_safe_scatter_act_trekker(df, "mean_strength",   "ACT vs TREKKER (Mean Node Strength)", "scatter_mean_strength.png")

# =========================
# 4) Descriptive summaries
# =========================
def summarize(df_in: pd.DataFrame, col: str) -> pd.DataFrame:
    if col not in df_in.columns:
        return pd.DataFrame()
    g = (
        df_in.groupby(["group","method"], dropna=False)[col]
        .agg(["count","mean","std","median"])
        .reset_index()
    )
    g.insert(0, "metric", col)
    return g

sum_frames = []
for mcol in ["mean_clustering","mean_strength","median_strength"]:
    s = summarize(df, mcol)
    if not s.empty:
        sum_frames.append(s)

summary_df = pd.concat(sum_frames, ignore_index=True) if sum_frames else pd.DataFrame()
if not summary_df.empty:
    outsum = os.path.join(OUTSUM, "group_summaries_clustering_strength.csv")
    summary_df.to_csv(outsum, index=False)
    print("\n=== Group summaries (mean ± SD) ===")
    print(summary_df)
    print(f"\nSaved summaries → {outsum}")
else:
    print("No summary tables written (no valid columns).")
