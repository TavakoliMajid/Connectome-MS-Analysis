# save as: analyze_metrics.py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# try to import scipy for stats; if missing, we’ll degrade gracefully
try:
    from scipy import stats
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

# =========================
# CONFIG
# =========================
BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes\\metrics"  # <-- set your path
OUTFIG = os.path.join(BASE, "figures")
OUTSUM = os.path.join(BASE, "summaries")
os.makedirs(OUTFIG, exist_ok=True)
os.makedirs(OUTSUM, exist_ok=True)

METRIC_FILES = {
    "global_efficiency": "global_efficiency_bct.csv",
    "char_path_length": "characteristic_path_length_bct.csv",
    "modularity_Q": "modularity_robust.csv",
}

TITLES = {
    "global_efficiency": "Global Efficiency",
    "char_path_length": "Characteristic Path Length",
    "modularity_Q": "Modularity Q",
}

# =========================
# 1) Load metric CSVs
# =========================
ge = pd.read_csv(os.path.join(BASE, METRIC_FILES["global_efficiency"]))[
    ["subject", "method", "group", "global_efficiency"]
]
cpl = pd.read_csv(os.path.join(BASE, METRIC_FILES["char_path_length"]))[
    ["subject", "method", "group", "char_path_length"]
]
mod = pd.read_csv(os.path.join(BASE, METRIC_FILES["modularity_Q"]))[
    ["subject", "method", "group", "modularity_Q", "n_nodes", "density", "used_density", "available_density"]
]

# =========================
# 2) Merge into master table
# =========================
df = ge.merge(cpl, on=["subject", "method", "group"], how="outer")
df = df.merge(mod, on=["subject", "method", "group"], how="outer")

# Replace inf/-inf with NaN (so seaborn skips them)
df = df.replace([np.inf, -np.inf], np.nan)

# Optional: quick QC flag (useful for CPL/modularity reporting)
# Flag rows with very small giant component (heuristic using n_nodes column if present)
if "n_nodes" in df.columns:
    med_n = df["n_nodes"].dropna().median()
    df["qc_gcc_ok"] = df["n_nodes"] >= 0.90 * med_n
else:
    df["qc_gcc_ok"] = True

master_path = os.path.join(BASE, "metrics_master.csv")
df.to_csv(master_path, index=False)
print("Merged table preview:")
print(df.head())
print(f"\nSaved master table → {master_path}\n")

# =========================
# 3) Boxplots (patients vs controls, colored by method)
# =========================
sns.set_context("talk")
sns.set_style("whitegrid")

metrics = ["global_efficiency", "char_path_length", "modularity_Q"]

for metric in metrics:
    if df[metric].dropna().empty:
        print(f"Skipping boxplot for {metric} due to no valid data")
        continue
    plt.figure(figsize=(8, 5))
    ax = sns.boxplot(data=df, x="group", y=metric, hue="method", palette="Set2", showfliers=False)
    # overlay raw points (no hue to avoid seaborn legend bug)
    sns.stripplot(data=df, x="group", y=metric, dodge=True, alpha=0.45, color="black", size=3)
    ax.set_title(f"{TITLES[metric]} by Group & Method")
    ax.set_xlabel("")
    ax.set_ylabel(TITLES[metric])
    ax.legend(title="Method")
    plt.tight_layout()
    outpath = os.path.join(OUTFIG, f"boxplot_{metric}.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved {outpath}")

# =========================
# 4) Scatter plots (ACT vs TREKKER, per subject)
# =========================
for metric in metrics:
    piv = df.pivot_table(index=["subject", "group"], columns="method", values=metric, aggfunc="mean")
    piv = piv.reset_index()
    # Check if pivot table has valid data and required columns
    if piv.empty or not all(col in piv.columns for col in ["ACT", "TREKKER"]):
        print(f"Skipping scatter for {metric} due to missing or invalid paired data (columns: {list(piv.columns)})")
        continue
    # Keep rows that have both ACT & TREKKER
    piv = piv.dropna(subset=["ACT", "TREKKER"], how="any")
    if piv.empty:
        print(f"Skipping scatter for {metric} due to no paired data after dropna")
        continue
    plt.figure(figsize=(6, 6))
    ax = sns.scatterplot(data=piv, x="ACT", y="TREKKER", hue="group", s=60)
    # diagonal
    vmin = np.nanmin([piv["ACT"].min(), piv["TREKKER"].min()])
    vmax = np.nanmax([piv["ACT"].max(), piv["TREKKER"].max()])
    pad = 0.05 * (vmax - vmin) if np.isfinite(vmax - vmin) else 1.0
    lims = [vmin - pad, vmax + pad]
    ax.plot(lims, lims, 'k--', alpha=0.6)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f"ACT vs TREKKER ({TITLES[metric]})")
    ax.set_xlabel("ACT"); ax.set_ylabel("TREKKER")
    plt.tight_layout()
    outpath = os.path.join(OUTFIG, f"scatter_{metric}.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved {outpath}")

# =========================
# 5) Group summaries & stats
# =========================
def summarize_by_group(df_in: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Mean ± SD by (group, method) for a metric."""
    g = (
        df_in.groupby(["group", "method"], dropna=False)[metric]
        .agg(["count", "mean", "std", "median"])
        .reset_index()
    )
    g["metric"] = metric
    return g[["metric", "group", "method", "count", "mean", "std", "median"]]

def mannwhitney_and_ttest_ind(df_in: pd.DataFrame, metric: str, method: str):
    """Patients vs Controls (independent) within a method."""
    a = df_in[(df_in["method"] == method) & (df_in["group"] == "patient")][metric].dropna()
    b = df_in[(df_in["method"] == method) & (df_in["group"] == "control")][metric].dropna()
    res = {"metric": metric, "method": method, "n_patient": len(a), "n_control": len(b)}
    if _HAS_SCIPY and len(a) > 1 and len(b) > 1:
        try:
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            res["mwU"] = u.statistic; res["mwP"] = u.pvalue
        except Exception:
            res["mwU"] = np.nan; res["mwP"] = np.nan
        try:
            t = stats.ttest_ind(a, b, equal_var=False)
            res["t"] = t.statistic; res["tp"] = t.pvalue
        except Exception:
            res["t"] = np.nan; res["tp"] = np.nan
    else:
        res["mwU"] = np.nan; res["mwP"] = np.nan; res["t"] = np.nan; res["tp"] = np.nan
    return res

def wilcoxon_and_paired_t(df_in: pd.DataFrame, metric: str):
    """ACT vs TREKKER (paired within subject), pooled across groups (and also per-group)."""
    records = []
    # pooled
    piv = df_in.pivot_table(index=["subject"], columns="method", values=metric, aggfunc="mean")
    pair = piv.dropna(subset=["ACT", "TREKKER"], how="any")
    res = {"metric": metric, "group": "ALL", "n_pairs": len(pair)}
    if _HAS_SCIPY and len(pair) > 1:
        try:
            w = stats.wilcoxon(pair["ACT"], pair["TREKKER"])
            res["wilcoxonW"] = w.statistic; res["wilcoxonP"] = w.pvalue
        except Exception:
            res["wilcoxonW"] = np.nan; res["wilcoxonP"] = np.nan
        try:
            t = stats.ttest_rel(pair["ACT"], pair["TREKKER"])
            res["t"] = t.statistic; res["tp"] = t.pvalue
        except Exception:
            res["t"] = np.nan; res["tp"] = np.nan
    else:
        res["wilcoxonW"] = np.nan; res["wilcoxonP"] = np.nan; res["t"] = np.nan; res["tp"] = np.nan
    records.append(res)

    # per-group
    for grp in ["patient", "control"]:
        pivg = df_in[df_in["group"] == grp].pivot_table(index=["subject"], columns="method", values=metric, aggfunc="mean")
        pairg = pivg.dropna(subset=["ACT", "TREKKER"], how="any")
        resg = {"metric": metric, "group": grp, "n_pairs": len(pairg)}
        if _HAS_SCIPY and len(pairg) > 1:
            try:
                w = stats.wilcoxon(pairg["ACT"], pairg["TREKKER"])
                resg["wilcoxonW"] = w.statistic; resg["wilcoxonP"] = w.pvalue
            except Exception:
                resg["wilcoxonW"] = np.nan; resg["wilcoxonP"] = np.nan
            try:
                t = stats.ttest_rel(pairg["ACT"], pairg["TREKKER"])
                resg["t"] = t.statistic; resg["tp"] = t.pvalue
            except Exception:
                resg["t"] = np.nan; resg["tp"] = np.nan
        else:
            resg["wilcoxonW"] = np.nan; resg["wilcoxonP"] = np.nan; resg["t"] = np.nan; resg["tp"] = np.nan
        records.append(resg)
    return pd.DataFrame.from_records(records)

# Prepare outputs
all_summaries = []
all_pv_pc = []
all_pv_acttrk = []

for metric in metrics:
    # Summary table
    summ = summarize_by_group(df, metric)
    all_summaries.append(summ)

    # Patients vs Controls (per method)
    for mth in ["ACT", "TREKKER"]:
        all_pv_pc.append(mannwhitney_and_ttest_ind(df, metric, mth))

    # ACT vs TREKKER (paired)
    all_pv_acttrk.append(wilcoxon_and_paired_t(df, metric))

summary_df = pd.concat(all_summaries, ignore_index=True)
pc_tests_df = pd.DataFrame(all_pv_pc)
acttrk_tests_df = pd.concat(all_pv_acttrk, ignore_index=True)

summary_df.to_csv(os.path.join(OUTSUM, "group_summaries.csv"), index=False)
pc_tests_df.to_csv(os.path.join(OUTSUM, "patients_vs_controls_tests.csv"), index=False)
acttrk_tests_df.to_csv(os.path.join(OUTSUM, "act_vs_trekker_paired_tests.csv"), index=False)

print("\n=== Group summaries (mean ± SD) ===")
print(summary_df)

print("\n=== Patients vs Controls (per method) ===")
print(pc_tests_df)

print("\n=== ACT vs TREKKER (paired) ===")
print(acttrk_tests_df)

print(f"\nSaved summaries to:\n  {os.path.join(OUTSUM, 'group_summaries.csv')}\n  {os.path.join(OUTSUM, 'patients_vs_controls_tests.csv')}\n  {os.path.join(OUTSUM, 'act_vs_trekker_paired_tests.csv')}")