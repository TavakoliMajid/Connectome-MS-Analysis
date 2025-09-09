# Connectome Analysis in Multiple Sclerosis (ACT vs TREKKER)

This project compares **Anatomically Constrained Tractography (ACT)** and **TREKKER** tractography pipelines for constructing structural brain connectomes in **Multiple Sclerosis (MS) patients vs healthy controls**.  
The pipeline reconstructs **≥1M streamlines per subject**, builds connectomes, and computes graph-theoretic metrics to assess group and method differences.



---

## 📂 Repository structure

```
connectome-ms-analysis/
│
├── bash/                         # tractography pipelines
│   ├── ACT_pipeline.sh            # ACT pipeline (MRtrix3)
│   └── TREKKER_pipeline.sh        # TREKKER pipeline
│
├── scripts/                       # Python metric analysis
│   ├── load_csv.py
│   ├── inspect_connectome.py
│   ├── global_efficiency.py
│   ├── characteristic_path_lengh.py
│   ├── modularity.py
│   ├── clustering_coefficient.py
│   ├── node_strength.py
│   ├── analyze_metrics.py
│   └── analyze_clustering_strength.py
│
├── results/                       # example outputs (not tracked here)
│   ├── figures/                   # PNG plots
│   ├── summaries/                 # CSV summary stats
│   └── metrics_master.csv
│
├
│
└── README.md
```

---

## 🧪 Pipelines

### ACT (MRtrix3)
Script: `bash/ACT_pipeline.sh`  
Steps:
- Convert DWI → `.mif`
- Estimate tissue responses (`dwi2response`)
- MSMT CSD → FODs
- Build 5TT + GMWMI mask
- **Tractography**: `tckgen -act -seed_gmwmi -select 1M`
- `tck2connectome` with `-assignment_radial_search 2.5 -zero_diagonal -symmetric`
- Output: `INsIDER_<ID>_ACT.csv`

**Installation**: Follow instructions on the official MRtrix3 ACT documentation:  
👉 [MRtrix3 ACT Installation Guide](https://github.com/MRtrix3/mrtrix3/blob/master/docs/quantitative_structural_connectivity/act.rst)

---

### TREKKER
Script: `bash/TREKKER_pipeline.sh`  
Steps:
- Convert DWI + build FODs
- Generate 5TT + GMWMI
- **Trekker tractography**: `trekker track ... --seed GMWMI`
- Convert `.vtk → .tck`
- `tck2connectome` with `-assignment_radial_search 2.5 -zero_diagonal -symmetric`
- Output: `INsIDER_<ID>_TREKKER.csv`

**Installation**: Follow instructions on the TREKKER GitHub page:  
👉 [TREKKER GitHub](https://github.com/dmritrekker/trekker)

---

Both pipelines use **GMWMI seeding**, identical assignment radius, and **streamline counts** for fair comparison.

---

## 📊 Metrics computed

Each script saves a CSV in `Processed_Connectomes/metrics/`:

- **Global Efficiency** → `global_efficiency_bct.csv`  
- **Characteristic Path Length (CPL)** → `characteristic_path_length_bct.csv`  
- **Modularity Q** → `modularity_robust.csv`  
- **Clustering Coefficient** → `clustering_coefficient_bct.csv`  
- **Node Strength** → `node_strength.csv`

Utility scripts:
- `inspect_connectome.py` → checks CSVs are square/valid  
- `load_csv.py` → merges ACT/TREKKER into tables  

---

 📈 Analysis & Visualization

 `analyze_metrics.py`
- Merges global_efficiency, CPL, modularity into one master table
- Cleans NaN/Inf values
- Generates:
  - **Boxplots** (patients vs controls, ACT vs TREKKER)
  - **Scatter plots** (ACT vs TREKKER per subject)
- Outputs:
  - `metrics_master.csv`
  - Group summaries (`group_summaries.csv`)
  - Stats (Mann–Whitney, t-tests, Wilcoxon)

 `analyze_clustering_strength.py`
- Specialized plots for clustering coefficient & node strength
- Boxplots & scatterplots
- Summaries → `group_summaries_clustering_strength.csv`

---

## ▶️ How to run

1. **Precompute connectomes**  
   Run ACT or TREKKER pipelines:
   
   bash/ACT_pipeline.sh
   bash/TREKKER_pipeline.sh
   

2. **Check and load connectomes**  
   
   python scripts/inspect_connectome.py
   python scripts/load_csv.py
   

3. **Compute metrics**  
   
   python scripts/global_efficiency.py
   python scripts/characteristic_path_lengh.py
   python scripts/modularity.py
   python scripts/clustering_coefficient.py
   python scripts/node_strength.py
   

4. **Aggregate & visualize**  
   
   python scripts/analyze_metrics.py
   python scripts/analyze_clustering_strength.py
   

Figures are saved in `metrics/figures/`, summaries in `metrics/summaries/`.

---

## 🖥️ System & Environment

- **ACT & TREKKER analysis**: Ubuntu 24.04 (WSL2)  
- **Python analysis**: Windows 11  

Hardware:  
- CPU: Intel i5-14600K (14 cores), RAM 32 GB DDR5  
- GPU: NVIDIA RTX 3050, CUDA 12.7  
- Storage: 1 TB NVMe SSD  

Software:  
- MRtrix3 (v3.0.4-git)  
- TREKKER (v1.0.0-rc5)  
- Python 3.10 (Conda env `image`)  

### Python deps (requirements.txt)
```
numpy
pandas
matplotlib
seaborn
networkx
bctpy
python-louvain
scipy
```

---

## 📚 References

- [MRtrix3: Anatomically Constrained Tractography](https://github.com/MRtrix3/mrtrix3/blob/master/docs/quantitative_structural_connectivity/act.rst)  
- [Trekker GitHub](https://github.com/dmritrekker/trekker)  
- [BCTpy: Brain Connectivity Toolbox](https://github.com/aestrivex/bctpy)   

---

## ✨ Author
Developed by **Majid Tavakoli (2025)**  *“ACT vs TREKKER connectome analysis in MS”*.
