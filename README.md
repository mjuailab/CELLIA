# CELLIA: Curated tissue-specific knowledge and evidence Integration for LLM-based cell type Annotation

**CELLIA** is a tissue-aware LLM workflow for cell type annotation in scRNA-seq data. It integrates statistical evidence from differential expression with curated tissue-specific marker knowledge to generate structured cell type predictions.

---

## :rocket: Overview

Conventional cell type annotation relies on manual curation or reference mapping, which is often slow and inconsistent across tissues.
**CELLIA** introduces a tissue-aware, knowledge-driven workflow that:

1. Performs differentially expressed gene (DEG) analysis.
2. Filters candidate marker genes using curated tissue-specific marker databases.
3. Constructs structured prompts combining tissue context and top-k markers.
4. Queries an LLM (e.g., GPT, Claude, Gemini) for cell type annotation.
5. Collects model outputs—predicted label, rationale, and evidence score.
6. Stores results in AnnData/JSON/CSV format for reproducible analysis.
7. Visualizes results through an interactive web interface (UMAP viewer, dot plots, explanation panels).

---

## :file_folder: Repository Structure

```
CELLIA/
│
├── cellia.py               # Core workflow
├── cellia_web.py           # Web interface
│
├── run_cellia.py           # Script to run LLM-based annotation only
├── run_cellia_web.py       # Script to run full workflow (annotation + web)
├── cellia_tutorial.ipynb   # Tutorial to CELLIA workflow
│
├── cellia_output/          # CELLIA annotation results (JSON / CSV)
│
├── dataset/                # Example datasets (.h5ad)
│   └── CRC.h5ad            # Example AnnData file used for testing
│
├── database/               # Curated tissue-specific marker gene resources
│   └── Marker_DB.csv
│
├── requirements.txt        # Dependencies
└── README.md               # Documentation
```

---

## :gear: Installation

```bash
git clone https://github.com/ssjiyeong/CELLIA.git
cd CELLIA

# ---------------------------------------------------------
# Option A: If you have Conda installed (recommended)
# ---------------------------------------------------------
# Create and activate a Conda environment
conda create -n cellia_env python=3.11 -y
source ~/.bashrc     # If you need
conda activate cellia_env

# ---------------------------------------------------------
# Option B: If you do NOT have Conda installed
# ---------------------------------------------------------
# Create and activate a Python virtual environment instead
# (Use this only when Conda is unavailable)
python -m venv .venv
source .venv/bin/activate

# ---------------------------------------------------------
# Install dependencies
# ---------------------------------------------------------
pip install -r requirements.txt
```

> Recommended: Python ≥ 3.11 with Scanpy and Clustered AnnData

---

## ⭐️ Basic Usage

### ► Command Line Usage
After installing dependencies and creating an environment,
you can run CELLIA in three ways:
```bash
export API_KEY="YOUR_API_KEY"
```
#### I. Run full wrokflow (Step A-E)
```bash
python run_cellia_web.py \
  --adata dataset/YourAnnData.h5ad \
  --tissue_db "PBMC|blood" \
  --tissue_type "human PBMCs" \
  --n_top_markers 15 \
  --api_key "YOUR_API_KEY" \
  --model "gpt-4.1-2025-04-14" \
  --port 8060 \
  --rationale_json cellia_output/gpt_explanations_db.json
```
Then open the web interface in your brower:
```text
http://localhost:port
```

#### II. LLM-based annotation workflow (Step A-D)
**(a) Major cell type annotation**
```bash
python run_cellia.py \
  --adata dataset/YourAnnData.h5ad \
  --tissue_db "lung" \
  --tissue_type "human lung tissue" \
  --n_top_markers 15 \
  --api_key "YOUR_API_KEY" \
  --model "gemini-2.5-flash" \
  --mode "major"
```
**(b) Subtype-level cell type annotation**
```bash
python run_cellia.py \
  --adata dataset/YourAnnData.h5ad \
  --tissue_db "PBMC|blood" \
  --subset_db "dendritic|DC" \
  --tissue_type "human PBMC" \
  --n_top_markers 15 \
  --api_key "YOUR_API_KEY" \
  --model "gpt-4.1-2025-04-14" \
  --mode "subset" \
  --parent_celltype "Dendritic cells"
```

#### III. Interactive interface only
Runs the CELLIA web interface using pre-computed results. \
This requires that the input AnnData and {LLM}_explanations_db.json already follow the CELLIA annotation output format.
```bash
python run_cellia_web_only.py
```

### ► Python script / Jupyter notebooks Usage

```python
from cellia import *
import scanpy as sc

adata = sc.read_h5ad("dataset/CRC.h5ad")

adata = cellia_run(
    adata,
    tissue_db="crc|colon",
    tissue_type="human colorectal cancer (CRC)",
    api_key="YOUR_API_KEY,
    n_top_markers=15,
    model="gpt-4.1-2025-04-14"
)
```
---

## :blue_book: Tutorial

A Jupyter notebook tutorial is provided in **cellia_tutorial.ipynb**. \
It shows the full CELLIA workflow with example data.

---

## 💻 Output Example

```json
{
    "1": {
        "cell_type": "Liver sinusoidal endothelial cell (LSEC)",
        "marker_explanations": {
            "FCN3": "Strongly and selectively expressed in LSECs; important for immune functions in the liver vasculature.",
            "CLEC4G": "Highly specific marker of human LSECs, mediating scavenger and cell adhesion roles.",
            "CLEC4M": "Also known as L-SIGN, specifically marks LSECs in the liver.",
            "LYVE1": "Expressed in liver sinusoidal endothelial cells involved in endocytic and scavenger functions.",
            "FCN2": "Associated with LSEC immune surveillance functions."
        }
    },
}
```

---

## :jigsaw: Citation

If you use **CELLIA** in your work, please cite:

---

## :mailbox_with_mail: Contact

**Author:** Jiyeong Shin \
**Email:** sssjiyeong@gmail.com

---

## :scroll: License