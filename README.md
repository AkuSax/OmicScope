# OmicScope: Interactive Spatial Transcriptomics Explorer

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Dash](https://img.shields.io/badge/dash-2.0+-orange.svg)](https://dash.plotly.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**OmicScope** is a lightweight, interactive dashboard for exploring Spatial Transcriptomics data. It bridges the gap between complex molecular matrices and intuitive visual analysis, allowing researchers to explore gene expression in the context of tissue histology.

<p align="center">
  <img src="assets/demo.gif" alt="OmicScope Demo" width="800">
</p>

### 🚀 [Try the Live Demo](https://omicscope.akulsaxena.me)
*See the tool in action on a human lymph node dataset.*

---

## 📦 Installation

OmicScope is designed to be installed locally, allowing you to analyze your own private **Scanpy** datasets.

Clone the repository and install it in editable mode:

```bash
git clone [https://github.com/AkuSax/OmicScope.git](https://github.com/AkuSax/OmicScope.git)
cd OmicScope
pip install -e .
````

-----

## 💻 Usage

You can use OmicScope directly within a Python script or Jupyter Notebook. It accepts any standard `AnnData` object (the standard format for single-cell and spatial data).

```python
import scanpy as sc
import omicscope

# 1. Load your dataset (e.g., 10x Genomics Visium)
# Ensure your object has 'spatial' in .obsm and 'leiden' in .obs
adata = sc.read_h5ad("path/to/my_sample.h5ad")

# 2. Launch the dashboard
omicscope.view(adata)
```

The dashboard will launch automatically at `http://127.0.0.1:8050`.

-----

## ✨ Key Features

  * **Dual View Mode:** Seamlessly toggle between **Spatial** (histological) and **UMAP** (transcriptional similarity) coordinates.
  * **Histology Overlay:** View high-resolution H\&E tissue images directly under your expression data.
  * **Interactive Exploration:** Zoom, pan, and hover to inspect individual spots/cells.
  * **Smart Gene Search:** Instantly visualize expression heatmaps for any gene in your dataset.
  * **Live Metadata:** Automatically fetches gene descriptions and summaries via the MyGene.info API.

-----

## 🧬 Data Requirements

OmicScope expects a standard `Scanpy` AnnData object with the following:

  * **`.obsm['spatial']`**: Spatial coordinates (for the tissue view).
  * **`.obsm['X_umap']`**: UMAP coordinates (for the cluster view).
  * **`.obs['leiden']`**: Cluster labels.
