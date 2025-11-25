# OmicScope: Interactive Spatial Transcriptomics Explorer

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Dash](https://img.shields.io/badge/dash-2.0+-orange.svg)](https://dash.plotly.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**OmicScope** is a lightweight, interactive dashboard for exploring Spatial Transcriptomics data. It bridges the gap between complex molecular matrices and intuitive visual analysis, allowing researchers to explore gene expression in the context of tissue histology.

### 🚀 [Live Web Demo](https://omicscope.akulsaxena.me)
Explore a pre-processed **Human Lymph Node** dataset (10x Genomics Visium) directly in your browser without installing anything.

---

## ✨ Key Features

* **Dual View Mode:** Seamlessly toggle between **Spatial** (histological) and **UMAP** (transcriptional similarity) coordinates.
* **Histology Overlay:** View high-resolution H&E tissue images directly under your expression data.
* **Interactive Exploration:** Zoom, pan, and hover to inspect individual spots/cells.
* **Smart Gene Search:** Instantly visualize expression heatmaps for any gene in your dataset.
* **Live Metadata:** Automatically fetches gene descriptions and summaries via MyGene.info API.

---

## 📦 Local Installation & Usage

OmicScope is designed to be installed locally so you can analyze your own private datasets.

### 1. Installation

Clone the repository and install the package in editable mode:

```bash
git clone [https://github.com/YOUR_USERNAME/OmicScope.git](https://github.com/YOUR_USERNAME/OmicScope.git)
cd OmicScope
pip install -e .
````

### 2\. Running on Your Data

You can use OmicScope directly within a Python script or Jupyter Notebook. It accepts any standard `Scanpy` AnnData object.

```python
import scanpy as sc
import omicscope

# 1. Load your data
adata = sc.read_h5ad("my_tumor_sample.h5ad")

# 2. Launch the dashboard
omicscope.view(adata)
```

*The dashboard will start locally at `http://127.0.0.1:8050`.*

-----

## 🏗️ Architecture

OmicScope uses a "dual-mode" architecture to support both lightweight web deployment and robust local analysis:

| Feature | **CLI Mode** (Local) | **Web Mode** (Demo) |
| :--- | :--- | :--- |
| **Data Input** | Live `AnnData` object in memory | Pre-computed `.parquet` file |
| **Computation** | Real-time (Scanpy) | Pre-calculated (Pandas) |
| **Use Case** | Research & Analysis | Public Portfolio Demo |

-----

## 📁 Repository Structure

  * `omicscope/`: Core package source code.
      * `core.py`: Main Dash application logic and callbacks.
  * `web_demo/`: Scripts specifically for the Heroku deployment.
      * `app.py`: Entry point that loads the demo parquet file.
  * `scripts/`: Utilities for data preparation.
      * `preprocess.py`: Converts raw 10x datasets into lightweight parquet files for the web demo.

-----

## 🧬 Data Source

The demo dataset is the **Human Lymph Node** Spatial Gene Expression dataset provided by 10x Genomics.
