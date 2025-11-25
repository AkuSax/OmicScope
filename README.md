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
git clone [https://github.com/AkuSax/OmicScope.git](https://github.com/AkuSax/OmicScope.git)
cd OmicScope
pip install -e .
