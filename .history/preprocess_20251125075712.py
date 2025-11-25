import scanpy as sc
import pandas as pd
import numpy as np
import os
from tqdm import tqdm  # Import the progress bar library
import warnings

# Suppress minor warnings for cleaner output
warnings.filterwarnings("ignore")

def preprocess_and_export(output_path="web_demo/data/demo_data.parquet"):
    """
    Loads a standard Visium dataset, processes it, finds markers,
    and exports a lightweight Parquet file for the web demo.
    """
    # 1. Set Scanpy verbosity to 3 (shows hints/progress for internal steps)
    sc.settings.verbosity = 3
    
    print("--- STEP 1: LOADING DATA ---")
    print("Downloading/Loading Visium dataset (V1_Human_Lymph_Node)...")
    print("(This may take 1-2 minutes if downloading for the first time)")
    
    # Check if data is already cached to manage expectations
    if os.path.exists("./data/V1_Human_Lymph_Node"):
        print("Found local cache.")
    
    # Load dataset
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
    adata.var_names_make_unique()
    print(f"Loaded successfully: {adata.shape[0]} cells, {adata.shape[1]} genes.")

    print("\n--- STEP 2: PREPROCESSING ---")
    # Basic filtering
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("\n--- STEP 3: CLUSTERING & EMBEDDING ---")
    # These steps are computationally intensive
    sc.pp.hvg(adata, n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    print("\n--- STEP 4: IDENTIFYING MARKERS ---")
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    
    # Extract top 5 genes per cluster
    top_genes = set()
    num_clusters = len(adata.obs['leiden'].unique())
    
    # Use TQDM to show progress for this loop
    print("Extracting top marker genes...")
    for i in tqdm(range(num_clusters), desc="Processing Clusters"):
        try:
            # Get top 5 genes for this cluster
            genes = sc.get.rank_genes_groups_df(adata, group=str(i)).head(5)['names'].tolist()
            top_genes.update(genes)
        except Exception as e:
            print(f"Warning processing cluster {i}: {e}")
    
    print(f"Identified {len(top_genes)} unique marker genes to preserve.")

    print("\n--- STEP 5: EXPORTING LITE DATA ---")
    # Initialize DataFrame with cell indices
    df = pd.DataFrame(index=adata.obs.index)

    # Add Coordinates
    # Visium data usually has 'spatial', but we add a fallback just in case
    if 'spatial' in adata.obsm:
        df['spatial_x'] = adata.obsm['spatial'][:, 0]
        df['spatial_y'] = adata.obsm['spatial'][:, 1]
    else:
        df['spatial_x'] = adata.obsm['X_umap'][:, 0]
        df['spatial_y'] = adata.obsm['X_umap'][:, 1]

    df['umap_1'] = adata.obsm['X_umap'][:, 0]
    df['umap_2'] = adata.obsm['X_umap'][:, 1]

    # Add Metadata
    df['cluster'] = adata.obs['leiden'].astype(str)

    # Add Gene Expression Counts with progress bar
    print("Flattening gene expression matrix...")
    for gene in tqdm(top_genes, desc="Exporting Genes"):
        if gene in adata.var_names:
            try:
                # Handle both sparse and dense matrices safely
                if pd.api.types.is_sparse(adata[:, gene].X):
                     data = adata[:, gene].X.toarray().flatten()
                elif hasattr(adata[:, gene].X, "toarray"):
                     data = adata[:, gene].X.toarray().flatten()
                else:
                     data = adata[:, gene].X.flatten()
                
                df[gene] = data
            except Exception as e:
                pass

    # Ensure directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save to parquet
    df.to_parquet(output_path, compression='snappy')
    
    file_size = os.path.getsize(output_path) / (1024*1024)
    print(f"\nSUCCESS: Data exported to {output_path}")
    print(f"Final File Size: {file_size:.2f} MB")

if __name__ == "__main__":
    preprocess_and_export()