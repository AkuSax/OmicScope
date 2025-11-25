import scanpy as sc
import pandas as pd
import numpy as np
import os

def preprocess_and_export(output_path="web_demo/data/demo_data.parquet"):
    """
    Loads a standard Visium dataset, processes it, finds markers,
    and exports a lightweight Parquet file for the web demo.
    """
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
    
    # Make variable names unique to prevent errors during export
    adata.var_names_make_unique()

    # Preprocessing
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Clustering
    print("Running PCA, Neighbors, UMAP, and Leiden...")
    sc.pp.hvg(adata, n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # Identify Marker Genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    
    # Extract top 5 genes per cluster
    top_genes = set()
    num_clusters = len(adata.obs['leiden'].unique())
    
    for i in range(num_clusters):
        genes = sc.get.rank_genes_groups_df(adata, group=str(i)).head(5)['names'].tolist()
        top_genes.update(genes)
    
    print(f"Identified {len(top_genes)} unique marker genes to preserve.")

    df = pd.DataFrame(index=adata.obs.index)
    # Check if spatial info exists, otherwise mock it (fallback safety)
    if 'spatial' in adata.obsm:
        df['spatial_x'] = adata.obsm['spatial'][:, 0]
        df['spatial_y'] = adata.obsm['spatial'][:, 1]
    else:
        print("Warning: No spatial info found. Using UMAP as spatial placeholder.")
        df['spatial_x'] = adata.obsm['X_umap'][:, 0]
        df['spatial_y'] = adata.obsm['X_umap'][:, 1]

    df['umap_1'] = adata.obsm['X_umap'][:, 0]
    df['umap_2'] = adata.obsm['X_umap'][:, 1]

    # Add Metadata (Clusters)
    df['cluster'] = adata.obs['leiden'].astype(str)

    # Add Gene Expression Counts (Only for the top genes)
    # We extract the normalized matrix for the specific genes
    # Check if sparse or dense
    for gene in top_genes:
        if gene in adata.var_names:
            try:
                # Handle sparse matrix extraction
                if pd.api.types.is_sparse(adata[:, gene].X):
                     data = adata[:, gene].X.toarray().flatten()
                elif hasattr(adata[:, gene].X, "toarray"):
                     data = adata[:, gene].X.toarray().flatten()
                else:
                     data = adata[:, gene].X.flatten()
                
                df[gene] = data
            except Exception as e:
                print(f"Could not export gene {gene}: {e}")

    # Export
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_parquet(output_path, compression='snappy')
    
    print("Done! File size:")
    print(f"{os.path.getsize(output_path) / (1024*1024):.2f} MB")

if __name__ == "__main__":
    preprocess_and_export()