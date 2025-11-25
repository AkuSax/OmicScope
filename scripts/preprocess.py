import scanpy as sc
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import warnings
from PIL import Image

warnings.filterwarnings("ignore")

def preprocess_and_export(output_dir="web_demo/data"):
    """
    Exports data and the histology image. 
    Scales the spot coordinates to match the image pixels.
    """
    sc.settings.verbosity = 3
    
    # Load dataset
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
    adata.var_names_make_unique()
    lib_id = list(adata.uns['spatial'].keys())[0]
    
    hires_image = adata.uns['spatial'][lib_id]['images']['hires']
    scale_factor = adata.uns['spatial'][lib_id]['scalefactors']['tissue_hires_scalef']
    
    print(f"Scale factor found: {scale_factor}")

    # Save Image
    # Convert numpy array (0-1) to integer (0-255) for PIL
    im_data = (hires_image * 255).astype(np.uint8)
    img = Image.fromarray(im_data)
    
    os.makedirs(output_dir, exist_ok=True)
    image_path = os.path.join(output_dir, "tissue.png")
    img.save(image_path)
    print(f"Saved histology image to {image_path}")

    print("\n--- STEP 3: PREPROCESSING ---")
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print("Running Clustering & Embedding...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    print("\n--- STEP 4: IDENTIFYING MARKERS ---")
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    
    # Extract top 5 genes per cluster
    top_genes = set()
    num_clusters = len(adata.obs['leiden'].unique())
    
    for i in tqdm(range(num_clusters), desc="Finding Markers"):
        try:
            genes = sc.get.rank_genes_groups_df(adata, group=str(i)).head(5)['names'].tolist()
            top_genes.update(genes)
        except Exception as e:
            print(f"Warning cluster {i}: {e}")

    print(f"Identified {len(top_genes)} unique marker genes.")

    print("\n--- STEP 5: EXPORTING DATA ---")
    # Initialize DataFrame with cell indices
    df = pd.DataFrame(index=adata.obs.index)

    # CRITICAL: Scale spatial coordinates to match the image
    df['spatial_x'] = adata.obsm['spatial'][:, 0] * scale_factor
    df['spatial_y'] = adata.obsm['spatial'][:, 1] * scale_factor

    df['umap_1'] = adata.obsm['X_umap'][:, 0]
    df['umap_2'] = adata.obsm['X_umap'][:, 1]
    df['cluster'] = adata.obs['leiden'].astype(str)

    # Robust Export Loop
    print("Flattening gene expression matrix...")
    for gene in tqdm(top_genes, desc="Exporting Genes"):
        if gene in adata.var_names:
            # Extract the column for the specific gene
            X_col = adata[:, gene].X
            
            # SAFE CHECK: Does it have a .toarray() method? (Works for all sparse types)
            if hasattr(X_col, "toarray"):
                data = X_col.toarray().flatten()
            else:
                # It's already a numpy array
                data = X_col.flatten()
            
            df[gene] = data

    # Save
    parquet_path = os.path.join(output_dir, "demo_data.parquet")
    df.to_parquet(parquet_path, compression='snappy')
    
    print(f"\nSUCCESS. Files saved to {output_dir}/")
    print(f"1. tissue.png")
    print(f"2. demo_data.parquet")

if __name__ == "__main__":
    preprocess_and_export()