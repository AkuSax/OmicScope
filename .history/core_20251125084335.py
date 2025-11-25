import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import requests
from PIL import Image
import os

# --- HELPER: DATA NORMALIZATION ---
def convert_anndata_to_df(adata):
    """
    CLI Mode Helper: Converts AnnData to DataFrame.
    NOTE: For CLI mode with images, we would need to extract the image dynamically.
    For this demo, we assume the parquet/image pair exists.
    """
    df = pd.DataFrame(index=adata.obs.index)
    # Just basic extraction for now
    if 'spatial' in adata.obsm:
        df['spatial_x'] = adata.obsm['spatial'][:, 0]
        df['spatial_y'] = adata.obsm['spatial'][:, 1]
    else:
        df['spatial_x'] = adata.obsm['X_umap'][:, 0]
        df['spatial_y'] = adata.obsm['X_umap'][:, 1]
    
    df['umap_1'] = adata.obsm['X_umap'][:, 0]
    df['umap_2'] = adata.obsm['X_umap'][:, 1]
    
    if 'leiden' in adata.obs:
        df['cluster'] = adata.obs['leiden'].astype(str)
    else:
        df['cluster'] = "1"
    return df

# --- HELPER: GENE LOOKUP API ---
def fetch_gene_info(gene_symbol):
    if not gene_symbol or gene_symbol == 'cluster':
        return "Select a gene to see details.", ""
    try:
        url = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=name,summary"
        response = requests.get(url, timeout=3)
        data = response.json()
        if 'hits' in data and len(data['hits']) > 0:
            hit = data['hits'][0]
            return hit.get('name', 'Name unavailable'), hit.get('summary', 'No summary available.')
        else:
            return "Unknown Gene", "No info found."
    except:
        return "Error", "Connection failed."

# --- MAIN APP GENERATOR ---
def create_app(data_input, mode='web'):
    
    # 1. LOAD DATA & IMAGE
    if mode == 'web':
        # data_input is the path to the parquet file
        df = pd.read_parquet(data_input)
        
        # Look for the image in the same folder
        data_dir = os.path.dirname(data_input)
        img_path = os.path.join(data_dir, "tissue.png")
        
        try:
            # Load image and get dimensions
            tissue_img = Image.open(img_path)
            img_width, img_height = tissue_img.size
            print(f"Loaded tissue image: {img_width}x{img_height}")
        except Exception as e:
            print(f"Warning: Could not load image at {img_path}. {e}")
            tissue_img = None
            img_width, img_height = 1000, 1000 # Fallback
            
    else:
        df = convert_anndata_to_df(data_input)
        tissue_img = None # CLI image support skipped for MVP
        img_width, img_height = 1000, 1000

    # 2. PREPARE DROPDOWN
    exclude_cols = ['spatial_x', 'spatial_y', 'umap_1', 'umap_2', 'cluster']
    available_genes = [c for c in df.columns if c not in exclude_cols]
    available_genes.sort()

    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])
    
    sidebar = dbc.Card(
        [
            html.H2("OmicScope", className="display-6"),
            html.Hr(),
            html.P("Exploratory Analysis Tool", className="text-muted"),
            html.Label("Color By:"),
            dcc.Dropdown(
                id='color-dropdown',
                options=[{'label': 'Cluster (Leiden)', 'value': 'cluster'}] + 
                        [{'label': g, 'value': g} for g in available_genes],
                value='cluster', clearable=False
            ),
            html.Br(),
            html.Label("View Coordinates:"),
            dbc.RadioItems(
                id='view-toggle',
                options=[
                    {'label': 'Spatial (Tissue)', 'value': 'spatial'},
                    {'label': 'UMAP (Embedding)', 'value': 'umap'}
                ],
                value='spatial', inline=True
            ),
            html.Hr(),
            html.Div(id='gene-info-container', children=[
                html.H5(id='gene-name-display', className="text-primary"),
                html.P(id='gene-summary-display', className="small text-muted", 
                       style={'max-height': '200px', 'overflow-y': 'auto'})
            ])
        ],
        body=True, style={"height": "100vh", "background-color": "#f8f9fa"}
    )

    main_content = dbc.Container(
        [dbc.Row([dbc.Col(dcc.Graph(id='main-scatter', style={"height": "90vh"}), width=12)], 
                 align="center", className="h-100")],
        fluid=True
    )

    app.layout = dbc.Container(
        dbc.Row([dbc.Col(sidebar, width=3, className="p-0"), dbc.Col(main_content, width=9)]),
        fluid=True
    )

    @app.callback(
        [Output('main-scatter', 'figure'),
         Output('gene-name-display', 'children'),
         Output('gene-summary-display', 'children')],
        [Input('color-dropdown', 'value'),
         Input('view-toggle', 'value')]
    )
    def update_graph_and_info(color_col, view_mode):
        # A. Info Panel
        if color_col == 'cluster':
            g_name, g_sum = "Leiden Clustering", "Unsupervised clustering of gene expression."
        else:
            g_name, g_sum = fetch_gene_info(color_col)

        # B. Graph Logic
        if view_mode == 'spatial':
            x_col, y_col = 'spatial_x', 'spatial_y'
            title = "Spatial Transcriptomics View"
        else:
            x_col, y_col = 'umap_1', 'umap_2'
            title = "UMAP Embedding"

        # Create Scatter
        if color_col == 'cluster':
            fig = px.scatter(df, x=x_col, y=y_col, color=color_col,
                             color_discrete_sequence=px.colors.qualitative.Plotly,
                             hover_data=['cluster'], title=title)
        else:
            fig = px.scatter(df, x=x_col, y=y_col, color=color_col,
                             color_continuous_scale=px.colors.sequential.Viridis,
                             hover_data=['cluster', color_col], 
                             title=f"{title}: {color_col}")
            
        # Refine Layout
        fig.update_layout(plot_bgcolor='white', dragmode='pan', legend_title_text='Group')
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        
        # --- IMAGE OVERLAY LOGIC ---
        if view_mode == 'spatial' and tissue_img:
            # Add the image as a background layer
            fig.add_layout_image(
                dict(
                    source=tissue_img,
                    xref="x",
                    yref="y",
                    x=0,
                    y=0,
                    sizex=img_width,
                    sizey=img_height,
                    sizing="stretch",
                    opacity=0.5, # Semi-transparent so we can see spots clearly
                    layer="below"
                )
            )
            
            # Configure axes to match image dimensions
            # Standard image coords: (0,0) is top-left.
            # Plotly default: (0,0) is bottom-left.
            # We must reverse Y to match the image coordinate system.
            fig.update_yaxes(autorange="reversed")
            
            # Lock the view to the image size
            fig.update_xaxes(range=[0, img_width])
            fig.update_yaxes(range=[img_height, 0]) # Note the reversed order for range
            
            # Make spots smaller so we can see the tissue underneath
            fig.update_traces(marker=dict(size=5, opacity=0.8))
            
        else:
            # For UMAP, we don't want reversed axes or images
            fig.update_yaxes(autorange=True)
            fig.update_traces(marker=dict(size=5))

        return fig, g_name, g_sum

    return app