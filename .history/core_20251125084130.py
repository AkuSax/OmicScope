import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import requests

# Data Normalization
def convert_anndata_to_df(adata):
    """
    Converts a live Scanpy AnnData object into the dataframe format 
    expected by the dashboard.
    """
    df = pd.DataFrame(index=adata.obs.index)
    
    # 1. Coordinates
    if 'spatial' in adata.obsm:
        df['spatial_x'] = adata.obsm['spatial'][:, 0]
        df['spatial_y'] = adata.obsm['spatial'][:, 1]
    else:
        df['spatial_x'] = adata.obsm['X_umap'][:, 0]
        df['spatial_y'] = adata.obsm['X_umap'][:, 1]
        
    df['umap_1'] = adata.obsm['X_umap'][:, 0]
    df['umap_2'] = adata.obsm['X_umap'][:, 1]
    
    # 2. Metadata
    if 'leiden' in adata.obs:
        df['cluster'] = adata.obs['leiden'].astype(str)
    else:
        df['cluster'] = "1"

    return df

# --- HELPER: GENE LOOKUP API ---
def fetch_gene_info(gene_symbol):
    """
    Queries MyGene.info for a gene description.
    """
    if not gene_symbol or gene_symbol == 'cluster':
        return "Select a gene to see details.", ""
    
    try:
        # Simple GET request to MyGene.info
        url = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=name,summary"
        response = requests.get(url, timeout=3) # 3s timeout to prevent hanging
        data = response.json()
        
        if 'hits' in data and len(data['hits']) > 0:
            hit = data['hits'][0]
            name = hit.get('name', 'Name unavailable')
            summary = hit.get('summary', 'No summary available for this gene.')
            return name, summary
        else:
            return "Unknown Gene", "Could not find information in database."
    except Exception as e:
        return "Error", "Could not connect to gene database."


# --- MAIN APP GENERATOR ---
def create_app(data_input, mode='web'):
    
    # 1. LOAD DATA
    if mode == 'web':
        df = pd.read_parquet(data_input)
    else:
        df = convert_anndata_to_df(data_input)

    # 2. PREPARE DROPDOWN OPTIONS
    exclude_cols = ['spatial_x', 'spatial_y', 'umap_1', 'umap_2', 'cluster']
    available_genes = [c for c in df.columns if c not in exclude_cols]
    available_genes.sort()

    # 3. INITIALIZE DASH APP
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])
    
    # 4. DEFINE LAYOUT
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
                value='cluster',
                clearable=False
            ),
            html.Br(),
            
            html.Label("View Coordinates:"),
            dbc.RadioItems(
                id='view-toggle',
                options=[
                    {'label': 'Spatial (Tissue)', 'value': 'spatial'},
                    {'label': 'UMAP (Embedding)', 'value': 'umap'}
                ],
                value='spatial',
                inline=True
            ),
            html.Br(),
            html.Small("Use your mouse to pan, zoom, and hover over cells.", className="text-muted"),
            
            html.Hr(),
            
            # --- NEW: GENE INFO CARD ---
            html.Div(id='gene-info-container', children=[
                html.H5(id='gene-name-display', className="text-primary"),
                html.P(id='gene-summary-display', className="small text-muted", 
                       style={'max-height': '200px', 'overflow-y': 'auto'})
            ])
        ],
        body=True,
        style={"height": "100vh", "background-color": "#f8f9fa"}
    )

    main_content = dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(dcc.Graph(id='main-scatter', style={"height": "90vh"}), width=12)
                ],
                align="center",
                className="h-100"
            )
        ],
        fluid=True
    )

    app.layout = dbc.Container(
        dbc.Row(
            [
                dbc.Col(sidebar, width=3, className="p-0"),
                dbc.Col(main_content, width=9)
            ]
        ),
        fluid=True
    )

    # 5. DEFINE CALLBACKS
    @app.callback(
        [Output('main-scatter', 'figure'),
         Output('gene-name-display', 'children'),
         Output('gene-summary-display', 'children')],
        [Input('color-dropdown', 'value'),
         Input('view-toggle', 'value')]
    )
    def update_graph_and_info(color_col, view_mode):
        # --- A. UPDATE INFO PANEL ---
        if color_col == 'cluster':
            gene_name = "Leiden Clustering"
            gene_summary = "Unsupervised clustering showing groups of cells with similar gene expression profiles."
        else:
            # Fetch real info
            gene_name, gene_summary = fetch_gene_info(color_col)

        # --- B. UPDATE GRAPH ---
        # Determine coordinates
        if view_mode == 'spatial':
            x_col, y_col = 'spatial_x', 'spatial_y'
            title = "Spatial Transcriptomics View"
        else:
            x_col, y_col = 'umap_1', 'umap_2'
            title = "UMAP Embedding"

        # Handle Categorical vs Continuous
        if color_col == 'cluster':
            fig = px.scatter(
                df, x=x_col, y=y_col, color=color_col,
                color_discrete_sequence=px.colors.qualitative.Plotly,
                hover_data=['cluster'],
                title=title
            )
        else:
            fig = px.scatter(
                df, x=x_col, y=y_col, color=color_col,
                color_continuous_scale=px.colors.sequential.Viridis,
                hover_data=['cluster', color_col],
                title=f"{title}: {color_col} Expression"
            )

        fig.update_layout(
            plot_bgcolor='white',
            dragmode='pan',
            legend_title_text='Group'
        )
        # Flip Y-axis for spatial (images typically have (0,0) at top-left)
        if view_mode == 'spatial':
            fig.update_yaxes(autorange="reversed") 
        
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        
        return fig, gene_name, gene_summary

    return app