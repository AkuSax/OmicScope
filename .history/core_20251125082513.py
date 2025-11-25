import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import numpy as np

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
        # Fallback if no spatial info
        df['spatial_x'] = adata.obsm['X_umap'][:, 0]
        df['spatial_y'] = adata.obsm['X_umap'][:, 1]
        
    df['umap_1'] = adata.obsm['X_umap'][:, 0]
    df['umap_2'] = adata.obsm['X_umap'][:, 1]
    
    # Metadata
    if 'leiden' in adata.obs:
        df['cluster'] = adata.obs['leiden'].astype(str)
    else:
        df['cluster'] = "1"

    return df

def create_app(data_input, mode='web'):
    """
    Args:
        data_input: Either a filepath (str) for Web mode, or AnnData for CLI mode.
        mode: 'web' or 'cli'
    """
    
    # 1. LOAD DATA
    if mode == 'web':
        # Load the parquet file we just generated
        df = pd.read_parquet(data_input)
    else:
        # CLI Mode: data_input is an AnnData object
        # Note: For the MVP CLI, we are doing a basic conversion.
        # Ideally, we would pass the whole adata and query it dynamically.
        df = convert_anndata_to_df(data_input)

    # 2. PREPARE DROPDOWN OPTIONS
    # Columns that are NOT coordinates or clusters are treated as Genes
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
                value='cluster', # Default value
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
            html.Small("Use your mouse to pan, zoom, and hover over cells.", className="text-muted")
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
        Output('main-scatter', 'figure'),
        [Input('color-dropdown', 'value'),
         Input('view-toggle', 'value')]
    )
    def update_graph(color_col, view_mode):
        # Determine coordinates
        if view_mode == 'spatial':
            x_col, y_col = 'spatial_x', 'spatial_y'
            title = "Spatial Transcriptomics View"
        else:
            x_col, y_col = 'umap_1', 'umap_2'
            title = "UMAP Embedding"

        # Handle Categorical (Cluster) vs Continuous (Gene)
        if color_col == 'cluster':
            fig = px.scatter(
                df, x=x_col, y=y_col, color=color_col,
                color_discrete_sequence=px.colors.qualitative.Plotly,
                hover_data=['cluster'],
                title=title
            )
        else:
            # Gene Expression (Continuous)
            fig = px.scatter(
                df, x=x_col, y=y_col, color=color_col,
                color_continuous_scale=px.colors.sequential.Viridis,
                hover_data=['cluster', color_col],
                title=f"{title}: {color_col} Expression"
            )

        # Clean up the graph aesthetics
        fig.update_layout(
            plot_bgcolor='white',
            dragmode='pan',
            legend_title_text='Group'
        )
        fig.update_yaxes(scaleanchor="x", scaleratio=1) # Keep aspect ratio for spatial
        
        return fig

    return app