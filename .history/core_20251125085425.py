import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import requests
from PIL import Image
import os

# --- HELPER: DATA NORMALIZATION ---
def convert_anndata_to_df(adata):
    # (CLI Helper remains the same for now)
    df = pd.DataFrame(index=adata.obs.index)
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
    # (API Helper remains the same)
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
        df = pd.read_parquet(data_input)
        data_dir = os.path.dirname(data_input)
        img_path = os.path.join(data_dir, "tissue.png")
        try:
            tissue_img = Image.open(img_path)
            img_width, img_height = tissue_img.size
        except:
            tissue_img = None
            img_width, img_height = 1000, 1000
            
    else:
        df = convert_anndata_to_df(data_input)
        tissue_img = None
        img_width, img_height = 1000, 1000

    # 2. PREPARE DROPDOWN
    exclude_cols = ['spatial_x', 'spatial_y', 'umap_1', 'umap_2', 'cluster']
    available_genes = [c for c in df.columns if c not in exclude_cols]
    available_genes.sort()

    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])
    
    # 3. DEFINE LAYOUT
    sidebar = dbc.Card(
        [
            html.H2("OmicScope", className="display-6"),
            html.Hr(),
            html.P("Exploratory Analysis Tool", className="text-muted"),
            
            dbc.Label("Color By:", html_for="color-dropdown", className="fw-bold"),
            dcc.Dropdown(
                id='color-dropdown',
                options=[{'label': 'Cluster (Leiden)', 'value': 'cluster'}] + 
                        [{'label': g, 'value': g} for g in available_genes],
                value='cluster', clearable=False
            ),
            html.Br(),
            
            dbc.Label("Layout View:", className="fw-bold"),
            dbc.RadioItems(
                id='view-toggle',
                options=[
                    {'label': 'Spatial (Tissue)', 'value': 'spatial'},
                    {'label': 'UMAP (Embedding)', 'value': 'umap'}
                ],
                value='spatial', inline=True
            ),
            
            # --- NEW: Image Toggle Checkbox ---
            html.Div(id="image-toggle-container", children=[
                html.Br(),
                dbc.Checkbox(
                    id='show-image-check',
                    label="Show Histology Image",
                    value=True, # Default ON
                ),
            ], style={'display': 'block'}), # Initially visible

            html.Hr(),
            dbc.Label("Gene Info:", className="fw-bold"),
            html.Div(id='gene-info-container', children=[
                html.H5(id='gene-name-display', className="text-primary"),
                html.P(id='gene-summary-display', className="small text-muted", 
                       style={'max-height': '300px', 'overflow-y': 'auto'})
            ])
        ],
        body=True, style={"height": "100vh", "background-color": "#f8f9fa"}
    )

    main_content = dbc.Container(
        [dbc.Row([
            dbc.Col(
                # --- UPDATE: Enable Scroll Zoom in Config ---
                dcc.Graph(id='main-scatter', 
                          style={"height": "95vh"},
                          config={'scrollZoom': True, 'displayModeBar': True}
                         ), 
                width=12
            )
        ], align="center", className="h-100")],
        fluid=True, className="p-0"
    )

    app.layout = dbc.Container(
        dbc.Row([dbc.Col(sidebar, width=3, className="p-0"), dbc.Col(main_content, width=9, className="p-0")]),
        fluid=True, className="g-0"
    )

    # 4. DEFINE CALLBACKS
    # Callback to hide/show the image checkbox based on view mode
    @app.callback(
        Output('image-toggle-container', 'style'),
        Input('view-toggle', 'value')
    )
    def toggle_checkbox_visibility(view_mode):
        if view_mode == 'spatial':
            return {'display': 'block'}
        return {'display': 'none'}


    # Main Graph Callback
    @app.callback(
        [Output('main-scatter', 'figure'),
         Output('gene-name-display', 'children'),
         Output('gene-summary-display', 'children')],
        [Input('color-dropdown', 'value'),
         Input('view-toggle', 'value'),
         Input('show-image-check', 'value')] # NEW INPUT
    )
    def update_graph_and_info(color_col, view_mode, show_image):
        # A. Info Panel
        if color_col == 'cluster':
            g_name, g_sum = "Leiden Clustering", "Unsupervised clustering grouping cells by similar expression profiles."
        else:
            g_name, g_sum = fetch_gene_info(color_col)

        # B. Graph Logic setup
        if view_mode == 'spatial':
            x_col, y_col = 'spatial_x', 'spatial_y'
            title = "Spatial Transcriptomics View"
            # Smaller, more transparent spots for spatial view so image shows through
            marker_size = 4
            marker_opacity = 0.65 
        else:
            x_col, y_col = 'umap_1', 'umap_2'
            title = "UMAP Embedding"
            # Larger, more solid spots for UMAP
            marker_size = 5
            marker_opacity = 0.9

        # C. Create Base Scatter Plot
        if color_col == 'cluster':
            fig = px.scatter(df, x=x_col, y=y_col, color=color_col,
                             color_discrete_sequence=px.colors.qualitative.Plotly,
                             hover_data=['cluster'])
        else:
            fig = px.scatter(df, x=x_col, y=y_col, color=color_col,
                             color_continuous_scale=px.colors.sequential.Viridis,
                             hover_data=['cluster', color_col])
            
            # --- UPDATE: Add Colorbar Title ---
            fig.update_layout(coloraxis_colorbar=dict(
                title="Log-Normalized<br>Expression",
                tickfont=dict(size=10)
            ))

        # D. Refine Layout & Traces
        fig.update_layout(
            title=dict(text=title, x=0.05, y=0.98),
            plot_bgcolor='white',
            margin=dict(l=20, r=20, t=50, b=20),
            legend_title_text='Group'
        )
        fig.update_traces(marker=dict(size=marker_size, opacity=marker_opacity))
        fig.update_yaxes(scaleanchor="x", scaleratio=1)

        # E. Image Overlay Logic (Only if Spatial AND Checkbox is checked)
        if view_mode == 'spatial' and tissue_img and show_image:
            fig.add_layout_image(
                dict(
                    source=tissue_img,
                    xref="x", yref="y",
                    x=0, y=0,
                    sizex=img_width, sizey=img_height,
                    sizing="stretch",
                    # --- UPDATE: High opacity for background image ---
                    opacity=1.0, 
                    layer="below"
                )
            )
            fig.update_yaxes(autorange="reversed", visible=False) 
            fig.update_xaxes(visible=False)
            fig.update_layout(
                xaxis=dict(range=[0, img_width], constrain='domain'),
                yaxis=dict(range=[img_height, 0], constrain='domain')
            )
        elif view_mode == 'spatial':
             # Spatial coordinates but image is turned off
             fig.update_yaxes(autorange="reversed", visible=True)
             fig.update_xaxes(visible=True)
        else:
            # UMAP view
            fig.update_yaxes(autorange=True, visible=True)
            fig.update_xaxes(visible=True)

        return fig, g_name, g_sum

    return app