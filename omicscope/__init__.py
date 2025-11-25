from .core import create_app

def view(adata, port=8050, debug=False):
    app = create_app(data_input=adata, mode='cli', dataset_title="Local Analysis")
    print(f"Starting OmicScope on http://127.0.0.1:{port}")
    app.run_server(debug=debug, port=port)