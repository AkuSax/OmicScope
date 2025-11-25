import sys
import os

# Ensure the root directory is in the path so we can import 'omicscope'
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import create_app

# Path to the parquet file we generated
DATA_PATH = os.path.join(os.path.dirname(__file__), 'data', 'demo_data.parquet')

# Initialize the app in 'web' mode
app = create_app(data_input=DATA_PATH, mode='web')

# Expose server for Gunicorn
server = app.server

if __name__ == "__main__":
    app.run_server(debug=True, port=8050)