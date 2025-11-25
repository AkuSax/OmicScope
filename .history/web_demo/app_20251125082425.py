import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core import create_app

DATA_PATH = os.path.join(os.path.dirname(__file__), 'data', 'demo_data.parquet')
app = create_app(data_input=DATA_PATH, mode='web')
server = app.server

if __name__ == "__main__":
    app.run(debug=True, port=8050)