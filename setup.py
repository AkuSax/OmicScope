from setuptools import setup, find_packages

setup(
    name="omicscope",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "dash",
        "dash-bootstrap-components",
        "pandas",
        "scanpy",
        "requests",
        "Pillow"
    ],
    packages=find_packages(exclude=["web_demo", "web_demo.*"]),
    author="Akul Saxena",
    description="Interactive EDA tool for Spatial Transcriptomics",
)