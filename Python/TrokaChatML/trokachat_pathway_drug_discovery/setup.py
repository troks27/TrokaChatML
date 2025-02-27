from setuptools import setup, find_packages

setup(
    name="trokachat_pathway_drug_discovery",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn",
        "kneed",
        "pybiomart",
        "tqdm",
        # Optional dependencies
        # "multixrank",
        # "chembl_downloader",
        # "networkx",
    ],
    extras_require={
        "full": [
            "multixrank",
            "chembl_downloader",
            "networkx",
        ],
    },
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for pathway and drug discovery analysis",
    keywords="bioinformatics, pathway, drug-discovery, multixrank",
    url="https://github.com/yourusername/trokachat_pathway_drug_discovery",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
