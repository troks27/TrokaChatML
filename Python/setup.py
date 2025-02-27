from setuptools import setup, find_packages

setup(
    name="trokachatml",
    version="0.0.0.9000",
    packages=find_packages(),
    install_requires=[
        # Add your dependencies here, for example:
        # "numpy>=1.20.0",
        # "pandas>=1.3.0",
        # "scikit-learn>=1.0.0",
        # "tensorflow>=2.7.0",
    ],
    author="Michael Troka",
    author_email="troks27@upenn.edu",
    description="Analysis Pipeline for Single-Cell RNAseq to Drug Target Discovery",
    long_description="""TrokaChatML is an integrated computational framework designed for analyzing single-cell RNA sequencing (scRNAseq) data to uncover perturbed cell--cell communications and identify potential therapeutic targets. The pipeline combines differential expression analysis, machine learning-based noise filtration and perturbed communication prediction, tensor decomposition, and knowledge graph--driven drug target discovery, providing a robust and scalable approach from raw data processing to actionable biological insights.""",
    long_description_content_type="text/markdown",
    url="https://github.com/miketroka/TrokaChatML",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    python_requires=">=3.7",
)