from setuptools import setup, find_packages

setup(
    name="trokachat_tensor_decomp",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "scikit-learn",
        "matplotlib",
        "seaborn",
        "lightgbm",
        "scikit-optimize",
        "scipy",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for cell-cell communication perturbation prediction",
    keywords="bioinformatics, cell-communication, machine-learning",
    url="https://github.com/yourusername/cell_comm_perturb",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
