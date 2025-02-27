"""
Pathway Drug Discovery Package

A package for analyzing pathway and drug data, determining significant factors,
mapping between species, and enriching drug information.
"""

from .data_processing import load_csv_file, save_to_excel, find_files
from .elbow_detection import find_elbow_point, process_elbow_data
from .gene_mapping import get_mouse_mapping, collect_gene_names, transform_gene_file
from .multixrank_utils import create_multixrank_config, run_multixrank_analysis
from .drug_enrichment import query_drug_info, enrich_drug_data, update_drug_info
from .visualization import plot_elbow_detection, plot_factor_heatmap
from .utils import print_directory_tree, ensure_directory_exists

__version__ = "0.1.0"
__author__ = "Your Name"