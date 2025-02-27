"""
Cell Communication Perturbation Prediction Package

This package provides tools for analyzing and predicting perturbations in
cell-cell communication pathways, particularly in the context of comparing
control and experimental conditions.
"""

# Import main functionality to make it easily accessible
from .data_processing import (
    process_pct_filter,
    preprocess_datasets,
    train_set_create,
    test_set_create,
    stratified_split_handle_rare_classes
)

from .visualization import (
    plot_venn_diagram,
    plot_feature_importance,
    plot_predicted_vs_actual,
    plot_residuals_distribution,
    plot_predicted_vs_actual_values,
    plot_qq,
    plot_standardized_residuals_vs_predicted
)

from .modeling import (
    create_lgb_dataset,
    perform_kfold_cv,
    save_model,
    load_model,
    evaluate_model,
    evaluate_model_with_plots,
    train_final_model,
    preprocess_test_dataset
)

# Define package info
__version__ = "0.1.0"
__author__ = "Your Name"
