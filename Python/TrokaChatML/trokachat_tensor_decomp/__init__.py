"""
Tensor Decomposition Subpackage

Tools for tensor decomposition, rank estimation, and analysis.
"""

from .data_loading import load_tensor_data, create_sparse_tensor, save_dimension_mappings
from .decomposition import perform_decomposition, save_decomposition_results
from .rank_estimation import calc_denoising_loss_parallel
from .visualization import plot_mse_boxplot

__all__ = [
    # Data loading
    'load_tensor_data',
    'create_sparse_tensor',
    'save_dimension_mappings',
    
    # Decomposition
    'perform_decomposition',
    'save_decomposition_results',
    
    # Rank estimation
    'calc_denoising_loss_parallel',
    
    # Visualization
    'plot_mse_boxplot'
]
