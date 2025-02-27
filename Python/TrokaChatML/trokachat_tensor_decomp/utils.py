"""
Utilities Module

Utility functions for tensor decomposition operations.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tensorly.cp_tensor import cp_to_tensor


def ensure_directory_exists(directory):
    """
    Create directory if it does not exist.
    
    Parameters:
    -----------
    directory : str
        Directory path to create
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created directory: {directory}")


def batch_process_noise_levels(sparse_tl_tensor, noise_levels, rank_max, n_rep, output_dir):
    """
    Process multiple noise levels to determine tensor rank.
    
    Parameters:
    -----------
    sparse_tl_tensor : tensorly.contrib.sparse.tensor
        The sparse tensor to analyze
    noise_levels : list
        List of noise levels to test
    rank_max : int
        Maximum rank to test
    n_rep : int
        Number of repetitions for each noise level
    output_dir : str
        Directory to save the results
        
    Returns:
    --------
    dict
        Dictionary with results for each noise level
    """
    from .rank_estimation import calc_denoising_loss_parallel
    from .visualization import plot_mse_boxplot
    
    ensure_directory_exists(output_dir)
    results = {}
    
    for noise in noise_levels:
        print(f"\nProcessing noise level: {noise}")
        
        # Run the denoising loss calculation with parallelization
        avg_errors, err_summary_df = calc_denoising_loss_parallel(
            sparse_tl_tensor, 
            rank_max=rank_max, 
            noise=noise, 
            n_rep=n_rep
        )
        
        # Save the DataFrame to a CSV file
        csv_filename = os.path.join(output_dir, f"denoising_mask_noise={noise}_rank_estim.csv")
        err_summary_df.to_csv(csv_filename, index=False)
        
        # Create and save the plot
        output_filename = os.path.join(output_dir, f"denoising_mask_noise={noise}_rank_estim_plot.pdf")
        plot_mse_boxplot(err_summary_df=err_summary_df, output_filename=output_filename)
        
        results[noise] = {
            'avg_errors': avg_errors,
            'err_summary_df': err_summary_df,
            'csv_filename': csv_filename,
            'plot_filename': output_filename
        }
    
    return results


def tensor_summary_statistics(sparse_tensor):
    """
    Calculate summary statistics for a sparse tensor.
    
    Parameters:
    -----------
    sparse_tensor : tensorly.contrib.sparse.tensor
        The sparse tensor to analyze
        
    Returns:
    --------
    dict
        Dictionary containing summary statistics
    """
    # Get the tensor as a dense array (be careful with large tensors)
    dense_tensor = sparse_tensor.todense()
    
    # Calculate statistics
    stats = {
        'shape': sparse_tensor.shape,
        'dimensions': len(sparse_tensor.shape),
        'total_elements': np.prod(sparse_tensor.shape),
        'non_zero_elements': sparse_tensor.nnz,
        'sparsity': 1 - (sparse_tensor.nnz / np.prod(sparse_tensor.shape)),
        'min': float(np.min(dense_tensor)),
        'max': float(np.max(dense_tensor)),
        'mean': float(np.mean(dense_tensor)),
        'median': float(np.median(dense_tensor)),
        'std': float(np.std(dense_tensor))
    }
    
    return stats


def compare_original_reconstructed(original_tensor, factors, sample_size=1000):
    """
    Compare original tensor values with reconstructed values from CP decomposition.
    
    Parameters:
    -----------
    original_tensor : tensorly.contrib.sparse.tensor
        The original sparse tensor
    factors : tuple
        (weights, factors_list) - The weights and factors of the decomposition
    sample_size : int, default=1000
        Number of non-zero elements to sample for comparison
        
    Returns:
    --------
    tuple
        (rmse, r2, comparison_df) - RMSE, R² score, and DataFrame with comparison
    """
    from sklearn.metrics import mean_squared_error, r2_score
    
    # Reconstruct the tensor
    reconstructed = cp_to_tensor(factors)
    
    # Get coordinates of non-zero elements
    non_zero_coords = np.array(original_tensor.coords).T
    
    # Sample coordinates if there are more than sample_size
    if len(non_zero_coords) > sample_size:
        indices = np.random.choice(len(non_zero_coords), sample_size, replace=False)
        sampled_coords = non_zero_coords[indices]
    else:
        sampled_coords = non_zero_coords
    
    # Get original and reconstructed values
    original_values = [original_tensor[tuple(coord)] for coord in sampled_coords]
    reconstructed_values = [reconstructed[tuple(coord)] for coord in sampled_coords]
    
    # Calculate error metrics
    rmse = np.sqrt(mean_squared_error(original_values, reconstructed_values))
    r2 = r2_score(original_values, reconstructed_values)
    
    # Create comparison DataFrame
    comparison_df = pd.DataFrame({
        'Original': original_values,
        'Reconstructed': reconstructed_values,
        'Difference': np.array(original_values) - np.array(reconstructed_values)
    })
    
    return rmse, r2, comparison_df
