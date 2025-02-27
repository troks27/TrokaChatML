"""
Decomposition Module

Functions for tensor decomposition and saving decomposition results.
"""

import os
import numpy as np
import pandas as pd
from tensorly.contrib.sparse.decomposition import non_negative_parafac
from tensorly.cp_tensor import cp_to_tensor


def perform_decomposition(sparse_tensor, rank, init='random', n_iter_max=1000, verbose=False):
    """
    Perform non-negative CP decomposition on a sparse tensor.
    
    Parameters:
    -----------
    sparse_tensor : tensorly.contrib.sparse.tensor
        The sparse tensor to decompose
    rank : int
        The rank of the decomposition
    init : str, default='random'
        Initialization method for the decomposition
    n_iter_max : int, default=1000
        Maximum number of iterations
    verbose : bool, default=False
        Whether to print progress information
        
    Returns:
    --------
    tuple
        (weights, factors) - The weights and factors of the decomposition
    """
    # Perform non-negative CP decomposition
    factors = non_negative_parafac(
        sparse_tensor, 
        rank=rank, 
        init=init, 
        verbose=verbose, 
        n_iter_max=n_iter_max
    )
    
    return factors


def save_decomposition_results(factors, output_dir):
    """
    Save tensor decomposition results (weights and factors) to CSV files.
    
    Parameters:
    -----------
    factors : tuple
        (weights, factors_list) - The weights and factors of the decomposition
    output_dir : str
        Directory to save the decomposition results
        
    Returns:
    --------
    dict
        Dictionary with paths to the saved files
    """
    weights, factors_list = factors
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    saved_files = {}
    
    # Convert weights to a dense array
    weights_array = np.array(weights.todense()).flatten()
    
    # Save the weights if they are non-zero
    if np.count_nonzero(weights_array) > 0:
        weights_file_path = os.path.join(output_dir, "weights.csv")
        pd.DataFrame(weights_array, columns=['Weight']).to_csv(weights_file_path, index=False)
        saved_files['weights'] = weights_file_path
    
    # Save the factor matrices
    factor_files = []
    for i, factor in enumerate(factors_list):
        try:
            # Convert factor to a dense array
            dense_factor = factor.todense()  # Convert COO to dense
            # Convert each dense factor to DataFrame and save
            df_factor = pd.DataFrame(dense_factor)
            factor_file_path = os.path.join(output_dir, f"factor_{i}.csv")
            df_factor.to_csv(factor_file_path, index=False)
            factor_files.append(factor_file_path)
        except Exception as e:
            print(f"Error saving factor {i}: {e}")
    
    saved_files['factors'] = factor_files
    
    return saved_files


def reconstruct_tensor(factors):
    """
    Reconstruct a tensor from its CP decomposition factors.
    
    Parameters:
    -----------
    factors : tuple
        (weights, factors_list) - The weights and factors of the decomposition
        
    Returns:
    --------
    numpy.ndarray
        The reconstructed tensor
    """
    weights, factors_list = factors
    
    # Use tensorly's cp_to_tensor function to reconstruct
    reconstructed = cp_to_tensor((weights, factors_list))
    
    return reconstructed
