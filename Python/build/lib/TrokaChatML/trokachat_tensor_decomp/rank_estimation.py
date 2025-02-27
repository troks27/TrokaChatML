"""
Rank Estimation Module

Functions for estimating the rank of a tensor through various methods,
including denoising loss calculation.
"""

import numpy as np
import random
import math
import sparse
import multiprocessing
import time
import pandas as pd
import os


def compute_error_for_rank(rank, masked_tensor, original_tensor, mask, sampled_locs, method="MSE", eps=1E-10):
    """
    Compute the error for a specific rank decomposition.
    
    Parameters:
    -----------
    rank : int
        The rank to test
    masked_tensor : sparse.COO
        The tensor with masked elements
    original_tensor : sparse.COO
        The original tensor
    mask : sparse.COO
        The mask indicating which elements were masked
    sampled_locs : list
        The locations of masked elements
    method : str, default="MSE"
        The error metric to use
    eps : float, default=1E-10
        Small value to avoid division by zero
    
    Returns:
    --------
    tuple
        (rank, error) - The rank and corresponding error
    """
    from tensorly.contrib.sparse.decomposition import non_negative_parafac
    
    try:
        # Perform non-negative PARAFAC decomposition on the masked tensor
        weights, factors = non_negative_parafac(masked_tensor, rank=rank, init='random', 
                                                verbose=False, n_iter_max=100, tol=1e-6)
        
        # Reconstruct the tensor
        from tensorly.cp_tensor import cp_to_tensor
        reconstructed = cp_to_tensor((weights, factors))
        
        # Get the values at the masked locations from both tensors
        original_values = [original_tensor[tuple(loc)] for loc in sampled_locs]
        reconstructed_values = [reconstructed[tuple(loc)] for loc in sampled_locs]
        
        # Calculate the error based on the selected method
        if method.upper() == "MSE":
            error = np.mean([(ov - rv)**2 for ov, rv in zip(original_values, reconstructed_values)])
        elif method.upper() == "RMSE":
            error = np.sqrt(np.mean([(ov - rv)**2 for ov, rv in zip(original_values, reconstructed_values)]))
        elif method.upper() == "MAE":
            error = np.mean([abs(ov - rv) for ov, rv in zip(original_values, reconstructed_values)])
        elif method.upper() == "MAPE":
            error = np.mean([abs((ov - rv) / (ov + eps)) for ov, rv in zip(original_values, reconstructed_values)]) * 100
        else:
            error = np.mean([(ov - rv)**2 for ov, rv in zip(original_values, reconstructed_values)])  # Default to MSE
            
        return rank, error
    except Exception as e:
        print(f"Error computing for rank {rank}: {e}")
        return rank, None


def calc_denoising_loss_parallel(sparse_tensor, rank_max=10, noise=0.3, n_rep=2, method="MSE", eps=1E-10, output_dir=None):
    """
    Calculate denoising loss to determine optimal tensor rank.
    
    Parameters:
    -----------
    sparse_tensor : tensorly.contrib.sparse.tensor
        The sparse tensor to analyze
    rank_max : int, default=10
        Maximum rank to test
    noise : float, default=0.3
        Proportion of non-zero elements to mask
    n_rep : int, default=2
        Number of repetitions for more robust estimation
    method : str, default="MSE"
        Error metric to use
    eps : float, default=1E-10
        Small value to avoid division by zero
    output_dir : str, optional
        Directory to save results
        
    Returns:
    --------
    tuple
        (avg_errors, err_summary_df) - Average errors by rank and detailed DataFrame
    """
    err_summary = []
    
    # Ensure we're working with a sparse tensor
    original_tensor = sparse_tensor.copy()
    non_zero_indices = np.array(sparse_tensor.coords).T
    print(f"Total non-zero elements: {len(non_zero_indices)}")
    
    for n in range(n_rep):
        print(f"\nStarting iteration {n+1}/{n_rep}")
        
        # Randomly select indices to mask
        sample_indices = random.sample(range(len(non_zero_indices)), math.floor(len(non_zero_indices) * noise))
        sampled_locs = non_zero_indices[sample_indices]
        print(f"Number of elements to mask: {len(sample_indices)}")
        
        # Create mask with the same sparsity structure
        mask_coords = []
        mask_data = []

        masked_data = sparse_tensor.data.copy()
        
        for idx in sample_indices:
            loc = tuple(non_zero_indices[idx])
            mask_coords.append(loc)
            mask_data.append(0)
            masked_data[idx] = 0  # Mask the tensor by setting selected elements to 0

        # Create new masked tensor and mask in sparse format
        masked_tensor = sparse.COO(coords=sparse_tensor.coords, data=masked_data, shape=sparse_tensor.shape)
        mask = sparse.COO(coords=np.array(mask_coords).T, data=np.array(mask_data), shape=sparse_tensor.shape)

        print("Masked tensor created")
        print(f"Masked tensor non-zero elements: {masked_tensor.nnz}")
        print(f"Mask non-zero elements: {mask.nnz}")
        
        # Use multiprocessing to compute error for each rank
        ranks = range(1, rank_max + 1)
        PROCESSES = max(1, multiprocessing.cpu_count() - 1)  # Ensure at least one process is used

        with multiprocessing.Pool(PROCESSES) as pool:
            results = [pool.apply_async(compute_error_for_rank, 
                                       (rank, masked_tensor, original_tensor, mask, sampled_locs, method, eps)) 
                      for rank in ranks]

            errors = []
            for r in results:
                rank, error = r.get()
                if error is not None:
                    errors.append((rank, error))
                    print(f"Rank {rank} completed. Error: {error}")
                else:
                    print(f"Error for rank {rank}")

        # Sort errors by rank to ensure order
        errors.sort(key=lambda x: x[0])
        ordered_errors = [error for _, error in errors]

        err_summary.append(ordered_errors)
        print(f"Completed iteration {n+1}/{n_rep}")

    # Create DataFrame with results
    err_summary_df = pd.DataFrame(err_summary, columns=[f'Rank_{i+1}' for i in range(rank_max)])
    
    # Save results if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        csv_filename = os.path.join(output_dir, f"denoising_mask_noise={noise}_rank_estim.csv")
        err_summary_df.to_csv(csv_filename, index=False)
        print(f"Results saved to {csv_filename}")
    
    print("Denoising loss calculation completed")
    return err_summary_df.mean(axis=0).tolist(), err_summary_df


def find_optimal_rank(err_summary_df, method='min_mean'):
    """
    Determine the optimal rank based on error metrics.
    
    Parameters:
    -----------
    err_summary_df : pandas.DataFrame
        DataFrame containing error metrics for different ranks
    method : str, default='min_mean'
        Method to determine optimal rank:
        'min_mean' - rank with minimum mean error
        'elbow' - rank at the elbow point of the error curve
        
    Returns:
    --------
    int
        Optimal rank
    """
    means = err_summary_df.mean()
    
    if method == 'min_mean':
        # Return rank with minimum mean error
        return np.argmin(means) + 1
    
    elif method == 'elbow':
        # Implement elbow method using second derivative
        means_array = means.values
        # Calculate first differences
        first_diff = np.diff(means_array)
        # Calculate second differences
        second_diff = np.diff(first_diff)
        # Find the index of the maximum second difference
        elbow_idx = np.argmax(np.abs(second_diff)) + 2  # +2 because of two diff operations
        
        return elbow_idx
    
    else:
        # Default to minimum mean
        return np.argmin(means) + 1
