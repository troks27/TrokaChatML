"""
Data Loading Module

Functions for loading tensor data from CSV files and converting to sparse tensors.
"""

import numpy as np
import pandas as pd
import sparse
from tensorly.contrib.sparse import tensor as tl_tensor
import os


def load_tensor_data(tensor_data_path, tensor_dimensions_path=None):
    """
    Load tensor data from CSV files.
    
    Parameters:
    -----------
    tensor_data_path : str
        Path to the CSV file containing tensor data
    tensor_dimensions_path : str, optional
        Path to the CSV file containing tensor dimensions
        
    Returns:
    --------
    pandas.DataFrame
        Loaded tensor data
    pandas.DataFrame or None
        Loaded tensor dimensions if path is provided, otherwise None
    """
    tensor_data = pd.read_csv(tensor_data_path)
    tensor_dimensions = None
    
    if tensor_dimensions_path:
        tensor_dimensions = pd.read_csv(tensor_dimensions_path)
    
    return tensor_data, tensor_dimensions


def create_sparse_tensor(tensor_data, dimension_columns=None):
    """
    Create a sparse tensor from tensor data.
    
    Parameters:
    -----------
    tensor_data : pandas.DataFrame
        DataFrame containing tensor data
    dimension_columns : list, optional
        List of column names corresponding to tensor dimensions.
        If None, defaults to ['Source', 'Target', 'Ligand_Receptor', 'Condition']
        
    Returns:
    --------
    tuple
        (sparse_tensor, tensor_shape, dim_mappings) - Sparse tensor, shape and dimension mappings
    """
    if dimension_columns is None:
        dimension_columns = ['Source', 'Target', 'Ligand_Receptor', 'Condition']
    
    # Get tensor shape based on unique values in each dimension
    tensor_shape = tuple(tensor_data[col].nunique() for col in dimension_columns)
    dim_mappings = []

    # Create mappings for each dimension based on the actual column names
    for col in dimension_columns:
        unique_values = tensor_data[col].unique()
        dim_mappings.append(dict(zip(unique_values, range(len(unique_values)))))

    # Initialize the tensor with zeros
    tensor_final = np.zeros(tensor_shape)

    # Populate the tensor with frequencies
    for _, row in tensor_data.iterrows():
        indices = tuple(dim_mappings[i][row[col]] for i, col in enumerate(dimension_columns))
        tensor_final[indices] = row['Freq']

    # Convert to sparse tensor
    sparse_coo = sparse.COO.from_numpy(tensor_final)
    sparse_tl_tensor = tl_tensor(sparse_coo)
    
    return sparse_tl_tensor, tensor_shape, dim_mappings


def save_dimension_mappings(dim_mappings, output_dir):
    """
    Save the mappings for each dimension to CSV files.
    
    Parameters:
    -----------
    dim_mappings : list
        List of dictionaries mapping dimension values to indices
    output_dir : str
        Directory to save the mapping files
        
    Returns:
    --------
    list
        Paths to the saved mapping files
    """
    os.makedirs(output_dir, exist_ok=True)
    mapping_files = []
    
    for i, mapping in enumerate(dim_mappings):
        mapping_df = pd.DataFrame(list(mapping.items()), columns=[f'Category_{i}', 'Index'])
        mapping_file_path = os.path.join(output_dir, f"mapping_{i}.csv")
        mapping_df.to_csv(mapping_file_path, index=False)
        mapping_files.append(mapping_file_path)
        
    return mapping_files
