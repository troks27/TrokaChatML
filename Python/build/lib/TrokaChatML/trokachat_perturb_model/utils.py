"""
Utilities Module

Utility functions for the cell communication perturbation package.
"""

import pandas as pd
import os
import numpy as np


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


def calculate_statistics(df, column):
    """
    Calculate basic statistics for a DataFrame column.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the data
    column : str
        Column name to calculate statistics for
        
    Returns:
    --------
    dict
        Dictionary containing statistics
    """
    stats = {
        'mean': df[column].mean(),
        'median': df[column].median(),
        'std': df[column].std(),
        'min': df[column].min(),
        'max': df[column].max(),
        'range': df[column].max() - df[column].min()
    }
    
    return stats


def compare_dataframes(df1, df2, key_columns, value_columns=None):
    """
    Compare two DataFrames and identify differences.
    
    Parameters:
    -----------
    df1 : pandas.DataFrame
        First DataFrame
    df2 : pandas.DataFrame
        Second DataFrame
    key_columns : list
        List of columns to use as keys for matching rows
    value_columns : list, optional
        List of columns to compare values, if None compares all columns
        
    Returns:
    --------
    tuple
        (in_df1_only, in_df2_only, different_values) - DataFrames with differences
    """
    # Set keys for comparison
    df1_key = df1.set_index(key_columns)
    df2_key = df2.set_index(key_columns)
    
    # Find keys that are only in df1 or only in df2
    df1_only_keys = df1_key.index.difference(df2_key.index)
    df2_only_keys = df2_key.index.difference(df1_key.index)
    
    # Extract rows that are unique to each DataFrame
    in_df1_only = df1.loc[df1.set_index(key_columns).index.isin(df1_only_keys)]
    in_df2_only = df2.loc[df2.set_index(key_columns).index.isin(df2_only_keys)]
    
    # If no value columns specified, compare all columns
    if value_columns is None:
        value_columns = [col for col in df1.columns if col not in key_columns]
    
    # Find common keys in both DataFrames
    common_keys = df1_key.index.intersection(df2_key.index)
    
    # For common keys, find rows where values differ
    different_values = []
    for key in common_keys:
        row1 = df1_key.loc[key][value_columns]
        row2 = df2_key.loc[key][value_columns]
        
        # Compare values and allow for floating point imprecision
        if not row1.equals(row2):
            if isinstance(row1, pd.Series) and isinstance(row2, pd.Series):
                non_matching_cols = [col for col in value_columns 
                                    if not np.isclose(row1[col], row2[col]) 
                                    if pd.api.types.is_numeric_dtype(row1[col]) 
                                    and pd.api.types.is_numeric_dtype(row2[col])
                                    else row1[col] != row2[col]]
                
                if non_matching_cols:
                    differences = {
                        'key': key,
                        'non_matching_columns': non_matching_cols,
                        'df1_values': row1[non_matching_cols].to_dict(),
                        'df2_values': row2[non_matching_cols].to_dict()
                    }
                    different_values.append(differences)
    
    return in_df1_only, in_df2_only, different_values


def batch_process(df, batch_size=1000, func=None, **kwargs):
    """
    Process a large DataFrame in batches to avoid memory issues.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to process
    batch_size : int, default=1000
        Size of each batch
    func : callable, default=None
        Function to apply to each batch
    **kwargs : dict
        Additional arguments to pass to func
        
    Returns:
    --------
    pandas.DataFrame
        Processed DataFrame
    """
    if func is None:
        return df
    
    result_dfs = []
    
    # Process DataFrame in batches
    for i in range(0, len(df), batch_size):
        batch = df.iloc[i:i+batch_size].copy()
        processed_batch = func(batch, **kwargs)
        result_dfs.append(processed_batch)
        
        print(f"Processed batch {i//batch_size + 1} of {(len(df) + batch_size - 1)//batch_size}")
    
    # Combine results
    return pd.concat(result_dfs, ignore_index=True)


def get_overlapping_genes(df1, df2, gene_column='Ligand'):
    """
    Get genes that appear in both DataFrames.
    
    Parameters:
    -----------
    df1 : pandas.DataFrame
        First DataFrame
    df2 : pandas.DataFrame
        Second DataFrame
    gene_column : str, default='Ligand'
        Column containing gene names
        
    Returns:
    --------
    set
        Set of genes present in both DataFrames
    """
    genes1 = set(df1[gene_column].unique())
    genes2 = set(df2[gene_column].unique())
    
    return genes1.intersection(genes2)
