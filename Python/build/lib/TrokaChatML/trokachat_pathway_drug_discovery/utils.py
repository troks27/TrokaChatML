"""
Utilities Module

Utility functions for file system operations and other helper functions.
"""

import os
import re
import shutil
import pandas as pd


def ensure_directory_exists(directory):
    """
    Create directory if it doesn't exist.
    
    Parameters:
    -----------
    directory : str
        Directory path to create
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created directory: {directory}")


def print_directory_tree(startpath, max_level=None, include_files=True):
    """
    Print the directory structure as a tree.
    
    Parameters:
    -----------
    startpath : str
        Root directory to start from
    max_level : int, optional
        Maximum directory level to print
    include_files : bool, default=True
        Whether to include files in the output
        
    Returns:
    --------
    str
        String representation of the directory tree
    """
    output = []
    
    for root, dirs, files in os.walk(startpath):
        # Calculate the indentation level
        level = root.replace(startpath, '').count(os.sep)
        
        # Stop if we've reached max_level
        if max_level is not None and level > max_level:
            continue
        
        # Print the directory name
        indent = ' ' * 4 * level
        output.append(f"{indent}{os.path.basename(root)}/")
        
        # Print each file in the directory with additional indentation
        if include_files:
            subindent = ' ' * 4 * (level + 1)
            for f in sorted(files):
                output.append(f"{subindent}{f}")
    
    tree_str = '\n'.join(output)
    print(tree_str)
    return tree_str


def copy_file_to_destination(src, dst, create_dirs=True):
    """
    Copy a file to a destination, optionally creating directories.
    
    Parameters:
    -----------
    src : str
        Source file path
    dst : str
        Destination file path
    create_dirs : bool, default=True
        Whether to create destination directories if they don't exist
        
    Returns:
    --------
    str
        Path to the copied file
    """
    if not os.path.exists(src):
        raise FileNotFoundError(f"Source file not found: {src}")
    
    if create_dirs:
        ensure_directory_exists(os.path.dirname(dst))
    
    shutil.copy2(src, dst)
    return dst


def find_files_with_pattern(directory, pattern, recursive=True):
    """
    Find files matching a regex pattern.
    
    Parameters:
    -----------
    directory : str
        Directory to search in
    pattern : str
        Regex pattern to match against filenames
    recursive : bool, default=True
        Whether to search recursively in subdirectories
        
    Returns:
    --------
    list
        List of matching file paths
    """
    matching_files = []
    
    for root, _, files in os.walk(directory):
        for filename in files:
            if re.search(pattern, filename):
                matching_files.append(os.path.join(root, filename))
        
        if not recursive:
            break  # Stop after first level if not recursive
    
    return matching_files


def batch_process_files(files, process_function, **kwargs):
    """
    Apply a processing function to multiple files.
    
    Parameters:
    -----------
    files : list
        List of file paths to process
    process_function : callable
        Function to apply to each file, should accept file_path as first argument
    **kwargs : dict
        Additional keyword arguments to pass to process_function
        
    Returns:
    --------
    dict
        Dictionary with results for each file
    """
    results = {}
    
    for file_path in files:
        try:
            result = process_function(file_path, **kwargs)
            results[file_path] = result
        except Exception as e:
            results[file_path] = {'error': str(e)}
    
    return results


def filter_dataframe(df, filter_dict):
    """
    Filter a DataFrame based on column value conditions.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to filter
    filter_dict : dict
        Dictionary with column names as keys and filter conditions as values.
        For each column, the value can be:
        - A scalar value: exact match
        - A list/tuple: values to include (OR condition)
        - A dict with 'min'/'max' keys: range filter
        - A function: applied to column values, should return boolean Series
        
    Returns:
    --------
    pandas.DataFrame
        Filtered DataFrame
    """
    filtered_df = df.copy()
    
    for column, condition in filter_dict.items():
        if column not in filtered_df.columns:
            continue
        
        if callable(condition):
            # Function condition
            filtered_df = filtered_df[condition(filtered_df[column])]
        elif isinstance(condition, (list, tuple, set)):
            # List/tuple of values (OR condition)
            filtered_df = filtered_df[filtered_df[column].isin(condition)]
        elif isinstance(condition, dict):
            # Range filter with min/max
            if 'min' in condition:
                filtered_df = filtered_df[filtered_df[column] >= condition['min']]
            if 'max' in condition:
                filtered_df = filtered_df[filtered_df[column] <= condition['max']]
        else:
            # Scalar value (exact match)
            filtered_df = filtered_df[filtered_df[column] == condition]
    
    return filtered_df


def clean_column_names(df):
    """
    Clean column names in a DataFrame (strip whitespace, replace spaces with underscores).
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to clean column names
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with cleaned column names
    """
    df = df.copy()
    df.columns = [col.strip().replace(' ', '_').lower() for col in df.columns]
    return df
