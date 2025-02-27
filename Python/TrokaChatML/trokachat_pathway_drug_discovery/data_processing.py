"""
Data Processing Module

Functions for loading, processing, and saving data files.
"""

import os
import re
import pandas as pd


def extract_factor_number(filename):
    """
    Extract factor number from a filename.
    
    Parameters:
    -----------
    filename : str
        Filename that may contain 'Factor X'
        
    Returns:
    --------
    int or None
        The factor number if found, otherwise None
    """
    match = re.search(r'Factor (\d+)', filename)
    return int(match.group(1)) if match else None


def load_csv_file(file_path, sep=',', header='infer'):
    """
    Load a CSV file into a pandas DataFrame.
    
    Parameters:
    -----------
    file_path : str
        Path to the CSV file
    sep : str, default=','
        Delimiter to use
    header : int, 'infer' or None, default='infer'
        Row to use as header
        
    Returns:
    --------
    pandas.DataFrame
        Loaded DataFrame
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    if sep == 'auto':
        # Try to infer separator
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if '\t' in first_line:
                sep = '\t'
            elif ',' in first_line:
                sep = ','
            elif ';' in first_line:
                sep = ';'
    
    return pd.read_csv(file_path, sep=sep, header=header)


def save_to_excel(data_frames, output_path, sheet_names=None):
    """
    Save multiple DataFrames to an Excel file with multiple sheets.
    
    Parameters:
    -----------
    data_frames : list
        List of pandas DataFrames to save
    output_path : str
        Path where to save the Excel file
    sheet_names : list, optional
        List of sheet names, if None will use default names
        
    Returns:
    --------
    str
        Path to the saved Excel file
    """
    # Ensure directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create default sheet names if not provided
    if sheet_names is None:
        sheet_names = [f"Sheet{i+1}" for i in range(len(data_frames))]
    
    # Ensure we have the same number of sheet names as data_frames
    if len(sheet_names) != len(data_frames):
        raise ValueError("Number of sheet names must match number of DataFrames")
    
    # Write to Excel
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        for df, sheet_name in zip(data_frames, sheet_names):
            # Excel has a 31 character limit for sheet names
            sheet_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    return output_path


def find_files(directory, patterns=None, extensions=None, recursive=True):
    """
    Find files in a directory that match given patterns or extensions.
    
    Parameters:
    -----------
    directory : str
        Directory to search in
    patterns : list, optional
        List of regex patterns to match against filenames
    extensions : list, optional
        List of file extensions to match (e.g., ['.csv', '.tsv'])
    recursive : bool, default=True
        Whether to search recursively in subdirectories
        
    Returns:
    --------
    list
        List of found file paths
    """
    found_files = []
    
    if recursive:
        walker = os.walk(directory)
    else:
        walker = [(directory, [], os.listdir(directory))]
    
    for root, _, files in walker:
        for file in files:
            # Check extensions
            if extensions and not any(file.endswith(ext) for ext in extensions):
                continue
            
            # Check patterns
            if patterns:
                if not any(re.search(pattern, file) for pattern in patterns):
                    continue
            
            found_files.append(os.path.join(root, file))
    
    return found_files


def merge_data_files(file_paths, on=None, how='inner'):
    """
    Merge multiple data files into a single DataFrame.
    
    Parameters:
    -----------
    file_paths : list
        List of file paths to merge
    on : str or list, optional
        Column(s) to join on
    how : str, default='inner'
        Type of merge to perform
        
    Returns:
    --------
    pandas.DataFrame
        Merged DataFrame
    """
    dfs = []
    for file_path in file_paths:
        # Infer separator based on file extension
        sep = '\t' if file_path.endswith('.tsv') else ','
        df = load_csv_file(file_path, sep=sep)
        dfs.append(df)
    
    # Start with the first DataFrame
    result = dfs[0]
    
    # Merge with the rest
    for df in dfs[1:]:
        result = pd.merge(result, df, on=on, how=how)
    
    return result
