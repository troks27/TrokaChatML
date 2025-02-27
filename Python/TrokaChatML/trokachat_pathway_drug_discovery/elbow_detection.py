"""
Improved Elbow Detection Module

Functions for finding optimal cutoff points using the elbow method.
Adjusted to match the original implementation's approach to rolling averages.
"""

import os
import numpy as np
import pandas as pd
from kneed import KneeLocator
import matplotlib.pyplot as plt


def find_elbow_point(values, curve='convex', direction='decreasing', window_size=10):
    """
    Find the elbow point in a list of values using the KneeLocator algorithm,
    properly applying a rolling window to smooth the data.
    
    Parameters:
    -----------
    values : array-like
        Values to find elbow point in
    curve : str, default='convex'
        Type of curve ('convex' or 'concave')
    direction : str, default='decreasing'
        Direction of curve ('decreasing' or 'increasing')
    window_size : int, default=10
        Size of rolling window for smoothing
        
    Returns:
    --------
    tuple
        (elbow_index, elbow_value, values, rolling_avg)
    """
    # Convert to pandas Series to use rolling function
    values_series = pd.Series(values)
    
    # Create x-data
    x_data = np.arange(1, len(values_series) + 1)
    
    # Calculate rolling average, dropping NaN values 
    rolling_avg = values_series.rolling(window=window_size).mean().dropna()
    
    # Adjust x_data to match the size of rolling_avg (after dropping NaNs)
    # The indices of rolling_avg correspond to x_data[window_size-1:]
    adjusted_x_data = x_data[window_size-1:]
    
    # Find elbow using KneeLocator
    kneedle = KneeLocator(adjusted_x_data, rolling_avg, curve=curve, direction=direction)
    
    # Get elbow index
    elbow_idx = kneedle.elbow
    
    # If no elbow found, use default
    if elbow_idx is None:
        elbow_idx = 0
        elbow_value = None
    else:
        # Map the elbow back to the original data index
        # We need to ensure we're getting the correct value from the original data
        elbow_value = values_series.iloc[elbow_idx - 1]
    
    return elbow_idx, elbow_value, values_series.values, rolling_avg.values


def process_elbow_data(df, value_column, window_size=10, curve='convex', direction='decreasing'):
    """
    Process DataFrame by sorting it by a value column and finding the elbow point.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to process
    value_column : str
        Column name to sort by and find elbow point
    window_size : int, default=10
        Size of rolling window for smoothing
    curve : str, default='convex'
        Type of curve for KneeLocator
    direction : str, default='decreasing'
        Direction of curve for KneeLocator
        
    Returns:
    --------
    tuple
        (elbow_idx, elbow_value, sorted_df, filtered_df)
    """
    # Check if value_column exists
    if value_column not in df.columns:
        raise ValueError(f"Column '{value_column}' not found in DataFrame")
    
    # Sort DataFrame by value_column
    sorted_df = df.sort_values(by=value_column, ascending=(direction != 'decreasing')).reset_index(drop=True)
    
    # Find elbow point
    elbow_idx, elbow_value, _, _ = find_elbow_point(
        sorted_df[value_column], 
        curve=curve, 
        direction=direction, 
        window_size=window_size
    )
    
    # Filter DataFrame based on elbow point
    if elbow_value is not None:
        if direction == 'decreasing':
            filtered_df = sorted_df[sorted_df[value_column] >= elbow_value]
        else:
            filtered_df = sorted_df[sorted_df[value_column] <= elbow_value]
    else:
        filtered_df = sorted_df
    
    return elbow_idx, elbow_value, sorted_df, filtered_df


def plot_elbow_detection(values, elbow_idx, window_size=10, title='Elbow Detection',
                        save_path=None, figsize=(10, 5), curve='convex', direction='decreasing'):
    """
    Plot the values with the detected elbow point and rolling average.
    
    Parameters:
    -----------
    values : array-like
        Values to plot
    elbow_idx : int
        Index of the detected elbow point
    window_size : int, default=10
        Size of rolling window used for smoothing
    title : str, default='Elbow Detection'
        Title for the plot
    save_path : str, optional
        Path to save the plot, if None plot is not saved
    figsize : tuple, default=(10, 5)
        Figure size as (width, height) in inches
    curve : str, default='convex'
        Type of curve ('convex' or 'concave')
    direction : str, default='decreasing'
        Direction of curve ('decreasing' or 'increasing')
        
    Returns:
    --------
    matplotlib.figure.Figure
        The figure object
    """
    # Convert to pandas Series
    values_series = pd.Series(values)
    
    # Create x-data
    x_data = np.arange(1, len(values_series) + 1)
    
    # Calculate rolling average, dropping NaN values
    rolling_avg = values_series.rolling(window=window_size).mean().dropna()
    adjusted_x_data = x_data[window_size-1:]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot original data
    ax.plot(x_data, values_series, label='Original Data', color='blue')
    
    # Plot rolling average
    ax.plot(adjusted_x_data, rolling_avg, label='Rolling Average', color='green')
    
    # Mark elbow point
    if elbow_idx > 0:
        ax.axvline(x=elbow_idx, color='red', linestyle='--', label='Elbow Point')
        # Adding a dot at the elbow point on the rolling average line
        y_value = rolling_avg.iloc[adjusted_x_data.tolist().index(elbow_idx)] if elbow_idx in adjusted_x_data else None
        if y_value is not None:
            ax.scatter([elbow_idx], [y_value], color='red', s=50, zorder=5)
    
    # Set labels and title
    ax.set_xlabel("Index (Sorted Order)")
    ax.set_ylabel("Value")
    ax.set_title(title)
    ax.legend()
    
    # Save figure if path provided
    if save_path is not None:
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def batch_process_files(file_paths, value_column, output_dir=None, window_size=10, 
                       curve='convex', direction='decreasing', plot=True):
    """
    Process multiple files to find elbow points and filter data.
    
    Parameters:
    -----------
    file_paths : list
        List of file paths to process
    value_column : str
        Column name to sort by and find elbow point
    output_dir : str, optional
        Directory to save output files, if None will use input file directory
    window_size : int, default=10
        Size of rolling window for smoothing
    curve : str, default='convex'
        Type of curve for KneeLocator
    direction : str, default='decreasing'
        Direction of curve for KneeLocator
    plot : bool, default=True
        Whether to generate and save plots
        
    Returns:
    --------
    dict
        Dictionary with results for each file
    """
    results = {}
    
    for file_path in file_paths:
        # Determine file type and load data
        sep = '\t' if file_path.endswith('.tsv') else ','
        df = pd.read_csv(file_path, sep=sep)
        
        try:
            # Process data
            elbow_idx, elbow_value, sorted_df, filtered_df = process_elbow_data(
                df, value_column, window_size, curve, direction
            )
            
            # Determine output directory
            if output_dir is None:
                output_dir = os.path.dirname(file_path)
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate output file paths
            basename = os.path.splitext(os.path.basename(file_path))[0]
            filtered_path = os.path.join(output_dir, f"{basename}_filtered.tsv")
            plot_path = os.path.join(output_dir, f"{basename}_elbow_plot.pdf")
            
            # Save filtered data
            filtered_df.to_csv(filtered_path, sep='\t', index=False)
            
            # Generate plot if requested
            if plot:
                fig = plot_elbow_detection(
                    sorted_df[value_column], elbow_idx, window_size,
                    title=f"Elbow Detection for {basename}",
                    save_path=plot_path,
                    curve=curve,
                    direction=direction
                )
            
            # Store results
            results[file_path] = {
                'elbow_idx': elbow_idx,
                'elbow_value': elbow_value,
                'filtered_rows': len(filtered_df),
                'total_rows': len(df),
                'filtered_path': filtered_path,
                'plot_path': plot_path if plot else None
            }
            
        except Exception as e:
            results[file_path] = {'error': str(e)}
    
    return results