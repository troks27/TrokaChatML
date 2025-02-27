"""
Visualization Module

Functions for visualizing tensor decomposition results and rank estimation.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


def plot_rank_errors(avg_errors, rank_max, title="Rank Determination", log_scale=True, save_path=None):
    """
    Plot average errors for different ranks.
    
    Parameters:
    -----------
    avg_errors : list
        List of average errors for each rank
    rank_max : int
        Maximum rank tested
    title : str, default="Rank Determination"
        Plot title
    log_scale : bool, default=True
        Whether to use logarithmic scale for y-axis
    save_path : str, optional
        Path to save the plot
        
    Returns:
    --------
    tuple
        (fig, ax) - Figure and axis objects
    """
    rank_list = list(range(1, rank_max + 1))
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.plot(rank_list, avg_errors, marker='o')
    
    if log_scale:
        ax.set_yscale('log')
    
    ax.set_xlabel("Rank")
    ax.set_ylabel("Loss")
    ax.set_title(title)
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig, ax


def plot_mse_boxplot(csv_filename=None, err_summary_df=None, output_filename=None):
    """
    Create a boxplot of MSE values for different ranks.
    
    Parameters:
    -----------
    csv_filename : str, optional
        Path to CSV file with rank error data
    err_summary_df : pandas.DataFrame, optional
        DataFrame with rank error data (alternative to csv_filename)
    output_filename : str, optional
        Path to save the output plot
        
    Returns:
    --------
    tuple
        (fig, ax) - Figure and axis objects
    """
    # Set the style
    plt.style.use("seaborn-v0_8-whitegrid")
    
    # Define color palette
    palette = sns.color_palette()
    box_color = 'black'
    mean_line_color = palette[0]  # Blue
    scatter_color = palette[7]    # Gray
    min_mse_line_color = palette[3]  # Red

    # Load the data
    if err_summary_df is None:
        if csv_filename is None:
            raise ValueError("Either csv_filename or err_summary_df must be provided")
        df = pd.read_csv(csv_filename)
    else:
        df = err_summary_df
    
    ranks = df.columns

    # Calculate the mean and standard error for each rank
    means = df.mean()
    std_errors = df.sem()

    # Identify the rank with the minimum mean MSE
    min_rank = np.argmin(means) + 1

    # Create a figure and axis for plotting
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    # Plot the boxplot for each rank
    bp = ax.boxplot([df[rank] for rank in ranks], positions=np.arange(1, len(ranks) + 1),
                    showfliers=False, patch_artist=True)

    # Customize boxplot colors
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=box_color)

    for patch in bp['boxes']:
        patch.set(facecolor='none')  # No fill for boxes

    # Add jittered scatter points for each rank
    jitter_strength = 0.1
    for i, rank in enumerate(ranks):
        jittered_x = np.random.normal(loc=i+1, scale=jitter_strength, size=len(df[rank]))
        ax.scatter(jittered_x, df[rank], alpha=0.4, color=scatter_color, s=15, label='_nolegend_')

    # Plot the line connecting the mean values at each rank
    ax.plot(np.arange(1, len(ranks) + 1), means, color=mean_line_color, label='Mean MSE', lw=2)

    # Add a vertical line at the rank with the minimum MSE
    ax.axvline(x=min_rank, color=min_mse_line_color, linestyle='--', linewidth=2,
               label=f'Minimum MSE at Rank {min_rank}')

    # Customize the plot
    ax.set_xlabel('Rank', fontsize=14, fontweight='bold')
    ax.set_ylabel('Mean Squared Error (MSE)', fontsize=14, fontweight='bold')
    ax.set_title('MSE vs Rank', fontsize=16, fontweight='bold')
    ax.set_xticks(np.arange(1, len(ranks) + 1))
    ax.set_xticklabels(np.arange(1, len(ranks) + 1))
    ax.tick_params(axis='both', which='major', labelsize=10)

    # Set y-axis limits and ticks based on data
    y_max = means.max() * 1.2
    y_ticks = np.linspace(0, y_max, 10)
    ax.set_ylim(0, y_max)
    ax.set_yticks(y_ticks)

    # Enhance the legend
    ax.legend(fontsize=10, loc='upper left')

    # Adjust layout and save the plot if output filename is provided
    plt.tight_layout()
    if output_filename:
        plt.savefig(output_filename, format='pdf', dpi=300, bbox_inches='tight')
    
    return fig, ax


def plot_factor_heatmap(factor_matrix, dimension_name=None, save_path=None):
    """
    Create a heatmap visualization of a factor matrix.
    
    Parameters:
    -----------
    factor_matrix : numpy.ndarray
        Factor matrix to visualize
    dimension_name : str, optional
        Name of the dimension for the plot title
    save_path : str, optional
        Path to save the heatmap
        
    Returns:
    --------
    tuple
        (fig, ax) - Figure and axis objects
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create heatmap
    heatmap = sns.heatmap(factor_matrix, cmap="YlGnBu", ax=ax)
    
    # Set title and labels
    title = f"Factor Matrix Heatmap" if dimension_name is None else f"Factor Matrix for {dimension_name}"
    ax.set_title(title)
    ax.set_xlabel("Component")
    ax.set_ylabel("Index")
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig, ax


def plot_component_distribution(factor_matrix, component_idx, dimension_name=None, save_path=None):
    """
    Plot the distribution of values in a specific component of a factor matrix.
    
    Parameters:
    -----------
    factor_matrix : numpy.ndarray
        Factor matrix
    component_idx : int
        Index of the component to plot
    dimension_name : str, optional
        Name of the dimension for the plot title
    save_path : str, optional
        Path to save the plot
        
    Returns:
    --------
    tuple
        (fig, ax) - Figure and axis objects
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Get the component values
    component_values = factor_matrix[:, component_idx]
    
    # Create histogram and density plot
    sns.histplot(component_values, kde=True, ax=ax)
    
    # Set title and labels
    title = f"Component {component_idx+1} Distribution"
    if dimension_name:
        title += f" ({dimension_name})"
    
    ax.set_title(title)
    ax.set_xlabel(f"Component {component_idx+1} Value")
    ax.set_ylabel("Frequency")
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig, ax
