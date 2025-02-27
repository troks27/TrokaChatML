"""
Visualization Module

Functions for creating visualizations and plots.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_elbow_detection(values, elbow_idx, window_size=10, title=None, 
                        save_path=None, figsize=(10, 6), min_index=3):
    """
    Plot values with elbow point detection.
    
    Parameters:
    -----------
    values : array-like
        Values to plot
    elbow_idx : int
        Index of elbow point
    window_size : int, default=10
        Size of rolling window for smoothing
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save plot
    figsize : tuple, default=(10, 6)
        Figure size
    min_index : int, default=3
        Minimum index to consider for elbow (for visualization)
        
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Convert to numpy array if not already
    values = np.array(values)
    
    # Create x values
    x_data = np.arange(1, len(values) + 1)
    
    # Calculate rolling average
    rolling_avg = pd.Series(values).rolling(window=window_size).mean().fillna(method='bfill').values
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot original data
    ax.plot(x_data, values, label='Original Data', color='blue')
    
    # Plot rolling average
    ax.plot(x_data, rolling_avg, label='Rolling Average', color='green')
    
    # Mark elbow point with vertical line and annotation
    if elbow_idx > 0 and elbow_idx < len(values):
        ax.axvline(x=elbow_idx, color='red', linestyle='--', 
                 label=f'Elbow Point at {elbow_idx}')
        
        # Add a clear marker for the elbow point
        ax.plot(elbow_idx, values[elbow_idx-1], 'ro', markersize=8)
        
        # Annotate the value
        ax.annotate(f'Value: {values[elbow_idx-1]:.2f}', 
                  xy=(elbow_idx, values[elbow_idx-1]),
                  xytext=(elbow_idx+5, values[elbow_idx-1]),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1.5))
    
    # Add labels and title
    ax.set_xlabel('Index (Sorted)')
    ax.set_ylabel('Value')
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Elbow Detection')
    
    # Add legend
    ax.legend(loc='upper right')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_factor_heatmap(factor_matrix, factor_name=None, row_labels=None, 
                       col_labels=None, save_path=None, figsize=(12, 8), cmap='viridis'):
    """
    Plot factor matrix as a heatmap.
    
    Parameters:
    -----------
    factor_matrix : array-like
        2D factor matrix to plot
    factor_name : str, optional
        Name of the factor
    row_labels : list, optional
        Labels for rows
    col_labels : list, optional
        Labels for columns
    save_path : str, optional
        Path to save plot
    figsize : tuple, default=(12, 8)
        Figure size
    cmap : str, default='viridis'
        Colormap to use
        
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure
    """
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    sns.heatmap(factor_matrix, annot=False, cmap=cmap, ax=ax, 
               xticklabels=col_labels, yticklabels=row_labels)
    
    # Add labels and title
    ax.set_xlabel('Component')
    ax.set_ylabel('Item')
    title = 'Factor Matrix Heatmap'
    if factor_name:
        title += f' for {factor_name}'
    ax.set_title(title)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_correlation_matrix(corr_matrix, labels=None, save_path=None, 
                          figsize=(10, 8), cmap='coolwarm', annotate=True):
    """
    Plot correlation matrix as a heatmap.
    
    Parameters:
    -----------
    corr_matrix : array-like
        2D correlation matrix to plot
    labels : list, optional
        Labels for rows and columns
    save_path : str, optional
        Path to save plot
    figsize : tuple, default=(10, 8)
        Figure size
    cmap : str, default='coolwarm'
        Colormap to use
    annotate : bool, default=True
        Whether to annotate cells with correlation values
        
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure
    """
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create mask for upper triangle
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    
    # Create heatmap
    sns.heatmap(corr_matrix, mask=mask, cmap=cmap, vmax=1, vmin=-1, center=0,
               square=True, annot=annotate, fmt='.2f', 
               xticklabels=labels, yticklabels=labels, ax=ax)
    
    # Add title
    ax.set_title('Correlation Matrix')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_boxplot(data, x=None, y=None, hue=None, title=None, 
               save_path=None, figsize=(10, 6), orient='v'):
    """
    Create a boxplot of the data.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        Data to plot
    x : str, optional
        Column to use for x-axis
    y : str, optional
        Column to use for y-axis
    hue : str, optional
        Column to use for grouping
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save plot
    figsize : tuple, default=(10, 6)
        Figure size
    orient : str, default='v'
        Orientation of boxplot ('v' for vertical, 'h' for horizontal)
        
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure
    """
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create boxplot
    if orient == 'h':
        sns.boxplot(x=y, y=x, hue=hue, data=data, ax=ax)
    else:
        sns.boxplot(x=x, y=y, hue=hue, data=data, ax=ax)
    
    # Add title
    if title:
        ax.set_title(title)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_network(nodes, edges, node_colors=None, node_sizes=None, 
               node_labels=None, title=None, save_path=None, figsize=(10, 8)):
    """
    Plot a network graph.
    
    Parameters:
    -----------
    nodes : list
        List of node IDs
    edges : list
        List of (source, target) tuples
    node_colors : dict, optional
        Dictionary mapping node IDs to colors
    node_sizes : dict, optional
        Dictionary mapping node IDs to sizes
    node_labels : dict, optional
        Dictionary mapping node IDs to labels
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save plot
    figsize : tuple, default=(10, 8)
        Figure size
        
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx is required for network plotting")
    
    # Create graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Set default node colors if not provided
    if node_colors is None:
        node_colors = {node: 'skyblue' for node in nodes}
    
    # Set default node sizes if not provided
    if node_sizes is None:
        node_sizes = {node: 300 for node in nodes}
    
    # Set default node labels if not provided
    if node_labels is None:
        node_labels = {node: str(node) for node in nodes}
    
    # Create node color and size lists in the same order as nodes
    node_color_list = [node_colors.get(node, 'skyblue') for node in G.nodes()]
    node_size_list = [node_sizes.get(node, 300) for node in G.nodes()]
    
    # Set layout
    pos = nx.spring_layout(G, seed=42)
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, node_color=node_color_list, 
                         node_size=node_size_list, alpha=0.8, ax=ax)
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10, ax=ax)
    
    # Add title
    if title:
        ax.set_title(title)
    
    # Remove axis
    plt.axis('off')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig
