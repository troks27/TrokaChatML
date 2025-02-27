"""
MultiXRank Utilities Module

Functions for running MultiXRank analysis on biological networks.
"""

import os
import yaml
import pandas as pd
import shutil

try:
    import multixrank as mx
    MULTIXRANK_AVAILABLE = True
except ImportError:
    MULTIXRANK_AVAILABLE = False
    print("Warning: multixrank package not available. Some functionality will be limited.")


def create_multixrank_config(config_data, output_path=None):
    """
    Create a MultiXRank configuration file.
    
    Parameters:
    -----------
    config_data : dict
        Dictionary with MultiXRank configuration
    output_path : str, optional
        Path to save the configuration file, if None returns the YAML as string
        
    Returns:
    --------
    str
        Path to saved configuration file or YAML string
    """
    # Create a custom YAML dumper
    class CustomDumper(yaml.Dumper):
        def increase_indent(self, flow=False, indentless=False):
            return super(CustomDumper, self).increase_indent(flow, False)
    
    # Add string presenter for proper formatting
    def str_presenter(dumper, data):
        if len(data) == 2 and data.isdigit():
            return dumper.represent_scalar('tag:yaml.org,2002:str', data, style="'")
        return dumper.represent_scalar('tag:yaml.org,2002:str', data)
    
    yaml.add_representer(str, str_presenter, Dumper=CustomDumper)
    
    # Convert to YAML
    yaml_str = yaml.dump(config_data, Dumper=CustomDumper, default_flow_style=False)
    
    # If output_path is provided, save to file
    if output_path:
        with open(output_path, 'w') as file:
            file.write(yaml_str)
        return output_path
    else:
        return yaml_str


def generate_seeds_file(input_file, output_file):
    """
    Generate a seeds file for MultiXRank from an input file.
    
    Parameters:
    -----------
    input_file : str
        Path to input file containing seed IDs
    output_file : str
        Path to save the seeds file
        
    Returns:
    --------
    str
        Path to the generated seeds file
    """
    # Determine file type and load data
    sep = '\t' if input_file.endswith('.tsv') else ','
    df = pd.read_csv(input_file, sep=sep, header=None)
    
    # Collect unique values from all columns
    unique_values = set()
    for column in df.columns:
        unique_values.update(df[column].unique())
    
    # Sort and save to file
    sorted_values = sorted(unique_values)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Write seeds file
    with open(output_file, 'w') as f:
        for value in sorted_values:
            if pd.notna(value):  # Skip NaN values
                f.write(f"{value}\n")
    
    print(f"Seeds file created at: {output_file}")
    print(f"Total number of unique seeds: {len(sorted_values)}")
    
    return output_file


def run_multixrank_analysis(config_file, working_dir, output_dir=None):
    """
    Run MultiXRank analysis using a configuration file.
    
    Parameters:
    -----------
    config_file : str
        Path to MultiXRank configuration file
    working_dir : str
        Working directory for MultiXRank
    output_dir : str, optional
        Directory to save results, if None uses 'output' subdirectory in working_dir
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with ranking results
    """
    if not MULTIXRANK_AVAILABLE:
        raise ImportError("The multixrank package is required for this function.")
    
    # Initialize MultiXRank
    multixrank_obj = mx.Multixrank(config=config_file, wdir=working_dir)
    
    # Run ranking algorithm
    ranking_df = multixrank_obj.random_walk_rank()
    
    # Export results if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        multixrank_obj.write_ranking(ranking_df, path=output_dir)
    
    return ranking_df


def process_factor(factor_folder, base_dir, config_file):
    """
    Process a factor folder to run MultiXRank analysis.
    
    Parameters:
    -----------
    factor_folder : str
        Path to factor folder
    base_dir : str
        Base directory for MultiXRank
    config_file : str
        Path to MultiXRank configuration file
        
    Returns:
    --------
    dict
        Dictionary with processing results
    """
    if not MULTIXRANK_AVAILABLE:
        raise ImportError("The multixrank package is required for this function.")
    
    # Get factor name from folder
    factor_name = os.path.basename(factor_folder)
    print(f"Processing {factor_name}...")
    
    # Define file paths
    factor_factor_src = os.path.join(factor_folder, 'factor_factor.tsv')
    factor_factor_dst = os.path.join(base_dir, 'multiplex', 'factor', 'factor_factor.tsv')
    factor_gene_src = os.path.join(factor_folder, 'factor_gene.tsv')
    factor_gene_dst = os.path.join(base_dir, 'bipartite', 'factor_gene.tsv')
    
    # Ensure destination directories exist
    os.makedirs(os.path.dirname(factor_factor_dst), exist_ok=True)
    os.makedirs(os.path.dirname(factor_gene_dst), exist_ok=True)
    
    # Copy files to their destinations
    shutil.copy(factor_factor_src, factor_factor_dst)
    shutil.copy(factor_gene_src, factor_gene_dst)
    
    # Generate seeds file
    seeds_file = os.path.join(base_dir, 'seeds.txt')
    generate_seeds_file(factor_factor_src, seeds_file)
    
    # Define output directory
    output_dir = os.path.join(base_dir, f'output_{factor_name}')
    
    # Run MultiXRank analysis
    ranking_df = run_multixrank_analysis(config_file, base_dir, output_dir)
    
    return {
        'factor_name': factor_name,
        'factor_factor_src': factor_factor_src,
        'factor_gene_src': factor_gene_src,
        'seeds_file': seeds_file,
        'output_dir': output_dir,
        'ranking_df': ranking_df
    }


def process_all_factors(factor_graphs_dir, base_dir, config_file):
    """
    Process all factor folders to run MultiXRank analysis.
    
    Parameters:
    -----------
    factor_graphs_dir : str
        Directory containing factor folders
    base_dir : str
        Base directory for MultiXRank
    config_file : str
        Path to MultiXRank configuration file
        
    Returns:
    --------
    dict
        Dictionary with processing results for each factor
    """
    results = {}
    
    # Find all factor folders
    for folder in os.listdir(factor_graphs_dir):
        if folder.startswith('Factor'):
            factor_folder = os.path.join(factor_graphs_dir, folder)
            try:
                result = process_factor(factor_folder, base_dir, config_file)
                results[folder] = result
            except Exception as e:
                results[folder] = {'error': str(e)}
    
    return results
