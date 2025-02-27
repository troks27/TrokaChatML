"""
Gene Mapping Module

Functions for mapping genes between species using BioMart.
"""

import os
import pandas as pd
import itertools
from pybiomart import Server


def get_mouse_mapping(human_genes):
    """
    Get mappings from human genes to mouse homologs using BioMart.
    
    Parameters:
    -----------
    human_genes : list
        List of human gene symbols to map
        
    Returns:
    --------
    dict
        Dictionary mapping human genes to mouse genes
    """
    print("Querying Ensembl via pybiomart for human-to-mouse gene mappings...")
    server = Server(host='http://www.ensembl.org')
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    
    # Query 1: Get human gene IDs and symbols
    df_human = dataset.query(
        attributes=['ensembl_gene_id', 'external_gene_name'],
        use_attr_names=True
    )
    
    # Query 2: Get human gene IDs and their mouse homolog associated gene names
    df_mouse = dataset.query(
        attributes=['ensembl_gene_id', 'mmusculus_homolog_associated_gene_name'],
        use_attr_names=True
    )
    
    # Merge the two dataframes on 'ensembl_gene_id'
    merged = pd.merge(df_human, df_mouse, on='ensembl_gene_id', how='left')
    # Rename columns for clarity
    merged.rename(columns={
        'external_gene_name': 'human_gene',
        'mmusculus_homolog_associated_gene_name': 'mouse_gene'
    }, inplace=True)
    
    # Build the mapping dictionary
    # We assume 'human_genes' is already uppercase from collect_gene_names
    mapping_dict = {gene: set() for gene in human_genes}
    for _, row in merged.iterrows():
        # Normalize the human_gene to uppercase
        if pd.isna(row['human_gene']):
            continue
        human_gene = row['human_gene'].strip().upper()
        mouse_gene_raw = row['mouse_gene']

        # If this 'human_gene' is in our target set
        if human_gene in mapping_dict:
            # If no mouse homolog is provided, default to the human gene
            if pd.isna(mouse_gene_raw) or mouse_gene_raw.strip() == "":
                mapping_dict[human_gene].add(human_gene)
            else:
                # Sometimes multiple homologs are returned as a comma-separated string
                for mg in mouse_gene_raw.split(','):
                    mg_stripped = mg.strip()
                    if mg_stripped:  # Only add non-empty strings
                        mapping_dict[human_gene].add(mg_stripped)
    
    # Ensure that every gene in our list has at least one mapping
    for gene in mapping_dict:
        if len(mapping_dict[gene]) == 0:
            mapping_dict[gene].add(gene)
    
    # Convert sets to lists
    mapping_dict = {gene: list(genes) for gene, genes in mapping_dict.items()}
    print(f"Mapping complete for {len(mapping_dict)} genes.")
    return mapping_dict


def collect_gene_names(file_paths, gene_columns):
    """
    Collect unique gene names from files.
    
    Parameters:
    -----------
    file_paths : list
        List of file paths to extract genes from
    gene_columns : list
        List of column indices (0-indexed) containing gene names
        
    Returns:
    --------
    list
        List of unique gene names
    """
    gene_set = set()
    for f in file_paths:
        df = pd.read_csv(f, sep='\t', header=None, dtype=str)
        # Drop rows that have NaN in any gene column
        df.dropna(subset=gene_columns, inplace=True)
        # Strip whitespace and convert to uppercase
        for col in gene_columns:
            df[col] = df[col].str.strip().str.upper()
        # Add to our global set
        for col in gene_columns:
            gene_set.update(df[col].unique())
    return list(gene_set)


def transform_gene_file(input_path, output_path, gene_columns, mapping_dict):
    """
    Transform a file by replacing gene names with their mappings.
    
    Parameters:
    -----------
    input_path : str
        Path to input file
    output_path : str
        Path to save transformed file
    gene_columns : list
        List of column indices (0-indexed) containing gene names
    mapping_dict : dict
        Dictionary mapping genes to their homologs
        
    Returns:
    --------
    str
        Path to the transformed file
    """
    df = pd.read_csv(input_path, sep='\t', header=None, dtype=str)
    # Drop rows that have NaN in any gene column
    df.dropna(subset=gene_columns, inplace=True)
    
    # Normalize each gene column: strip + uppercase
    for col in gene_columns:
        df[col] = df[col].str.strip().str.upper()
    
    unique_rows = set()
    for _, row in df.iterrows():
        # For each gene column, get the list of mapped gene symbols
        mapped_lists = []
        for col in gene_columns:
            gene = row[col]
            # Retrieve the mapped mouse genes from our dictionary; fallback is the same gene
            mapped = mapping_dict.get(gene, [gene])
            mapped_lists.append(mapped)
        
        # Cartesian product across columns that need replacement
        for combination in itertools.product(*mapped_lists):
            new_row = list(row)  # row is a Series; convert to list
            for i, col_index in enumerate(gene_columns):
                new_row[col_index] = combination[i]
            # Add as a tuple to avoid duplicates
            unique_rows.add(tuple(new_row))
    
    # Convert our unique set of rows back into a DataFrame
    new_df = pd.DataFrame(list(unique_rows))
    # Optional: final safety measure to remove duplicates
    new_df.drop_duplicates(inplace=True)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    new_df.to_csv(output_path, sep='\t', header=False, index=False)
    
    return output_path


def map_all_gene_files(multiplex_gene_files, bipartite_gene_files, output_dir, 
                       multiplex_gene_columns, bipartite_gene_columns):
    """
    Map all gene files from human to mouse.
    
    Parameters:
    -----------
    multiplex_gene_files : list
        List of multiplex gene file paths
    bipartite_gene_files : list
        List of bipartite gene file paths
    output_dir : str
        Directory to save mapped files
    multiplex_gene_columns : list
        List of column indices for genes in multiplex files
    bipartite_gene_columns : list
        List of column indices for genes in bipartite files
        
    Returns:
    --------
    dict
        Dictionary with mapping results
    """
    from .utils import ensure_directory_exists
    
    # 1. Collect all unique gene names from both sets of files (convert to uppercase)
    all_gene_names = set(collect_gene_names(multiplex_gene_files, multiplex_gene_columns))
    all_gene_names.update(collect_gene_names(bipartite_gene_files, bipartite_gene_columns))
    all_gene_names = list(all_gene_names)
    print("Total unique human genes to map:", len(all_gene_names))

    # 2. Get the mapping dictionary using pybiomart
    mapping_dict = get_mouse_mapping(all_gene_names)

    # 3. Process each file and write to the corresponding output path
    results = {'multiplex': {}, 'bipartite': {}}
    
    #    (a) Multiplex gene files
    for in_path in multiplex_gene_files:
        # Determine output path (mirroring input folder structure)
        rel_path = os.path.relpath(in_path, os.path.dirname(output_dir))
        out_path = os.path.join(output_dir, rel_path)
        ensure_directory_exists(os.path.dirname(out_path))
        
        # Transform the file
        transformed_path = transform_gene_file(in_path, out_path, multiplex_gene_columns, mapping_dict)
        results['multiplex'][in_path] = transformed_path

    #    (b) Bipartite gene files
    for in_path in bipartite_gene_files:
        rel_path = os.path.relpath(in_path, os.path.dirname(output_dir))
        out_path = os.path.join(output_dir, rel_path)
        ensure_directory_exists(os.path.dirname(out_path))
        
        # Transform only the specified columns
        transformed_path = transform_gene_file(in_path, out_path, bipartite_gene_columns, mapping_dict)
        results['bipartite'][in_path] = transformed_path

    return results
