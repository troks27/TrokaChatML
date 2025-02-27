"""
Data Processing Module

Functions for loading, processing, and filtering cell communication data.
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit


def process_pct_filter(main_df, reference_df, pvalue_cutoff=0.05):
    """
    Process and filter DataFrame based on percentage cutoff and p-values.
    
    Parameters:
    -----------
    main_df : pandas.DataFrame
        The main DataFrame to process
    reference_df : pandas.DataFrame
        Reference DataFrame containing p-values
    pvalue_cutoff : float, default=0.05
        P-value threshold for filtering
        
    Returns:
    --------
    pandas.DataFrame
        Processed DataFrame with additional columns and filters
    """
    def process_csv(main_df):
        # Create the ligand_pct_adj_pval column (example: fill it with NaN or some default value)
        main_df['ligand_pct_adj_pval'] = float('nan')
        
        # Find the row with the most receptor subunits
        def count_subunits(receptor):
            return len(receptor.split('_'))
        
        main_df['receptor_subunit_count'] = main_df['Receptor'].apply(count_subunits)
        max_subunits = main_df['receptor_subunit_count'].max()
        
        # Create new columns for each receptor subunit based on the maximum number of subunits
        for i in range(1, max_subunits + 1):
            main_df[f'receptor_{i}_pct_adj_pval'] = pd.Series([float('nan')] * len(main_df), dtype='object')
        
        # Split receptor subunits and assign values to new columns
        for index, row in main_df.iterrows():
            subunits = row['Receptor'].split('_')
            for i, subunit in enumerate(subunits):
                main_df.at[index, f'receptor_{i+1}_pct_adj_pval'] = subunit
        
        # Drop the temporary column used for counting subunits
        main_df.drop(columns=['receptor_subunit_count'], inplace=True)
        
        return main_df
    
    def update_pvalues(main_df, reference_df):
        # Update ligand_pct_adj_pval based on matching Ligand and Source with gene and cluster
        for idx, row in main_df.iterrows():
            gene = row['Ligand']
            cluster = row['Source']
            match = reference_df[(reference_df['gene'] == gene) & (reference_df['cluster'] == cluster)]
            if not match.empty:
                main_df.at[idx, 'ligand_pct_adj_pval'] = match['adjusted_p_value'].values[0]
    
        # Update receptor_n_pct_adj_pval based on matching Receptor subunits and Target with gene and cluster
        for idx, row in main_df.iterrows():
            receptor_subunits = row['Receptor'].split('_')
            cluster = row['Target']
            
            for i, subunit in enumerate(receptor_subunits):
                match = reference_df[(reference_df['gene'] == subunit) & (reference_df['cluster'] == cluster)]
                if not match.empty:
                    main_df.at[idx, f'receptor_{i+1}_pct_adj_pval'] = match['adjusted_p_value'].values[0]
        return main_df
    
    def apply_pvalue_filter(main_df, pvalue_cutoff=pvalue_cutoff):
        # Ensure the ligand_pct_adj_pval column is numeric and round to 4 decimal places
        main_df['ligand_pct_adj_pval'] = pd.to_numeric(main_df['ligand_pct_adj_pval'], errors='coerce').round(4)
        
        # Identify the receptor columns based on the maximum number of receptor subunits
        receptor_cols = [col for col in main_df.columns if col.startswith('receptor_') and col.endswith('_pct_adj_pval')]
        
        # Ensure all receptor columns are numeric and round to 4 decimal places
        for col in receptor_cols:
            main_df[col] = pd.to_numeric(main_df[col], errors='coerce').round(4)
        
        # Function to check if all relevant p-values are below the cutoff
        def check_pvalues(row):
            pvalues = [row['ligand_pct_adj_pval']] + [row[col] for col in receptor_cols]
            pvalues = [p for p in pvalues if pd.notna(p)]  # Remove NaN values from consideration
            return 1 if all(p <= pvalue_cutoff for p in pvalues) else 0
        
        # Apply the function to each row to create the pct_filter column
        main_df['pct_filter'] = main_df.apply(check_pvalues, axis=1)
        
        return main_df
    
    main_df_1 = process_csv(main_df)
    main_df_1 = update_pvalues(main_df_1, reference_df)
    main_df_1 = apply_pvalue_filter(main_df_1, pvalue_cutoff)

    return main_df_1


def train_set_create(exp_df, control_df):
    """
    Create training dataset by filtering out communications that are only in the control condition.
    
    Parameters:
    -----------
    exp_df : pandas.DataFrame
        Experimental condition DataFrame
    control_df : pandas.DataFrame
        Control condition DataFrame
        
    Returns:
    --------
    pandas.DataFrame
        Filtered control DataFrame for training
    """
    # Identify communications that are only in the control condition
    control_only = control_df[~control_df[['Ligand', 'Receptor', 'Source', 'Target']].apply(tuple, axis=1).isin(
        exp_df[['Ligand', 'Receptor', 'Source', 'Target']].apply(tuple, axis=1))]

    # Exclude control-only communications from the control DataFrame
    control_df_filtered = control_df[~control_df.index.isin(control_only.index)]

    # Print the shape of the resulting DataFrame
    print(f"Final DataFrame shape: {control_df_filtered.shape[0]} rows, {control_df_filtered.shape[1]} columns")
    
    return control_df_filtered


def preprocess_datasets(df):
    """
    Preprocess datasets by dropping unnecessary columns and converting categorical columns.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to preprocess
        
    Returns:
    --------
    pandas.DataFrame
        Preprocessed DataFrame
    """
    df = df.drop(['Predicted Communication Score', 'Residual', 'z_score', 'p_value'], axis=1)
    
    categorical_columns = ['Ligand', 'Receptor', 'Source', 'Target', 'Pathway Label',
                          'source_target', 'ligand_receptor', 'source_ligand', 'target_receptor']
    
    for col in categorical_columns:
        df[col] = df[col].astype('category')
    
    return df


def test_set_create(exp_df, control_df):
    """
    Create test dataset by identifying rows based on specific conditions.
    
    Parameters:
    -----------
    exp_df : pandas.DataFrame
        Experimental condition DataFrame
    control_df : pandas.DataFrame
        Control condition DataFrame
        
    Returns:
    --------
    pandas.DataFrame
        Final DataFrame for testing
    """
    # Step 1: Merge the DataFrames to find common and unique rows
    merged_df = exp_df.merge(
        control_df, 
        on=['Ligand', 'Receptor', 'Source', 'Target'], 
        suffixes=('_exp', '_control'), 
        how='outer',
        indicator=True
    )
    
    # Step 2: Identify rows based on the adjusted_p_value conditions
    # Condition 1: Common rows with adjusted_p_value < 0.05 in either DataFrame
    common_condition = (
        (merged_df['_merge'] == 'both') & 
        (
            ((merged_df['adjusted_p_value_exp'] < 0.05) & (merged_df['pct_filter_exp'] == 1)) | 
            ((merged_df['adjusted_p_value_control'] < 0.05) & (merged_df['pct_filter_control'] == 1))
        )
    )
    
    # Create a new column to indicate if both, or one or the other condition has 1 as the value
    merged_df['condition_status'] = 'None'
    
    # Both conditions are met
    both_conditions = (
        (merged_df['_merge'] == 'both') & 
        ((merged_df['adjusted_p_value_exp'] < 0.05) & (merged_df['pct_filter_exp'] == 1)) & 
        ((merged_df['adjusted_p_value_control'] < 0.05) & (merged_df['pct_filter_control'] == 1))
    )
    merged_df.loc[both_conditions, 'condition_status'] = 'Both'
    
    # Only the experimental condition is met
    exp_only_condition = (
        (merged_df['_merge'] == 'both') & 
        ((merged_df['adjusted_p_value_exp'] < 0.05) & (merged_df['pct_filter_exp'] == 1)) & 
        ~((merged_df['adjusted_p_value_control'] < 0.05) & (merged_df['pct_filter_control'] == 1))
    )
    merged_df.loc[exp_only_condition, 'condition_status'] = 'Exp Only'
    
    # Only the control condition is met
    control_only_condition = (
        (merged_df['_merge'] == 'both') & 
        ~((merged_df['adjusted_p_value_exp'] < 0.05) & (merged_df['pct_filter_exp'] == 1)) & 
        ((merged_df['adjusted_p_value_control'] < 0.05) & (merged_df['pct_filter_control'] == 1))
    )
    merged_df.loc[control_only_condition, 'condition_status'] = 'Control Only'

    # Condition 2: Rows only in exp_df with adjusted_p_value < 0.05
    exp_only_rows = (
        (merged_df['_merge'] == 'left_only') & 
        (merged_df['adjusted_p_value_exp'] < 0.05) & 
        (merged_df['pct_filter_exp'] == 1)
    )
    
    # Condition 3: Rows only in control_df with adjusted_p_value < 0.05
    control_only_rows = (
        (merged_df['_merge'] == 'right_only') & 
        (merged_df['adjusted_p_value_control'] < 0.05) & 
        (merged_df['pct_filter_control'] == 1)
    )
    
    # Collect the rows to keep based on conditions
    rows_to_keep = merged_df[common_condition | exp_only_rows | control_only_rows]
    
    # Ensure condition_status is retained
    exp_rows_to_keep = exp_df.merge(
        rows_to_keep[['Ligand', 'Receptor', 'Source', 'Target', 'condition_status']],
        on=['Ligand', 'Receptor', 'Source', 'Target'],
        how='inner'
    )
    
    control_rows_to_keep = control_df.merge(
        rows_to_keep[['Ligand', 'Receptor', 'Source', 'Target', 'condition_status']],
        on=['Ligand', 'Receptor', 'Source', 'Target'],
        how='inner'
    )
    
    # Step 4: Concatenate the extracted rows
    final_df = pd.concat([exp_rows_to_keep, control_rows_to_keep], ignore_index=True)
    
    # Step 5: Remove duplicate rows
    final_df = final_df.drop_duplicates(subset=['Ligand', 'Receptor', 'Source', 'Target'])
    
    # Step 6: Update Condition column for rows found in both exp and control
    common_rows = final_df[['Ligand', 'Receptor', 'Source', 'Target']].apply(tuple, axis=1).isin(
        merged_df[common_condition][['Ligand', 'Receptor', 'Source', 'Target']].apply(tuple, axis=1)
    )
    
    final_df.loc[common_rows, 'Condition'] = 'both'
    
    print(f"Final DataFrame shape: {final_df.shape[0]} rows, {final_df.shape[1]} columns")
    return final_df


def stratified_split_handle_rare_classes(df, target_column, test_size=0.1, threshold=2, random_state=42):
    """
    Perform a stratified split of the DataFrame, handling rare classes by combining them into an 'other' category.

    Parameters:
    -----------
    df : pandas.DataFrame
        The input DataFrame
    target_column : str
        The name of the target column to stratify by
    test_size : float, default=0.1
        The proportion of the dataset to include in the test split
    threshold : int, default=2
        The minimum number of instances a class should have to not be combined into 'other'
    random_state : int, default=42
        Random seed for reproducibility

    Returns:
    --------
    tuple
        (train_df, test_df) - Training and test DataFrames
    """
    def handle_single_instance_classes(df, target_column, threshold=2):
        # Combine rare classes into 'other'
        class_counts = df[target_column].value_counts()
        df[target_column] = df[target_column].apply(
            lambda x: x if class_counts[x] >= threshold else 'other'
        )
        return df

    # Handle single-instance classes
    df = handle_single_instance_classes(df, target_column, threshold)

    # Proceed with the stratified split
    X = df.drop(columns=[target_column])
    y = df[target_column]

    # Separate 'other' instances
    other_instances = df[df[target_column] == 'other']
    df = df[df[target_column] != 'other']

    # Update X and y after removing 'other' instances
    X = df.drop(columns=[target_column])
    y = df[target_column]

    # Perform the stratified split
    splitter = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)

    for train_index, test_index in splitter.split(X, y):
        train_df = df.iloc[train_index]
        test_df = df.iloc[test_index]

    # Split the 'other' class instances between train and test sets
    num_other_instances = len(other_instances)

    if num_other_instances > 0:
        if num_other_instances == 1:
            # If there's only one 'other' instance, add it to the training set
            train_df = pd.concat([train_df, other_instances])
        else:
            # Calculate the split
            half_split = num_other_instances // 2
            train_other = other_instances.iloc[:half_split + (num_other_instances % 2)]
            test_other = other_instances.iloc[half_split + (num_other_instances % 2):]

            # Add to train and test sets
            train_df = pd.concat([train_df, train_other])
            if not test_other.empty:
                test_df = pd.concat([test_df, test_other])

    # Reset the index for train and test sets
    train_df = train_df.reset_index(drop=True)
    test_df = test_df.reset_index(drop=True)

    return train_df, test_df


def load_data(control_path, exp_path, control_pct_path, exp_pct_path, 
              control_condition="NG", exp_condition="DIAB"):
    """
    Load data from CSV files and add condition labels.
    
    Parameters:
    -----------
    control_path : str
        Path to control prediction CSV file
    exp_path : str
        Path to experimental prediction CSV file
    control_pct_path : str
        Path to control results CSV file
    exp_pct_path : str
        Path to experimental results CSV file
    control_condition : str, default="NG"
        Label for control condition
    exp_condition : str, default="DIAB"
        Label for experimental condition
        
    Returns:
    --------
    tuple
        (control_df, exp_df, control_pct, exp_pct) - Loaded DataFrames
    """
    # Load prediction data
    control_df = pd.read_csv(control_path)
    control_df['Condition'] = control_condition

    exp_df = pd.read_csv(exp_path)
    exp_df['Condition'] = exp_condition

    # Load PCT data
    control_pct = pd.read_csv(control_pct_path)
    control_pct['Condition'] = control_condition
    control_pct['adjusted_p_value'].fillna(1, inplace=True)

    exp_pct = pd.read_csv(exp_pct_path)
    exp_pct['Condition'] = exp_condition
    exp_pct['adjusted_p_value'].fillna(1, inplace=True)
    
    return control_df, exp_df, control_pct, exp_pct


def filter_with_residuals(exp_final, test_residuals, rmse_threshold=0.81, output_path=None):
    """
    Filter data based on residuals and save to CSV if output_path is provided.
    
    Parameters:
    -----------
    exp_final : pandas.DataFrame
        DataFrame with experimental data
    test_residuals : pandas.Series or numpy.ndarray
        Residuals from model prediction
    rmse_threshold : float, default=0.81
        Threshold for filtering based on RMSE
    output_path : str, optional
        Path to save filtered data, if None doesn't save
        
    Returns:
    --------
    pandas.DataFrame
        Filtered DataFrame
    """
    # Add the residuals back to the experimental data
    exp_final['Residuals'] = test_residuals  # Ensure the index matches

    # Apply the filtering conditions
    filtered_exp_final = exp_final[
        ~(
            ((exp_final['Condition'] == 'LS') & (exp_final['Residuals'] < -rmse_threshold)) |
            ((exp_final['Condition'] == 'NL') & (exp_final['Residuals'] > rmse_threshold)) |
            ((exp_final['condition_status'] == 'Control Only') & (exp_final['Residuals'] > rmse_threshold)) |
            ((exp_final['condition_status'] == 'Exp Only') & (exp_final['Residuals'] < -rmse_threshold))
        )
    ]

    # Remove the specified columns
    columns_to_remove = ['ligand_pct_adj_pval', 'pct_filter', 'condition_status'] + [
        col for col in filtered_exp_final.columns if col.startswith('receptor_')
    ]
    filtered_exp_final = filtered_exp_final.drop(columns=columns_to_remove)

    # Save the updated DataFrame to a CSV file if output_path is provided
    if output_path:
        filtered_exp_final.to_csv(output_path, index=False)
        print(f"Residuals added, filtered, columns removed, and data saved to '{output_path}'")
    
    return filtered_exp_final
