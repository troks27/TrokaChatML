# data_processing.py
import h5py
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

def list_h5_groups(filepath):
    with h5py.File(filepath, 'r') as f:
        groups = list(f.keys())
    return groups

def load_h5_batches(filepath, group_names):
    data_dict = {}
    with h5py.File(filepath, 'r') as f:
        for group in group_names:
            dset = f[group]['data'][:]
            df = pd.DataFrame(dset)
            df['Ligand'] = df['Ligand'].str.decode("utf-8")
            df['Receptor'] = df['Receptor'].str.decode("utf-8")
            df['Pathway Label'] = df['Pathway Label'].str.decode("utf-8")
            df['group'] = group
            if group not in data_dict:
                data_dict[group] = []
            data_dict[group].append(df)
        for group in data_dict:
            data_dict[group] = pd.concat(data_dict[group], ignore_index=True)
    return data_dict

def preprocess_datasets(data_dict):
    for key in data_dict:
        df = data_dict[key]
        df = df.drop(columns=['group'])
        epsilon = 1e-6
        df['logTF Communication Score'] = np.log(df['Communication Score'] + epsilon)
        df['source_target'] = df['Source'].astype(str) + '_' + df['Target'].astype(str)
        df['ligand_receptor'] = df['Ligand'].astype(str) + '_' + df['Receptor'].astype(str)
        df['source_ligand'] = df['Source'].astype(str) + '_' + df['Ligand'].astype(str)
        df['target_receptor'] = df['Target'].astype(str) + '_' + df['Receptor'].astype(str)
        categorical_columns = ['Ligand', 'Receptor', 'Source', 'Target', 'Pathway Label',
                               'source_target', 'ligand_receptor', 'source_ligand', 'target_receptor']
        for col in categorical_columns:
            df[col] = df[col].astype('category')
        data_dict[key] = df
    return data_dict

def preprocess_and_combine_datasets(data_dict):
    combined_data = pd.DataFrame()
    for df in data_dict.values():
        combined_data = pd.concat([combined_data, df], ignore_index=True)
    categorical_columns = ['Ligand', 'Receptor', 'Source', 'Target', 'Pathway Label',
                           'source_target', 'ligand_receptor', 'source_ligand', 'target_receptor']
    cols = categorical_columns
    combined_data.is_copy = False
    combined_data[cols] = combined_data[cols].astype('category')
    return combined_data

def preprocess_test_data(test_groups, filename, batch_size, n_jobs=9):
    def load_and_preprocess_wrapper(groups):
        data_dict = load_h5_batches(filename, groups)
        data_dict = preprocess_datasets(data_dict)
        return preprocess_and_combine_datasets(data_dict)

    with Parallel(n_jobs=n_jobs) as parallel:
        combined_data_batches = parallel(delayed(load_and_preprocess_wrapper)(test_groups[start:start + batch_size]) for start in range(0, len(test_groups), batch_size))
        combined_data = pd.concat(combined_data_batches, ignore_index=True)

    excluded_columns = ['logTF Communication Score', 'Communication Score', 'Uniqueness Score']
    feature_columns = [col for col in combined_data.columns if col not in excluded_columns]
    cat_indices = [i for i, col in enumerate(feature_columns) if combined_data[col].dtype.name == 'category']
    X_test = combined_data[feature_columns]
    y_test = combined_data['logTF Communication Score']
    
    cols = feature_columns
    X_test.is_copy = False
    X_test[cols] = X_test[cols].astype('category')
    return X_test, y_test

def preprocess_real_data(df):
    epsilon = 1e-6
    df['logTF Communication Score'] = np.log(df['Communication Score'] + epsilon)
    df['source_target'] = df['Source'].astype(str) + '_' + df['Target'].astype(str)
    df['ligand_receptor'] = df['Ligand'].astype(str) + '_' + df['Receptor'].astype(str)
    df['source_ligand'] = df['Source'].astype(str) + '_' + df['Ligand'].astype(str)
    df['target_receptor'] = df['Target'].astype(str) + '_' + df['Receptor'].astype(str)
    categorical_columns = ['Ligand', 'Receptor', 'Source', 'Target', 'Pathway Label',
                           'source_target', 'ligand_receptor', 'source_ligand', 'target_receptor']
    for col in categorical_columns:
        df[col] = df[col].astype('category')
    excluded_columns = ['logTF Communication Score', 'Communication Score', 'Uniqueness Score']
    feature_columns = [col for col in df.columns if col not in excluded_columns]
    X_real = df[feature_columns]
    y_real = df['logTF Communication Score']
    return X_real, y_real