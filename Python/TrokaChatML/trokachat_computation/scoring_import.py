import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from .scoring_pipeline_utils import (
    load_h5_file, save_h5_file, process_single_file, extract_sample_name, extract_sample_name_xlsx,
    preprocess_counts_file, create_lr_dict, find_matching_lr_pairs, group_results, load_h5_file_groups, save_h5_file_group, process_group
)

class TrokaChat_import:
    def __init__(self, input_path: str, file_prefixes: list, output_path: str, pct_min: float):
        self.input_path = input_path
        self.file_prefixes = file_prefixes
        self.output_path = output_path
        self.pct_min = pct_min
        self.frame = None
        self.file_data = {}  # To store data corresponding to each file

    def load_h5_file(self, filepath):
        print(f"Loading HDF5 file: {filepath}")
        data_dict = {}
        with h5py.File(filepath, 'r') as f:
            for group in f.keys():
                for dset in f[group].keys():
                    ds_data = f[group][dset][:]
                    df = pd.DataFrame(ds_data)
                    df['group'] = group
                    df['gene'] = df['gene'].str.decode("utf-8")
                    if group not in data_dict:
                        data_dict[group] = []
                    data_dict[group].append(df)
        for group in data_dict:
            data_dict[group] = pd.concat(data_dict[group], ignore_index=True)
        return data_dict

    def load_data(self):
        frames = []
        for prefix in self.file_prefixes:
            file_path = f"{self.input_path}/{prefix}"
            if prefix.endswith('.csv'):
                print(f"Loading file: {prefix}")
                df = pd.read_csv(file_path)
                self.file_data[prefix] = df
            elif prefix.endswith('.h5'):
                self.file_data[prefix] = self.load_h5_file(file_path)
            
            # Adjust condition assignment for null distribution files:
            if prefix.startswith('nulldist'):
                # Use the second part of the filename to name the sample,
                # so that "nulldist_NG Cre minus_vs_NG Cre minus.h5" becomes "nulldistavgs_NG Cre minus"
                condition = 'nulldistavgs_' + prefix.split('_')[1]
            else:
                condition = prefix.split('_')[0]
            
            # Assign the sample column and collect frames
            if isinstance(self.file_data[prefix], dict):
                for key in self.file_data[prefix]:
                    self.file_data[prefix][key]['sample'] = condition
                    frames.append(self.file_data[prefix][key])
            else:
                self.file_data[prefix]['sample'] = condition
                frames.append(self.file_data[prefix])
        self.frame = pd.concat(frames, ignore_index=True)

    def filter_transform_data(self):
        self.frame = self.frame[self.frame["pct.1"] > self.pct_min]
        self.frame['gene2'] = self.frame['gene']
        self.frame['avg_log2FC'] = (
            (self.frame['avg_log2FC'] - self.frame['avg_log2FC'].min()) *
            (self.frame['avg_log2FC'].max() + abs(self.frame['avg_log2FC'].min()) + 1 -
            (self.frame['avg_log2FC'].min() + abs(self.frame['avg_log2FC'].min()) + 1)) /
            (self.frame['avg_log2FC'].max() - self.frame['avg_log2FC'].min()) +
            (self.frame['avg_log2FC'].min() + abs(self.frame['avg_log2FC'].min()) + 1)
        )
        self.frame.drop(columns=['gene'], inplace=True)
        self.frame.rename(columns={'gene2': 'gene'}, inplace=True)

    def plot_data(self):
        plt.hist(self.frame["avg_log2FC"], bins=25, density=True, alpha=0.6, color='g')
        plt.show()

    def describe_data(self):
        print(self.frame.describe())

    def export_data(self):
        for condition, df in self.frame.groupby('sample'):
            if 'nulldistavgs_' in condition:
                self.save_h5_file(condition, df)
            else:
                df.to_csv(f"{self.output_path}/{condition} DEGs.csv", index=False)

    def save_h5_file(self, condition, df):
        filepath = f"{self.output_path}/{condition}.h5"
        group = condition.split('_')[1]
        print(f"Saving HDF5 file: {filepath}")
        with h5py.File(filepath, 'w') as hf:
            for i, (name, group_df) in enumerate(df.groupby('group'), start=1):
                group_name = group_df['group'].unique()
                group_df.drop(group_df.columns[group_df.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)
                group_df = group_df.dropna(axis=1)
                group_df = group_df[['p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'gene', 'cluster']].copy()
                
                # Convert 'gene' column to fixed-length string
                group_df['gene'] = group_df['gene'].apply(lambda x: x.encode('utf-8'))

                str_dt = h5py.string_dtype(encoding='utf-8')
                
                # Define the dtype for the structured array
                dt = np.dtype([
                    ('p_val', '<f8'),
                    ('avg_log2FC', '<f8'),
                    ('pct.1', '<f8'),
                    ('pct.2', '<f8'),
                    ('p_val_adj', '<f8'),
                    ('gene', str_dt),  # Fixed-length string of 9 characters
                    ('cluster', '<f8')
                ])
                
                # Convert the DataFrame to a structured array
                data_array = np.array([tuple(row) for row in group_df.to_records(index=False)], dtype=dt)
                
                # Create the group and dataset
                group = hf.create_group(group_name[0])
                group.create_dataset('data', data=data_array)

    def import_data(self):
        self.load_data()
        self.filter_transform_data()
        self.plot_data()
        self.describe_data()
        self.export_data()
