import pandas as pd
import numpy as np
import h5py
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice
from .scoring_functions import calculate_scores_parallel

def load_h5_file(filepath):
    print(f"Loading HDF5 file: {filepath}")
    data_dict = {}
    with h5py.File(filepath, 'r') as f:
        for group in f.keys():
            #print(f"Found group: {group}")
            data_dict[group] = []
            for dset in f[group].keys():
                #print(f"Found dataset: {dset} in group: {group}")
                ds_data = f[group][dset][:]
                df = pd.DataFrame(ds_data)
                df['group'] = group
                if 'gene' in df.columns:
                    df['gene'] = df['gene'].str.decode("utf-8")
                data_dict[group].append(df)
    for group in data_dict:
        data_dict[group] = pd.concat(data_dict[group], ignore_index=True)
    
    
    return data_dict

def load_h5_file_groups(filepath):
    print(f"Loading HDF5 file: {filepath}")
    with h5py.File(filepath, 'r') as f:
        for group in f.keys():
            data_dict = []
            for dset in f[group].keys():
                ds_data = f[group][dset][:]
                df = pd.DataFrame(ds_data)
                df['group'] = group
                if 'gene' in df.columns:
                    df['gene'] = df['gene'].str.decode("utf-8")
                data_dict.append(df)
            data_dict = pd.concat(data_dict, ignore_index=True)
            yield group, data_dict


def save_h5_file(filepath, permutations_dict):
    with h5py.File(filepath, 'w') as hf:
        for perm_name, df in permutations_dict.items():
            group = hf.create_group(perm_name)
            
            # Define variable-length string data type
            str_dt = h5py.string_dtype(encoding='utf-8')
            
            # Create the dataset with variable-length string dtype for specific columns
            dt = np.dtype([
                ('Ligand', str_dt),
                ('Receptor', str_dt),
                ('Source', '<f8'),
                ('Target', '<f8'),
                ('Communication Score', '<f8'),
                ('Uniqueness Score', '<f8'),
                ('Pathway Label', str_dt)
            ])
            
            # Convert DataFrame to a structured array
            data_array = np.array([tuple(row) for row in df.to_records(index=False)], dtype=dt)
            
            # Create dataset within the group
            group.create_dataset('data', data=data_array)

def save_h5_file_group(output_path, group_name, df):
    with h5py.File(output_path, 'a') as hf:
        group = hf.create_group(group_name)
        str_dt = h5py.string_dtype(encoding='utf-8')
        dt = np.dtype([
            ('Ligand', str_dt),
            ('Receptor', str_dt),
            ('Source', '<f8'),
            ('Target', '<f8'),
            ('Communication Score', '<f8'),
            ('Uniqueness Score', '<f8'),
            ('Pathway Label', str_dt)
        ])
        data_array = np.array([tuple(row) for row in df.to_records(index=False)], dtype=dt)
        group.create_dataset('data', data=data_array)



def process_group(group_name, df, lr_db_df, general_filepath, output_path, min_clus_percent, sample_name_from_R, counts_file, import_folder):
    sample_name = extract_sample_name(group_name)
    return process_single_file(group_name, df, sample_name, lr_db_df, general_filepath, output_path, min_clus_percent, sample_name_from_R, counts_file, import_folder)





def process_single_file(group_name, df, sample_name, lr_db_df, general_filepath, output_path, min_clus_percent, sample_name_from_R, counts_file, import_folder):
    #print(f"Processing group: {group_name}, sample name: {sample_name}")

    deg_data = df.loc[df['avg_log2FC'] > 0]
    deg_data['pct.2'] = deg_data['pct.2'].replace(0, 0.001)
    deg_data = deg_data.reset_index(drop=True)

    new_counts_list, cluster_percent, counts_filtered = preprocess_counts_file(min_clus_percent, general_filepath, import_folder, counts_file)
    lr_dict = create_lr_dict(lr_db_df)

    matching_pairs_df = find_matching_lr_pairs(deg_data, lr_dict, cluster_percent, sample_name, counts_filtered, general_filepath)

    #print(f"matching_pairs_df columns: {matching_pairs_df.columns}")

    if 'Source' not in matching_pairs_df.columns or 'Target' not in matching_pairs_df.columns:
        raise KeyError("Columns 'Source' or 'Target' are missing from matching_pairs_df")

    matching_pairs_df['Source'] = matching_pairs_df['Source'].astype(int).astype(str)
    matching_pairs_df['Target'] = matching_pairs_df['Target'].astype(int).astype(str)

    new_counts_list = [str(int(x)) for x in new_counts_list if x.isdigit()]
    matching_pairs_df.set_index('Target', inplace=True)
    matching_pairs_df = matching_pairs_df.loc[matching_pairs_df.index.isin(new_counts_list)].reset_index()
    matching_pairs_df.set_index('Source', inplace=True)
    matching_pairs_df = matching_pairs_df.loc[matching_pairs_df.index.isin(new_counts_list)].reset_index()
    matching_pairs_df = matching_pairs_df[['Ligand', 'Receptor', 'Source', 'Target', 'Communication Score', 'Uniqueness Score']]
    matching_pairs_df = matching_pairs_df.loc[matching_pairs_df["Uniqueness Score"] > 0]
    matching_pairs_df = pd.merge(matching_pairs_df, lr_db_df[['ligand', 'receptor', 'annotation']],
                                 left_on=['Ligand', 'Receptor'],
                                 right_on=['ligand', 'receptor'],
                                 how='left')
    matching_pairs_df.drop(['ligand', 'receptor'], axis=1, inplace=True)
    matching_pairs_df = matching_pairs_df.rename(columns={'annotation': 'Pathway Label'})

    return group_name, matching_pairs_df



def extract_sample_name(group_name):
    if 'permutation' in group_name:
        this = group_name.split('_')[0] # If 'permutation' is present, split by underscore and return the first part.
        return this
    else:
        # If 'permutation' is not present, split by space and return the first part.
        that = group_name.split(' ')[0]
        return that

def extract_sample_name_xlsx(group_name):
    thing = group_name.split(' ')[0]
    return thing


def preprocess_counts_file(min_clus_percent, general_filepath, import_folder, counts_file):
    counts = pd.read_csv(f'{general_filepath}{import_folder}/{counts_file}.csv').fillna(0)
    if counts.columns[0] == 'sample':
        counts.set_index('sample', inplace=True)
    total_counts = counts.sum(axis=1)
    counts_percent = counts.div(total_counts, axis=0)
    valid_clusters = counts_percent.columns[(counts_percent > min_clus_percent).any(axis=0)]
    cluster_percent = counts_percent.loc[:, valid_clusters]
    counts_filtered = counts.loc[:, valid_clusters]
    new_counts_list = valid_clusters.tolist()
    return new_counts_list, cluster_percent, counts_filtered

def create_lr_dict(lr_db_df):
    lr_dict = {}
    for _, row in lr_db_df.iterrows():
        ligand = row['ligand']
        receptor = row['receptor']
        subunits = [row[f'subunit_{i}'] for i in range(1, 5) if pd.notnull(row[f'subunit_{i}'])]
        if ligand not in lr_dict:
            lr_dict[ligand] = {}
        lr_dict[ligand][receptor] = subunits
    return lr_dict

def find_matching_lr_pairs(deg_data, lr_dict, cluster_percent, sample_name, counts_filtered, general_filepath):
    rows_list = []
    for ligand, receptors in lr_dict.items():
        if ligand in deg_data['gene'].values:  # Only process ligands present in the DEG data
            ligand_clusters = deg_data.loc[deg_data['gene'] == ligand, 'cluster'].values

            # Iterate over receptors and their associated subunits for a given ligand
            for receptor, subunits in receptors.items():
                receptor_subunits_data = deg_data[deg_data['gene'].isin(subunits)]
                receptor_clusters = receptor_subunits_data['cluster'].unique()

                # Iterate over unique receptor clusters
                for receptor_cluster in receptor_clusters:
                    # Check if all subunits of the receptor come from the same cluster
                    if set(subunits).issubset(receptor_subunits_data[receptor_subunits_data['cluster'] == receptor_cluster]['gene'].values):
                        # Iterate over ligand clusters
                        for ligand_cluster in ligand_clusters:
                            # Create a new row for the matching pair
                            new_row = {
                                'Ligand': ligand,
                                'Receptor': receptor,
                                'Subunits': subunits,
                                'Source': ligand_cluster,
                                'Target': receptor_cluster
                            }
                            # Append the new row data to the list
                            rows_list.append(new_row)

    matching_pairs_df = pd.DataFrame(rows_list)
    if matching_pairs_df.empty:
        # Ensure the DataFrame has the necessary columns, even if it's empty
        matching_pairs_df = pd.DataFrame(columns=['Ligand', 'Receptor', 'Subunits', 'Source', 'Target', 'Communication Score', 'Uniqueness Score'])
        return matching_pairs_df

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(calculate_scores_parallel, [(row, deg_data, cluster_percent, sample_name, counts_filtered, general_filepath) for _, row in matching_pairs_df.iterrows()]))

    matching_pairs_df['Communication Score'] = [res[0] for res in results]
    matching_pairs_df['Uniqueness Score'] = [res[1] for res in results]

    #print(f"Final matching pairs dataframe shape: {matching_pairs_df.shape}")
    #print(f"matching_pairs_df columns after calculation: {matching_pairs_df.columns}")

    return matching_pairs_df




def group_results(results):
    xlsx_files = {}
    h5_groups = {}
    for group_name, dfs in results.items():
        final_df = pd.concat(dfs, ignore_index=True)
        if group_name.endswith('.xlsx'):
            xlsx_files[group_name] = final_df
        else:
            base_name = group_name.split('_permutation')[0]
            if base_name not in h5_groups:
                h5_groups[base_name] = {}
            h5_groups[base_name][group_name] = final_df
    return xlsx_files, h5_groups