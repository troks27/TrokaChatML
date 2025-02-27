import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing
import time
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice
from .scoring_pipeline_utils import load_h5_file_groups, process_group, group_results, save_h5_file, extract_sample_name, process_single_file
from .scoring_functions import calculate_scores_parallel
from .config import Config


class TrokaChat_Computation:
    def __init__(self, import_folder, general_filepath, output_path, files_list, species, sample_name_from_R, counts_file, min_clus_percent):
        self.import_folder = import_folder
        self.general_filepath = general_filepath
        self.output_path = output_path
        self.files_list = files_list
        self.species = species
        self.sample_name_from_R = sample_name_from_R
        self.counts_file = counts_file
        self.min_clus_percent = min_clus_percent

    def process_files(self):
        start_time = time.time()
        config = Config(self.species)
        lr_db_df = config.load_reference()
        num_cores = multiprocessing.cpu_count()

        # Ensure output directories and files
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            futures = []
            for file_name in self.files_list:
                if file_name.endswith('.h5'):
                    filepath = f"{self.general_filepath}{self.import_folder}/{file_name}"
                    group_generator = load_h5_file_groups(filepath)

                    while True:
                        batch = list(islice(group_generator, num_cores))
                        if not batch:
                            break
                        for group_name, data in batch:
                            futures.append(executor.submit(
                                process_group, group_name, data, lr_db_df, self.general_filepath, self.output_path, 
                                self.min_clus_percent, self.sample_name_from_R, self.counts_file, self.import_folder
                            ))
                elif file_name.endswith('.csv'):
                    group_name = file_name
                    df = pd.read_csv(f"{self.general_filepath}{self.import_folder}/{file_name}")
                    sample_name = extract_sample_name(group_name)
                    futures.append(executor.submit(
                        process_single_file, group_name, df, sample_name, lr_db_df, self.general_filepath, 
                        self.output_path, self.min_clus_percent, self.sample_name_from_R, 
                        self.counts_file, self.import_folder
                    ))

            results = {}
            with tqdm(total=len(futures), desc="Processing groups") as pbar:
                for future in as_completed(futures):
                    group_name, processed_df = future.result()
                    if group_name not in results:
                        results[group_name] = []
                    results[group_name].append(processed_df)
                    pbar.update(1)

            # Save the results correctly
            xlsx_files, h5_groups = group_results(results)

            # Export Excel files (from CSV inputs):
            for file_name, df_list in xlsx_files.items():
                if not isinstance(df_list, list):
                    df_list = [df_list]
                df = pd.concat(df_list, ignore_index=True)
                # Remove any spurious parts if needed
                base_name = file_name.split('.')[0]
                output_filepath = os.path.join(self.output_path, f"{base_name}_allpathways.xlsx")
                df.to_excel(output_filepath, index=False)
                print(f"Saved Excel file: {output_filepath}")

            # Export HDF5 groups (for null distribution data):
            for base_name, group_dfs in h5_groups.items():
                # Only add "NULLDIST_" if the base name doesn’t already include it
                if not base_name.startswith("NULLDIST_"):
                    base_name = f"NULLDIST_{base_name}"
                output_filepath = os.path.join(self.output_path, f"{base_name}_allpathways.h5")
                save_h5_file(output_filepath, group_dfs)
                print(f"Saved HDF5 file: {output_filepath}")

        end_time = time.time()  # Record the end time
        total_time = end_time - start_time  # Calculate the total duration
        print(f"Total time taken: {total_time} seconds")
