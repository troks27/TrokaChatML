"""
Drug Enrichment Module

Functions for enriching drug data with information from ChEMBL.
"""

import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# Try to import ChEMBL downloader
try:
    import chembl_downloader
    CHEMBL_AVAILABLE = True
except ImportError:
    CHEMBL_AVAILABLE = False
    print("Warning: chembl_downloader package not available. Drug enrichment functionality will be limited.")


def test_chembl_connection():
    """
    Test connection to the ChEMBL database.
    
    Returns:
    --------
    bool
        True if connection is successful, False otherwise
    """
    if not CHEMBL_AVAILABLE:
        print("chembl_downloader package is not installed.")
        return False
    
    try:
        with chembl_downloader.cursor() as cursor:
            cursor.execute("SELECT 1")
            result = cursor.fetchone()
            if result[0] == 1:
                print("Successfully connected to the ChEMBL database.")
                return True
            else:
                print("Connected to the database, but received unexpected result.")
                return False
    except Exception as e:
        print(f"Failed to connect to the ChEMBL database: {str(e)}")
        return False


def query_drug_info(drug_id, by_chembl_id=True):
    """
    Query drug information from ChEMBL.
    
    Parameters:
    -----------
    drug_id : str
        ChEMBL ID or drug name to query
    by_chembl_id : bool, default=True
        Whether drug_id is a ChEMBL ID (True) or a drug name (False)
        
    Returns:
    --------
    tuple
        (drug_id, drug_info) where drug_info is a dictionary with drug information
    """
    if not CHEMBL_AVAILABLE:
        return drug_id, {'error': 'chembl_downloader package not available'}
    
    try:
        with chembl_downloader.cursor() as cursor:
            if by_chembl_id:
                query = """
                SELECT DISTINCT
                    md.pref_name AS drug_name,
                    dm.mechanism_of_action,
                    di.mesh_heading AS indication
                FROM
                    molecule_dictionary md
                LEFT JOIN drug_mechanism dm ON md.molregno = dm.molregno
                LEFT JOIN drug_indication di ON md.molregno = di.molregno
                WHERE
                    md.chembl_id = ?
                """
                param = drug_id
            else:
                query = """
                SELECT DISTINCT
                    md.pref_name AS drug_name,
                    dm.mechanism_of_action,
                    di.mesh_heading AS indication
                FROM
                    molecule_dictionary md
                LEFT JOIN drug_mechanism dm ON md.molregno = dm.molregno
                LEFT JOIN drug_indication di ON md.molregno = di.molregno
                WHERE
                    md.pref_name = ?
                """
                param = drug_id
            
            cursor.execute(query, (param,))
            results = cursor.fetchall()
        
        if results:
            mechanisms = list(set(row[1] for row in results if row[1]))
            indications = list(set(row[2] for row in results if row[2]))
            return drug_id, {
                'Mechanism': mechanisms[0] if mechanisms else '',
                'Clinical Trials': '; '.join(indications)
            }
        else:
            return drug_id, None
    except Exception as e:
        return drug_id, {'error': str(e)}


def query_multiple_drugs(drugs, by_chembl_id=True, max_workers=10):
    """
    Query information for multiple drugs in parallel.
    
    Parameters:
    -----------
    drugs : list
        List of drug IDs or names to query
    by_chembl_id : bool, default=True
        Whether drugs contains ChEMBL IDs (True) or drug names (False)
    max_workers : int, default=10
        Maximum number of worker threads
        
    Returns:
    --------
    dict
        Dictionary mapping drug IDs/names to drug information
    """
    if not CHEMBL_AVAILABLE:
        return {drug: {'error': 'chembl_downloader package not available'} for drug in drugs}
    
    drug_info = {}
    max_workers = min(max_workers, len(drugs))
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_drug = {executor.submit(query_drug_info, drug, by_chembl_id): drug for drug in drugs}
        for future in tqdm(future_to_drug, total=len(drugs), desc="Querying drug information"):
            drug, result = future.result()
            if result:
                drug_info[drug] = result
    
    return drug_info


def enrich_drug_data(file_paths, drug_column='node', by_chembl_id=False, max_workers=10):
    """
    Enrich drug data in CSV/TSV files with information from ChEMBL.
    
    Parameters:
    -----------
    file_paths : list
        List of file paths to enrich
    drug_column : str, default='node'
        Column name containing drug IDs/names
    by_chembl_id : bool, default=False
        Whether drug column contains ChEMBL IDs (True) or drug names (False)
    max_workers : int, default=10
        Maximum number of worker threads
        
    Returns:
    --------
    dict
        Dictionary with results for each file
    """
    if not CHEMBL_AVAILABLE:
        return {'error': 'chembl_downloader package not available'}
    
    results = {}
    
    # Collect all unique drugs from files
    all_drugs = set()
    for file_path in file_paths:
        sep = '\t' if file_path.endswith('.tsv') else ','
        df = pd.read_csv(file_path, sep=sep)
        if drug_column in df.columns:
            all_drugs.update(df[drug_column].dropna().unique())
    
    # Query information for all drugs at once
    all_drugs = list(all_drugs)
    drug_info = query_multiple_drugs(all_drugs, by_chembl_id, max_workers)
    
    # Update each file with drug information
    for file_path in file_paths:
        try:
            # Load file
            sep = '\t' if file_path.endswith('.tsv') else ','
            df = pd.read_csv(file_path, sep=sep)
            
            # Skip if drug column not found
            if drug_column not in df.columns:
                results[file_path] = {'error': f"Column '{drug_column}' not found"}
                continue
            
            # Add mechanism and clinical trial columns
            df['Mechanism'] = df[drug_column].map(lambda x: drug_info.get(x, {}).get('Mechanism', ''))
            df['Clinical Trials'] = df[drug_column].map(lambda x: drug_info.get(x, {}).get('Clinical Trials', ''))
            
            # Save updated file
            output_path = file_path.replace('.tsv', '_updated.tsv').replace('.csv', '_updated.csv')
            df.to_csv(output_path, sep=sep, index=False)
            
            results[file_path] = {
                'drugs_found': sum(df[drug_column].isin(drug_info.keys())),
                'total_drugs': df[drug_column].nunique(),
                'output_path': output_path
            }
        except Exception as e:
            results[file_path] = {'error': str(e)}
    
    return results


def download_chembl_database(target_dir=None):
    """
    Download and extract the ChEMBL SQLite database.
    
    Parameters:
    -----------
    target_dir : str, optional
        Directory to save the database, if None will use default location
        
    Returns:
    --------
    str
        Path to the downloaded database
    """
    if not CHEMBL_AVAILABLE:
        raise ImportError("chembl_downloader package is required for this function.")
    
    try:
        if target_dir:
            return chembl_downloader.download_extract_sqlite(dir_path=target_dir)
        else:
            return chembl_downloader.download_extract_sqlite()
    except Exception as e:
        raise Exception(f"Failed to download ChEMBL database: {str(e)}")




import os
import pandas as pd
from tqdm import tqdm
import traceback
import concurrent.futures
import chembl_downloader

def update_drug_info(base_path, folders, tsv_filename, drug_column):
    """
    Update TSV files with drug mechanism and clinical trial information from ChEMBL.

    Parameters:
    - base_path (str): Root directory containing the folders.
    - folders (list): List of folder names containing TSV files.
    - tsv_filename (str): Name of the TSV file in each folder.
    - drug_column (str): Column name in TSV files containing drug names (ChEMBL IDs).
    """

    # Download/extract the ChEMBL SQLite file if necessary
    chembl_downloader.download_extract_sqlite()

    # Step 1: Find all TSV files
    tsv_files = []
    for folder in folders:
        file_path = os.path.join(base_path, folder, tsv_filename)
        if os.path.exists(file_path):
            tsv_files.append(file_path)

    # Step 2: Extract unique drug names from all TSV files
    all_drugs = set()
    for file in tsv_files:
        df = pd.read_csv(file, sep='\t')
        all_drugs.update(df[drug_column].unique())
    all_drugs = list(all_drugs)

    # Step 3: Define the drug query function using chembl_downloader.cursor()
    def query_drug(drug):
        try:
            with chembl_downloader.cursor() as cursor:
                query = """
                SELECT DISTINCT
                    md.pref_name AS drug_name,
                    dm.mechanism_of_action,
                    di.mesh_heading AS indication
                FROM
                    molecule_dictionary md
                LEFT JOIN drug_mechanism dm ON md.molregno = dm.molregno
                LEFT JOIN drug_indication di ON md.molregno = di.molregno
                WHERE
                    md.chembl_id = ?
                """
                cursor.execute(query, (drug,))
                results = cursor.fetchall()

            if results:
                mechanisms = list(set(row[1] for row in results if row[1]))
                indications = list(set(row[2] for row in results if row[2]))
                return drug, {
                    'Mechanism': '; '.join(mechanisms) if mechanisms else '',
                    'Clinical Trials': '; '.join(indications) if indications else ''
                }
            else:
                print(f"No results found in the database for {drug}")
                return drug, None
        except Exception as e:
            print(f"Error querying {drug}: {str(e)}")
            print(f"Traceback: {traceback.format_exc()}")
            return drug, None

    # Query drugs in parallel
    drug_info = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=11) as executor:
        future_to_drug = {executor.submit(query_drug, drug): drug for drug in all_drugs}
        for future in tqdm(concurrent.futures.as_completed(future_to_drug),
                           total=len(all_drugs), desc="Querying drugs"):
            drug, result = future.result()
            if result:
                drug_info[drug] = result
            else:
                print(f"No results or error for {drug}")

    # Step 4: Update each TSV file with new drug information
    def update_tsv_file(file):
        df = pd.read_csv(file, sep='\t')
        df['Mechanism'] = df[drug_column].map(lambda x: drug_info.get(x, {}).get('Mechanism', ''))
        df['Clinical Trials'] = df[drug_column].map(lambda x: drug_info.get(x, {}).get('Clinical Trials', ''))
        updated_file = file.replace('.tsv', '_updated.tsv')
        df.to_csv(updated_file, sep='\t', index=False)
        return updated_file

    # Update TSV files in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=11) as executor:
        future_to_file = {executor.submit(update_tsv_file, file): file for file in tsv_files}
        for future in tqdm(concurrent.futures.as_completed(future_to_file),
                           total=len(tsv_files), desc="Updating TSV files"):
            updated_file = future.result()
            print(f"Updated {updated_file}")
