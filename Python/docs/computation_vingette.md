# TrokaChat Package: Quick Start Guide

This notebook demonstrates how to use the TrokaChat package to analyze cell-cell communication in single-cell RNA sequencing data.

## Installation

```
pip install trokachat_computation
```

## Step 1: Import and Preprocess Data

The first step uses `TrokaChat_import` to process input files:

```python
import trokachat_computation as tc
from trokachat_computation import TrokaChat_import

# Initialize the importer with file paths
importer = TrokaChat_import(
    input_path='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/csv1',
    file_prefixes=['NG_vs_NG.csv', 'DIAB_vs_NG.csv', 
                   'nulldist_NG_vs_NG.h5', 'nulldist_DIAB_vs_NG.h5'],
    output_path='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/DEG Output',
    pct_min=0.0  # Minimum percentage filter
)

# Process and export the data
importer.import_data()
```

**What happens during import:**
- Files are loaded from the input directory
- Data is filtered based on the minimum percentage threshold
- Processed data is exported to the output directory

## Step 2: Compute Communication Scores

After preprocessing, use `TrokaChat_Computation` to calculate communication scores:

```python
from trokachat_computation import TrokaChat_Computation

# Initialize the computation object
calculator = TrokaChat_Computation(
    import_folder='DEG Output',
    general_filepath='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/',
    output_path='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/HERE/',
    files_list=['NG DEGs.csv', 'DIAB DEGs.csv', 
                'nulldistavgs_NG.h5', 'nulldistavgs_DIAB.h5'],
    species='mouse',  # Determines reference database
    sample_name_from_R='sample',
    counts_file='counts',
    min_clus_percent=0  # Minimum cluster representation
)

# Run the computation
calculator.process_files()
```

**What happens during computation:**
- Reference database is loaded based on the species
- Ligand-receptor pairs are identified across cell clusters
- Communication scores are calculated for each pair
- Results are exported as Excel and H5 files

## Output Files

The analysis produces two types of output:
- Excel files (from CSV inputs): `{condition}_allpathways.xlsx`
- H5 files (from null distribution data): `NULLDIST_{condition}_allpathways.h5`

Each file contains:
- Ligand-receptor pairs
- Source and target cell clusters
- Communication and uniqueness scores
- Pathway annotations

## Complete Example

Here's the full workflow in a single code block:

```python
import trokachat_computation as tc
from trokachat_computation import TrokaChat_import, TrokaChat_Computation

# Step 1: Import and preprocess data
importer = TrokaChat_import(
    input_path='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/csv1',
    file_prefixes=['NG_vs_NG.csv', 'DIAB_vs_NG.csv', 
                  'nulldist_NG_vs_NG.h5', 'nulldist_DIAB_vs_NG.h5'],
    output_path='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/DEG Output',
    pct_min=0.0
)
importer.import_data()

# Step 2: Calculate communication scores
calculator = TrokaChat_Computation(
    import_folder='DEG Output',
    general_filepath='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/',
    output_path='/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/HERE/',
    files_list=['NG DEGs.csv', 'DIAB DEGs.csv', 
               'nulldistavgs_NG.h5', 'nulldistavgs_DIAB.h5'],
    species='mouse',
    sample_name_from_R='sample',
    counts_file='counts',
    min_clus_percent=0
)
calculator.process_files()
```