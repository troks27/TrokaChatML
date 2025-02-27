from trokachat_computation.scoring_import import TrokaChat_import
from trokachat_computation.scoring_compute import TrokaChat_Computation

def main_import():
    filename = TrokaChat_import(
        input_path='/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat 4.15.24/TrokaChat/csv1_TEST',
        file_prefixes=['NL_vs_NL.csv', 'LS_vs_NL.csv', 'nulldistavgs_NL_vs_NL.h5', 'nulldistavgs_LS_vs_NL.h5'],
        output_path='/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat 4.15.24/TrokaChat/DEG Output',
        pct_min=0.0
    )
    filename.import_data()

def main_computation():
    TrokaObject = TrokaChat_Computation(
        import_folder='DEG Output',
        general_filepath='/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat 4.15.24/TrokaChat/',
        output_path='/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat 4.15.24/TrokaChat/HERE/',
        files_list=['NL DEGs.xlsx', 'LS DEGs.xlsx', 'nulldistavgs_NL.h5', 'nulldistavgs_LS.h5'],
        species='human',
        sample_name_from_R='sample',
        counts_file='counts',
        min_clus_percent=0
    )
    TrokaObject.process_files()

if __name__ == "__main__":
    main_import()
    main_computation()