import pkg_resources
import pandas as pd

class Config:
    def __init__(self, species):
        self.species = species

    def load_reference(self):
        if self.species == 'human':
            resource_path = 'data/human_ALL_PATHWAYS_UPDATED.xlsx'
        elif self.species == 'mouse':
            resource_path = 'data/All_Signaling_Pathways.xlsx'
        else:
            raise ValueError("Invalid species. Expected 'human' or 'mouse'.")

        # Use pkg_resources to get the full path to the resource file
        reference_path = pkg_resources.resource_filename(__name__, resource_path)
        lr_db_df = pd.read_excel(reference_path, sheet_name=0)
        return lr_db_df