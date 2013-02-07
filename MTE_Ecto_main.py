from __future__ import division

import pandas as pd

import MTE_params

def create_species_replicateID_linker(data):
    replicates = list(set(data['replicate']))
    replicate_species_linker = []
    
    for replicate in replicates:
        replicate_data = data.ix[data['replicate'] == replicate]
        replicate_species = list(set(replicate_data['Species']))
        replicate_class = list(set(replicate_data['Class']))
        replicate_info = [replicate, replicate_species, replicate_class]
        replicate_species_linker.append(replicate_info)
        replicate_species_link = pd.DataFrame(replicate_species_linker, 
                                              columns=["replicateID", "species",
                                                       "class"])
    return replicate_species_link

"""Main Code"""
metabolic_params = MTE_params.get_class_MTE_params("Class_metabolicrates_Makrievadata.csv", 
                                                   "gillooly_fish.csv", 
                                                   "Whiteetal_Amphibiandata.csv")

SizeTemp_data = pd.read_csv("MTEEcto_data.csv")
replicate_species_DF = create_species_replicateID_linker(SizeTemp_data)
unique_replicates = set(replicate_species_DF['replicateID'])


    
    


