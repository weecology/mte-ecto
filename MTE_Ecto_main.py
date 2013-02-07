from __future__ import division

import pandas as pd

import math

import MTE_params

def create_species_replicateID_linker(data):
    """finds all unique study identifiers and returns pandas dataframe that links
    unique replicate ID with the species name and taxonomic class of the organism 
    used in that study"""
    
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

def get_metabolic_rate(size, temp, exponent, Ea):
    """function that calculates metabolic rate used Gillooly et al eqn"""
    
    metabolic_rate = (size ** exponent) * math.exp(Ea/(0.00086*(temp + 273.15)))
    return metabolic_rate

"""Main Code"""
metabolic_params = MTE_params.get_class_MTE_params("Class_metabolicrates_Makrievadata.csv", 
                                                   "gillooly_fish.csv", 
                                                   "Whiteetal_Amphibiandata.csv")

SizeTemp_data = pd.read_csv("MTEEcto_data.csv")
replicate_species_DF = create_species_replicateID_linker(SizeTemp_data)
unique_replicates = set(replicate_species_DF['replicateID'])
unique_classes = set(SizeTemp_data['Class'])

for index_class in unique_classes:
    class_sizedata = SizeTemp_data.ix[SizeTemp_data['Class'] == index_class]
    class_MTEparams = metabolic_params.ix[metabolic_params['Class'] == index_class]
    class_replicate = set(class_sizedata['replicate'])
    
    for index_replicate in class_replicate:
        replicate_sizedata = class_sizedata.ix[class_sizedata['replicate'] == 
                                         index_replicate]
        replicate_temps = set(replicate_sizedata['temp'])
        for current_temp in replicate_temps:
            initial_size = replicate_sizedata['mass'][replicate_sizedata['temp']
                                                      == current_temp]
            repeater = len(replicate_temps)
            no_sizechange = [get_metabolic_rate(initial_size, temp, 
                                                float(class_MTEparams['Exponent']), 
                                                float(class_MTEparams['Ea'])) 
                             for temp in replicate_temps]
                              
            
        

        
        


    
    


