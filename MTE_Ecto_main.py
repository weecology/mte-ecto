from __future__ import division

import pandas as pd

import math

import MTE_params

def create_species_studyID_linker(data):
    """finds all unique study identifiers and returns pandas dataframe that links
    unique study ID with the species name and taxonomic class of the organism 
    used in that study"""
    
    studies = list(set(data['studyID']))
    study_species_linker = []
    
    for study in studies:
        study_data = data.ix[data['studyID'] == study]
        study_species = list(set(study_data['Species']))
        study_class = list(set(study_data['Class']))
        study_info = [study, study_species, study_class]
        study_species_linker.append(study_info)
        study_species_link = pd.DataFrame(study_species_linker, 
                                          columns=["studyID", "species", "class"])
    return study_species_link

def get_metabolic_rate(size, temp, class_parameters):
    """function that calculates metabolic rate used Gillooly et al eqn"""
    extract_exponent = list(class_parameters['Exponent'])
    exponent = float(class_parameters['Exponent'])
    Ea = float(class_parameters['Ea'])
    metabolic_rate = float(math.pow(size,exponent) * math.exp((-Ea/(0.000086*(temp + 273.15)))))
    return metabolic_rate

def get_Class_size_temp_data(index_class, data):
    """subsets main dataframe table and returns just the data for a given 
    taxonomic class"""
    class_sizedata = data.ix[data['Class'] == index_class]
    class_study = set(class_sizedata['studyID'])
    return class_study, class_sizedata
    
def get_SizeTemp_metabolism(data, parameters):
    """Calculates metabolic rate for each paired size and temperature in the data
    and returns a dataframe"""
    output = []
    for row_index, row in data.iterrows():
        metabolic_rate = get_metabolic_rate(row['mass'], row['temp'], parameters)
        row_output = [metabolic_rate, float(row['temp']), float(row['mass'])]
        output.append(row_output)
    output_df = pd.DataFrame(output, columns=["Metabolic", "Temperature", "Size"])
    return output_df

def get_tempchange_metabolism(baseline_size, data, parameters):
    all_temps = set(data['temp'])
    Tchange_MRs = [(get_metabolic_rate(baseline_size, temp, parameters),
                    temp, baseline_size) for temp in all_temps]
    return Tchange_MRs

#    baseline_size = float(data['mass'][data['temp'] == index_temp])
#    Initial_MR = get_metabolic_rate(baseline_size, index_temp, parameters)
#    return Initial_MR
"""Main Code"""
metabolic_params = MTE_params.get_class_MTE_params("Class_metabolicrates_Makrievadata.csv", 
                                                   "gillooly_fish.csv", 
                                                   "Whiteetal_Amphibiandata.csv")

SizeTemp_data = pd.read_csv("MTEEcto_data.csv")
study_species_DF = create_species_studyID_linker(SizeTemp_data)
unique_studys = set(study_species_DF['studyID'])
unique_classes = set(SizeTemp_data['Class'])
Q_Tchange_results = []
for index_class in unique_classes:
    class_studies, class_sizedata = get_Class_size_temp_data(index_class, SizeTemp_data)
    class_MTEparams = metabolic_params.ix[metabolic_params['Class'] == index_class] 
    for index_study in class_studies:
        study_sizedata = class_sizedata.ix[class_sizedata['studyID'] == index_study]
        study_sizes = set(study_sizedata['mass'])
        SizeTemp_MR = get_SizeTemp_metabolism(study_sizedata, class_MTEparams)
        for current_size in study_sizes:
            Tchange_MRs = get_tempchange_metabolism(current_size,study_sizedata, 
                                                    class_MTEparams)
            for i in range(len(Tchange_MRs)):
                for j in range(i+1, len(Tchange_MRs)):
                    baseline_record = Tchange_MRs[i]
                    comparative_record = Tchange_MRs[j]
                    temp_change = comparative_record[1] - baseline_record[1] 
                    if temp_change > 0:
                        Q_Tchange = math.pow(comparative_record[0]/baseline_record[0],
                                             4/temp_change)
                        Q_Tchange_output= [index_study, temp_change, Q_Tchange]
                        Q_Tchange_results.append(Q_Tchange_output)
    
                    
                    
            
            
            
    
            
            
            
            
            
                                                
           
            
        

        
        


    
    


