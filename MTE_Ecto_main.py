"""Main code for calculating MTE slope and Ea and Q4 values for the 
MTE Ecto Project using raw data on size-temperature relationships and 
metabolic rate, body size/temperature"""

from __future__ import division

import csv
import numpy as np
import pandas as pd
import statsmodels.api as sm

def calculate_InvK(temperature):
    """takes temperature in Celsius and converts to 1/kelvin"""
    Inverse_Kelvin = 1/(273.15 + temperature)
    return Inverse_Kelvin

def convert_grams_to_kilograms(grams):
    """converts weights in grams to kilograms"""
    Mass_kg = grams/1000
    return Mass_kg

def extract_analysis_data(temperature_celsius, mass_grams):
    """converts to celsius and gram data to 1/Kelvin and kilograms and returns 
    formatted data"""
    InvK = calculate_InvK(temperature_celsius)
    Mass_kg = convert_grams_to_kilograms(mass_grams)
    LogMass_kg = np.log(Mass_kg)
    return InvK, LogMass_kg, Mass_kg

def get_metabolic_params(InvK, LogMass, Metabolic_rate):
    """conducts multiple regression to get fitted mass exponent (exponent) 
    and activation energy (Ea)"""
    analysis_data = pd.DataFrame(zip(InvK, LogMass_kg, Metabolic_rate),
                                 columns=['InvK', 'LogMass_kg', 'Metabolic_rate'])
    results = sm.OLS.from_formula("Metabolic_rate ~ LogMass_kg + InvK", 
                                  analysis_data).fit()
    Intercept, exponent, temp_slope = results.params
    Ea = abs(temp_slope) * 0.000086
    return exponent, Ea
 
def make_class_parameters_list(param_list, current_class, exponent, Ea):
    """adds analysis parameters to list along with the Class designation"""
    class_params = [current_class, exponent, Ea]
    param_list.append(class_params)
    return param_list

"""Main Body"""
metabolic_rate_filename = "Class_metabolicrates_Makrievadata.csv"

MR_data = pd.read_csv(metabolic_rate_filename)
Class_list = set(MR_data['Class'])
metabolic_params = []

for current_class in Class_list:
    class_data = MR_data.ix[MR_data['Class'] == current_class]
    InvK, LogMass_kg, Mass_kg = extract_analysis_data(class_data['TC '], 
                                                             class_data['Mg '])
    Metabolic_rate = np.log(class_data['qWkg'] * Mass_kg)
    exponent, Ea = get_metabolic_params(InvK, LogMass_kg, Metabolic_rate)
    metabolic_params = make_class_parameters_list(metabolic_params, 
                                                  current_class, exponent, Ea)

metabolic_rate_amphibians = "Whiteetal_Amphibiandata.csv"

amphi_data = pd.read_csv(metabolic_rate_amphibians)
InvK, LogMass_kg, Mass_kg = extract_analysis_data(amphi_data['TC'], amphi_data['Mg'])
Metabolic_rate = np.log(amphi_data['Watts'])
exponent, Ea = get_metabolic_params(InvK, LogMass_kg, Metabolic_rate)
metabolic_params = make_class_parameters_list(metabolic_params, 'Amphibians', 
                                              exponent, Ea)

metabolic_rate_fish = "gillooly_fish.csv"
fish_data = pd.read_csv(metabolic_rate_fish)
LogMass_kg = np.log(convert_grams_to_kilograms(fish_data['Mg']))
exponent, Ea = get_metabolic_params(fish_data['invK'], LogMass_kg, 
                                    np.log(fish_data['W']))
metabolic_params = make_class_parameters_list(metabolic_params, 'Actinoperygii', 
                                              exponent, Ea)

metabolic_params = make_class_parameters_list(metabolic_params, 'Insecta', 0.75,
                                              0.62)

other_classes = ['Gastropoda','Eurotatoria', 'Entognatha']
for current_class in other_classes:
    metabolic_params = make_class_parameters_list(metabolic_params, 
                                                  current_class, 0.75, 0.63)    
    

                    








    
    
    
    
    






