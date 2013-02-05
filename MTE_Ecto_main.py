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
    InvK = calculate_InvK(temperature_celsius)
    Mass_kg = convert_grams_to_kilograms(mass_grams)
    LogMass_kg = np.log(Mass_kg)
    return InvK, LogMass_kg, Mass_kg

def get_metabolic_params(InvK, LogMass, Metabolic_rate):
    analysis_data = pd.DataFrame(zip(InvK, LogMass_kg, Metabolic_rate),
                                 columns=['InvK', 'LogMass_kg', 'Metabolic_rate'])
    results = sm.OLS.from_formula("Metabolic_rate ~ LogMass_kg + InvK", analysis_data).fit()
    Intercept, exponent, temp_slope = results.params
    Ea = abs(temp_slope) * 0.000086
    return exponent, Ea, temp_slope
 
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
    class_params = [current_class, exponent, Ea]
    metabolic_params.append(class_params)

metabolic_rate_amphibians = "Whiteetal_Amphibiandata.csv"

amphi_data = pd.read_csv(metabolic_rate_amphibians)
InvK, LogMass_kg, Mass_kg = extract_analysis_data(amphi_data['TC'], amphi_data['Mg'])
Metabolic_rate = np.log(amphi_data['Watts'])
exponent, Ea, temp_slope = get_metabolic_params(InvK, LogMass_kg, Metabolic_rate)




    
    
    
    
    






