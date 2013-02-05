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
           
size_temp_filename = "MTEEcto_data.csv"
metabolic_rate_filename = "Class_metabolicrates_Makrievadata.csv"
ST_data = pd.read_csv(size_temp_filename)
MR_data = pd.read_csv(metabolic_rate_filename)
Class_list = set(MR_data['Class'])
for current_class in Class_list:
    class_data = MR_data.ix[MR_data['Class'] == current_class]
    InvK = calculate_InvK(class_data['TC '])
    Mass_kg = convert_grams_to_kilograms(class_data['Mg '])
    LogMass_kg = np.log(Mass_kg)
    Metabolic_rate = np.log(class_data['qWkg'] * Mass_kg)
    analysis_data = pd.DataFrame(zip(InvK, LogMass_kg, Metabolic_rate), 
                                 columns=['InvK', 'LogMass_kg', 'Metabolic_rate'])
    results = sm.OLS.from_formula("Metabolic_rate ~ LogMass_kg + InvK", analysis_data).fit()
    print current_class
    print results.summary()
    






