#############
## LIBRARIES
#############

import xarray as xr
import netCDF4
import numpy as np
import datetime
import pandas as pd
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import seaborn as sns
import cartopy.feature as cfeature
import os

##################################
## Paths, lists, maps
##################################

data_dir = 'D:/MERCURY_ANU/DATA_ANALYSIS/v14_6_0/Preind_outputs/New_diagnostic_island/'
years = ['2018'] ### The idea is to run 5 years from 2018 to 2022, but the results show that it is practically the same one or five years
perturbation_types = ['NO_PERTUBS', 'OCEAN_LESS', 'OCEAN_MORE', 'SSA_CONSTANT', 'SSA_ZERO', 'WIND_LESS', 'WIND_MORE', 'OX_ALL'] # OX_ALL means Br, BrO and OH have been increased 20% on 50S to 60S band  
data_types_all = ['DryDep', 'WetLossLS', 'WetLossConv', 'MercuryChem', 'MercuryEmis', 'MercuryOcean', 'SpeciesConc', 'StateMet', 'Budget']

perturbation_type_map = {
    'NO_PERTUBS': 'NP',
    'OCEAN_LESS': 'OL',
    'OCEAN_MORE': 'OM',
    'SSA_CONSTANT': 'SSAC',
    'SSA_ZERO': 'SSA0',
    'WIND_LESS': 'WL',
    'WIND_MORE': 'WM',
    'OX_ALL': 'OX'
}

data_type_map = {
    'DryDep': 'DD',  # <<<<
    'WetLossLS': 'WLLS',  # <<<<
    'WetLossConv': 'WLConv', # <<<<
    'MercuryChem': 'HgChem',  # <<<<
    'MercuryEmis': 'HgEmis',
    'MercuryOcean': 'HgOcean',
    'SpeciesConc': 'SpCon',
    'StateMet': 'Met',  # <<<<
    'Budget': 'Budget'
}

##########################################
### FUNCTIONS
##########################################

#### FUNCTION: get_dataset

def get_dataset(year, type_name, data_type):
    """
    Loading a dataset xarray specific according the provided year, pertubation_type, and data_type
    
    Args:
        year (str): data year (es. '2018').
        type_name (str): Simulation type (es. 'NO_PERTUBS').
        data_type (str): Data type  (es. 'DryDep').

    Returns:
        xarray.Dataset: The loaded dataset, o None 
    """
    # Building the path
    path_key = f"{perturbation_type_map[type_name]}_{year}_{data_type_map[data_type]}"
    full_path = os.path.join(data_dir, year, type_name, f'GEOSChem.{data_type}.{year}*.nc4')
    
    print(f"Files for {path_key}...")
    file_list = sorted(glob.glob(full_path))
    
    if not file_list:
        print(f"No files found {path_key}.")
        return None
    
    print(f"Loading {len(file_list)} file for {path_key}...")
    try:
        ds = xr.open_mfdataset(file_list, combine='by_coords', parallel=True)
        print("Loading Completed")
        return ds
    except Exception as e:
        print(f"Loading error for the dataset in {path_key}: {e}")
        return None



##########################################
### CONSTANTS AND CONVERSION FACTORS 
##########################################

# TIME

s_in_month = 2.628e6
s_in_yr = 3.154e7 

#MASS
g_kg = 1e3
ug_g = 1e6
ug_kg = 1e9

#AREA amd Volume
cm2_m2 = 1e4
cm3_m3 = 1e6

#Chemistry 
MW_Hg = 200.59      ## <<< Mercury Molecular Weight
NA = 6.023e23       ## <<< Avogadro's number

# CONVERSION FACTOR which don't need inquiring the sim outputs
cf_units_dd_month = (MW_Hg/NA) * ug_g * cm2_m2 * s_in_month    ### <<<  from molecule/cm2*s   to micrograms/m2*month
cf_units_dd_year  = (MW_Hg/NA) * ug_g * cm2_m2 * s_in_yr       ### <<<  from molecule/cm2*s   to micrograms/m2*year

conv = (MW_Hg/NA) * ug_g * 1e6 * s_in_yr
print(conv)

##########################################
### Plotting purposes
##########################################

# Defining maps' colors and line style 

colors = {
    'NP': 'black', 'WM': 'blue', 'WL': 'lightblue',
    'OM': 'green', 'OL': 'lightgreen', 'SSA0': 'red', 'SSAC': 'yellow', 'OX': 'orange'
}

linestyles = {
    'NP': '--', 'WM': '-', 'WL': '-',
    'OM': '-', 'OL': '-', 'SSA0': '-', 'SSAC': '-', 'OX': '-'
}

legend_map = {
    'NP': 'No Perturbations',
    'WM': 'More_Wind',
    'WL': 'Less_Wind',
    'OM': 'More_Ocean',
    'OL': 'Less_Ocean',
    'SSA0': 'Sea Salt Aerosols (Zero)',
    'SSAC': 'Sea Salt Aerosols (Constant)',
    'OX': 'Oxidants Br,BrO,OH +20%'
}

month_names = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']



##########################################
### Chemical species lists 
##########################################

hg2_species_list_sp = [
    'SpeciesConcVV_HgCl2',
    'SpeciesConcVV_HgOHOH',
    'SpeciesConcVV_HgOHBrO',
    'SpeciesConcVV_HgOHClO',
    'SpeciesConcVV_HgOHHO2',
    'SpeciesConcVV_HgOHNO2',
    'SpeciesConcVV_HgClOH',
    'SpeciesConcVV_HgClBr',
    'SpeciesConcVV_HgClBrO',
    'SpeciesConcVV_HgClClO',
    'SpeciesConcVV_HgClHO2',
    'SpeciesConcVV_HgClNO2',
    'SpeciesConcVV_HgBr2',
    'SpeciesConcVV_HgBrOH',
    'SpeciesConcVV_HgBrClO',
    'SpeciesConcVV_HgBrBrO',
    'SpeciesConcVV_HgBrHO2',
    'SpeciesConcVV_HgBrNO2'
]

hg2_species_list_DD = [    ##<<< They are 18 species
    'DryDep_HgCl2',
    'DryDep_HgOHOH',
    'DryDep_HgOHBrO',
    'DryDep_HgOHClO',
    'DryDep_HgOHHO2',
    'DryDep_HgOHNO2',
    'DryDep_HgClOH',
    'DryDep_HgClBr',
    'DryDep_HgClBrO',
    'DryDep_HgClClO',
    'DryDep_HgClHO2',
    'DryDep_HgClNO2',
    'DryDep_HgBr2',
    'DryDep_HgBrOH',
    'DryDep_HgBrClO',
    'DryDep_HgBrBrO',
    'DryDep_HgBrHO2',
    'DryDep_HgBrNO2'
]

hg2_species_list_DD_island = [    ##<<< They are 18 species
    'DryDepIsland_HgCl2',
    'DryDepIsland_HgOHOH',
    'DryDepIsland_HgOHBrO',
    'DryDepIsland_HgOHClO',
    'DryDepIsland_HgOHHO2',
    'DryDepIsland_HgOHNO2',
    'DryDepIsland_HgClOH',
    'DryDepIsland_HgClBr',
    'DryDepIsland_HgClBrO',
    'DryDepIsland_HgClClO',
    'DryDepIsland_HgClHO2',
    'DryDepIsland_HgClNO2',
    'DryDepIsland_HgBr2',
    'DryDepIsland_HgBrOH',
    'DryDepIsland_HgBrClO',
    'DryDepIsland_HgBrBrO',
    'DryDepIsland_HgBrHO2',
    'DryDepIsland_HgBrNO2'
]

## There are 18 from Hg(2+), 2 from Hg(P) ad 1 from Hg(0) => 21 variables for the total one 

##########################################
### 
##########################################






##########################################
### 
##########################################





##########################################
### 
##########################################





##########################################
### 
##########################################





##########################################
### 
##########################################





##########################################
### 
##########################################





##########################################
### 
##########################################





##########################################
### 
##########################################





##########################################
### 
##########################################





##########################################
### 
##########################################
