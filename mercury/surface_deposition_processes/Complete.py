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
###  DATA ANALYSIS 
###
##########################################

# Quindi da datasets_by_year['NP_2018_DD'] puoi recuperare DryDep_Hg0 per esempio
# Quindi da datasets_by_year['NP_2018_SpCon'] puoi recuperare SpeciesConcVV_Hg0 per esempio
# Quindi da datasets_by_year['NP_2018_HgChem'] puoi recuperare Hg2GasToSSA per esempio

# A bit slow, but should be good

data_types = ['SpeciesConc', 'DryDep', 'MercuryChem', 'StateMet']
data_types_WD = ['WetLossLS', 'WetLossConv'] # Poiché questo blocco di codice si concentra solo su WetLoss

datasets_by_year = {}
datasets_by_year_WD = {}

for y in years:
    
    for t in perturbation_types:
        
        for d in data_types:
            
            # Build the key dynamically
            path_key = f"{perturbation_type_map[t]}_{y}_{data_type_map[d]}"
            
            # Load the dataset on-demand
            dataset = get_dataset(year=y, type_name=t, data_type=d)
            
            # Save the dataset in the main dictionary only if it's valid
            if dataset is not None:
                datasets_by_year[path_key] = dataset

        for e in data_types_WD:
            
            # Costruisci la chiave in modo dinamico
            path_key_WD = f"{perturbation_type_map[t]}_{y}_{data_type_map[e]}"

            # Carica il dataset on-demand
            dataset_WD = get_dataset(year=y, type_name=t, data_type=e)

            # Salva il dataset nel dizionario principale solo se è valido
            if dataset_WD is not None:
                datasets_by_year_WD[path_key_WD] = dataset_WD

print("\n--- Loading completed ---")
print(f"{len(datasets_by_year)} datasets have been loaded.")
print("Dictionary keys are:")
print(datasets_by_year.keys())

print("\n--- Loading completed ---")
print(f"{len(datasets_by_year_WD)} datasets have been loaded.")
print("Dictionary keys are:")
print(datasets_by_year_WD.keys())

##########################################
### CREATING NEW VARIABLES
##########################################

# Iterate over each key-value pair in your dictionary

for key, ds in datasets_by_year.items():
    
    # Only process datasets that end with '_DD'
    
    if key.endswith('_DD'):
        print(f"Processing dataset for key: {key}")
    
    ##### 1. Calculate DryDep_Total
        # We find all variables starting with 'DryDep_' but not 'DryDepVel_'
        
        all_drydep_vars = [var for var in ds.data_vars if 'DryDep_' in var and 'DryDepVel_' not in var] 
        
        # Sum all selected DryDep DataArrays. xarray handles alignment automatically.
        ds['DryDep_Total'] = ds[all_drydep_vars].to_array().sum('variable', skipna=False)
    
        # Add a descriptive attribute
        ds['DryDep_Total'].attrs['long_name'] = 'Total Dry Deposition of all species'
        
        print(f"  - DryDep_Total created. It is the sum of {len(all_drydep_vars)} variables.")

        all_drydep_vars_I = [var for var in ds.data_vars if 'DryDepIsland_' in var and 'DryDepVelIsland_' not in var] 
        
        # Sum all selected DryDep DataArrays. xarray handles alignment automatically.
        ds['DryDepIsland_Total'] = ds[all_drydep_vars_I].to_array().sum('variable', skipna=False)
    
        # Add a descriptive attribute
        ds['DryDepIsland_Total'].attrs['long_name'] = 'Total Dry Deposition of all species on the island'
        
        print(f"  - DryDepIsland_Total created. It is the sum of {len(all_drydep_vars)} variables.")
    
    ##### 2. Calculate DryDep_Hg2
        # Create a new DataArray by summing the specified species.
        # Note: `ds[hg2_species_list]` returns a new Dataset with only these variables
        # The .to_array().sum() converts the Dataset to a DataArray and sums over the new dimension.
        
        valid_hg2_species = [var for var in hg2_species_list_DD if var in ds.data_vars]
        if valid_hg2_species:
            ds['DryDep_Hg2'] = ds[valid_hg2_species].to_array().sum('variable', skipna=False)
            ds['DryDep_Hg2'].attrs['long_name'] = 'Total Dry Deposition of Hg(2+) species'
            print(f"  - DryDep_Hg2 created. Summed {len(valid_hg2_species)} species.")
        else:
            print("  - Warning: No valid Hg2+ species found in this dataset. Skipping DryDep_Hg2.")

        valid_hg2_species_I = [var for var in hg2_species_list_DD_island if var in ds.data_vars]
        if valid_hg2_species_I:
            ds['DryDepIsland_Hg2'] = ds[valid_hg2_species_I].to_array().sum('variable', skipna=False)
            ds['DryDepIsland_Hg2'].attrs['long_name'] = 'Total Dry Deposition of Hg(2+) species on the island'
            print(f"  - DryDepIsland_Hg2 created. Summed {len(valid_hg2_species_I)} species.")
        else:
            print("  - Warning: No valid Hg2+ species found in this dataset. Skipping DryDepIsland_Hg2.")
    
    ##### 3. Calculate DryDep_HgP
        # Find all variables that end with 'P'
        hgp_species = [var for var in ds.data_vars if var.endswith('P') and 'DryDep_' in var and 'DryDepVel_' not in var]
        
        # Sum the selected variables
        if hgp_species:
            ds['DryDep_HgP'] = ds[hgp_species].to_array().sum('variable', skipna=False)
            ds['DryDep_HgP'].attrs['long_name'] = 'Total Dry Deposition of Particulate species'
            print(f"  - DryDep_HgP created. Summed {len(hgp_species)} particulate species.")
        else:
            print("  - Warning: No particulate species found in this dataset. Skipping DryDep_HgP.")

        # Find all variables that end with 'P'
        hgp_species_I = [var for var in ds.data_vars if var.endswith('P') and 'DryDepIsland_' in var and 'DryDepVelIsland_' not in var]
        
        # Sum the selected variables
        if hgp_species_I:
            ds['DryDepIsland_HgP'] = ds[hgp_species_I].to_array().sum('variable', skipna=False)
            ds['DryDepIsland_HgP'].attrs['long_name'] = 'Total Dry Deposition of Particulate species on the island'
            print(f"  - DryDepIsland_HgP created. Summed {len(hgp_species_I)} particulate species.")
        else:
            print("  - Warning: No particulate species found in this dataset. Skipping DryDep_HgP.")
    
        print("------------------------------------------")

# Find the first key ending in '_DD' to print a valid example
first_dd_key = next((key for key in datasets_by_year.keys() if key.endswith('_DD')), None)

if first_dd_key:
    print(f"\nFinal variables for a processed dataset ({first_dd_key}):")
    print(list(datasets_by_year[first_dd_key].variables))
    print(f"\nExample of new variable `DryDep_Total`:")
    print(datasets_by_year[first_dd_key]['DryDep_Total'])
else:
    print("\nNo '_DD' datasets were found to process or print.")


##########################################
### MEANS CALCS 
##########################################

# loops: Datarrays for each perturbation type 
data_arrays_by_type_DD = {u: [] for u in perturbation_type_map.values()}  # <<< total
data_arrays_by_type_DD_Hg2 = {u: [] for u in perturbation_type_map.values()}
data_arrays_by_type_DD_HgP = {u: [] for u in perturbation_type_map.values()}
data_arrays_by_type_DD_Hg0 = {u: [] for u in perturbation_type_map.values()}

data_arrays_by_type_DD_I = {u: [] for u in perturbation_type_map.values()}  # <<< total
data_arrays_by_type_DD_Hg2_I = {u: [] for u in perturbation_type_map.values()}
data_arrays_by_type_DD_HgP_I = {u: [] for u in perturbation_type_map.values()}
data_arrays_by_type_DD_Hg0_I = {u: [] for u in perturbation_type_map.values()}

# Dictionary to save the final monthly and annual mean values 

monthly_means_by_type_DD = {}
annual_means_by_type_DD = {}

monthly_means_by_type_DD_Hg2 = {}
annual_means_by_type_DD_Hg2 = {}

monthly_means_by_type_DD_HgP = {}
annual_means_by_type_DD_HgP = {}

monthly_means_by_type_DD_Hg0 = {}
annual_means_by_type_DD_Hg0  = {}

###ISLAND
monthly_means_by_type_DD_I = {}
annual_means_by_type_DD_I = {}

monthly_means_by_type_DD_Hg2_I = {}
annual_means_by_type_DD_Hg2_I = {}

monthly_means_by_type_DD_HgP_I = {}
annual_means_by_type_DD_HgP_I = {}

monthly_means_by_type_DD_Hg0_I = {}
annual_means_by_type_DD_Hg0_I  = {}

#### STANDARD DEVIATION

monthly_std_by_type_DD = {}
annual_std_by_type_DD = {}

monthly_std_by_type_DD_Hg2 = {}
annual_std_by_type_DD_Hg2 = {}

monthly_std_by_type_DD_HgP = {}
annual_std_by_type_DD_HgP = {}

monthly_std_by_type_DD_Hg0 = {}
annual_std_by_type_DD_Hg0 = {}

###ISLAND
monthly_std_by_type_DD_I = {}
annual_std_by_type_DD_I = {}

monthly_std_by_type_DD_Hg2_I = {}
annual_std_by_type_DD_Hg2_I = {}

monthly_std_by_type_DD_HgP_I = {}
annual_std_by_type_DD_HgP_I = {}

monthly_std_by_type_DD_Hg0_I = {}
annual_std_by_type_DD_Hg0_I = {}


###########################################################################
# First loop: Datarrays for each perturbation type 

for t, t_key in perturbation_type_map.items():
    
    for y in years:
        
        full_key_DD = f"{t_key}_{y}_DD"

        #### DRYDEP 
        
        if full_key_DD in datasets_by_year:
            data_array_DD = datasets_by_year[full_key_DD]['DryDep_Total']
            data_array_DD_I = datasets_by_year[full_key_DD]['DryDepIsland_Total']
            data_arrays_by_type_DD[t_key].append(data_array_DD)
            data_arrays_by_type_DD_I[t_key].append(data_array_DD_I)
        else:
            print(f"Warning: The key {full_key_DD} not found .")
        
        if full_key_DD in datasets_by_year:
            data_array_DD_Hg2 = datasets_by_year[full_key_DD]['DryDep_Hg2']
            data_array_DD_Hg2_I = datasets_by_year[full_key_DD]['DryDepIsland_Hg2']
            data_arrays_by_type_DD_Hg2[t_key].append(data_array_DD_Hg2)
            data_arrays_by_type_DD_Hg2_I[t_key].append(data_array_DD_Hg2_I)
        else:
            print(f"Warning: The key {full_key_DD} not found .")
        
        if full_key_DD in datasets_by_year:
            data_array_DD_HgP = datasets_by_year[full_key_DD]['DryDep_HgP']
            data_array_DD_HgP_I = datasets_by_year[full_key_DD]['DryDepIsland_HgP']
            data_arrays_by_type_DD_HgP[t_key].append(data_array_DD_HgP)
            data_arrays_by_type_DD_HgP_I[t_key].append(data_array_DD_HgP_I)
        else:
            print(f"Warning: The key {full_key_DD} not found .")
            
        if full_key_DD in datasets_by_year:
            data_array_DD_Hg0 = datasets_by_year[full_key_DD]['DryDep_Hg0']
            data_array_DD_Hg0_I = datasets_by_year[full_key_DD]['DryDepIsland_Hg0']
            data_arrays_by_type_DD_Hg0[t_key].append(data_array_DD_Hg0)
            data_arrays_by_type_DD_Hg0_I[t_key].append(data_array_DD_Hg0_I)
        else:
            print(f"Warning: The key {full_key_DD} not found .")

#################################################################################
# Second loop: Calculating the monthly and annual mean for each perturbation type 

###### DryDep_Total

for t_key, data_arrays_list_DD in data_arrays_by_type_DD.items():
    
    if data_arrays_list_DD:
        
        combined_data_DD = xr.concat(data_arrays_list_DD, dim='time')
        
        # Calcolo della media annuale (DD)
        
        annual_mean_dd = combined_data_DD.mean(dim='time', skipna=True)
        annual_std_dd = combined_data_DD.std(dim='time', skipna=True)
        
        annual_means_by_type_DD[t_key] = annual_mean_dd
        annual_std_by_type_DD[t_key] = annual_std_dd
        
        print(f"DD - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD)
        
        monthly_mean_dd = combined_data_DD.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd = combined_data_DD.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD[t_key] = monthly_mean_dd
        monthly_std_by_type_DD[t_key] = monthly_std_dd
        
        print(f"DD - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD - Nessun dato trovato per {t_key}.")

###### DryDep_Total on th eIsland

for t_key, data_arrays_list_DD_I in data_arrays_by_type_DD_I.items():
    
    if data_arrays_list_DD_I:
        
        combined_data_DD_I = xr.concat(data_arrays_list_DD_I, dim='time')
        
        # Calcolo della media annuale (DD_I)
        
        annual_mean_dd_i = combined_data_DD_I.mean(dim='time', skipna=True)
        annual_std_dd_i = combined_data_DD_I.std(dim='time', skipna=True)
        
        annual_means_by_type_DD_I[t_key] = annual_mean_dd_i
        annual_std_by_type_DD_I[t_key] = annual_std_dd_i
        
        print(f"DD_I - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD_I)
        
        monthly_mean_dd_i = combined_data_DD_I.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd_i = combined_data_DD_I.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD_I[t_key] = monthly_mean_dd_i
        monthly_std_by_type_DD_I[t_key] = monthly_std_dd_i
        
        print(f"DD_I - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD_I - Nessun dato trovato per {t_key}.")

###### DryDep_Hg2

for t_key, data_arrays_list_DD_Hg2 in data_arrays_by_type_DD_Hg2.items():
    
    if data_arrays_list_DD_Hg2:
        
        combined_data_DD_Hg2 = xr.concat(data_arrays_list_DD_Hg2, dim='time')
        
        # Calcolo della media annuale (DD_Hg2)
        
        annual_mean_dd_Hg2 = combined_data_DD_Hg2.mean(dim='time', skipna=True)
        annual_std_dd_Hg2 = combined_data_DD_Hg2.std(dim='time', skipna=True)
        
        annual_means_by_type_DD_Hg2[t_key] = annual_mean_dd_Hg2
        annual_std_by_type_DD_Hg2[t_key] = annual_std_dd_Hg2
        
        print(f"DD_Hg2 - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD_HG2)
        
        monthly_mean_dd_Hg2 = combined_data_DD_Hg2.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd_Hg2 = combined_data_DD_Hg2.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD_Hg2[t_key] = monthly_mean_dd_Hg2
        monthly_std_by_type_DD_Hg2[t_key] = monthly_std_dd_Hg2
        
        print(f"DD_Hg2 - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD_Hg2 - Nessun dato trovato per {t_key}.")

###### DryDep_Hg2 on the island

for t_key, data_arrays_list_DD_Hg2_I in data_arrays_by_type_DD_Hg2_I.items():
    
    if data_arrays_list_DD_Hg2_I:
        
        combined_data_DD_Hg2_I = xr.concat(data_arrays_list_DD_Hg2_I, dim='time')
        
        # Calcolo della media annuale (DD_Hg2_I)
        
        annual_mean_dd_Hg2_i = combined_data_DD_Hg2_I.mean(dim='time', skipna=True)
        annual_std_dd_Hg2_i = combined_data_DD_Hg2_I.std(dim='time', skipna=True)
        
        annual_means_by_type_DD_Hg2_I[t_key] = annual_mean_dd_Hg2_i
        annual_std_by_type_DD_Hg2_I[t_key] = annual_std_dd_Hg2_i
        
        print(f"DD_Hg2_I - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD_Hg2_I)
        
        monthly_mean_dd_Hg2_i = combined_data_DD_Hg2_I.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd_Hg2_i = combined_data_DD_Hg2_I.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD_Hg2_I[t_key] = monthly_mean_dd_Hg2_i
        monthly_std_by_type_DD_Hg2_I[t_key] = monthly_std_dd_Hg2_i
        
        print(f"DD_Hg2_I - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD_Hg2_I - Nessun dato trovato per {t_key}.")

###### DryDep_HgP

for t_key, data_arrays_list_DD_HgP in data_arrays_by_type_DD_HgP.items():
    
    if data_arrays_list_DD_HgP:
        
        combined_data_DD_HgP = xr.concat(data_arrays_list_DD_HgP, dim='time')
        
        # Calcolo della media annuale (DD_HgP)
        
        annual_mean_dd_HgP = combined_data_DD_HgP.mean(dim='time', skipna=True)
        annual_std_dd_HgP = combined_data_DD_HgP.std(dim='time', skipna=True)
        
        annual_means_by_type_DD_HgP[t_key] = annual_mean_dd_HgP
        annual_std_by_type_DD_HgP[t_key] = annual_std_dd_HgP
        
        print(f"DD_HgP - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD_HgP)
        
        monthly_mean_dd_HgP = combined_data_DD_HgP.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd_HgP = combined_data_DD_HgP.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD_HgP[t_key] = monthly_mean_dd_HgP
        monthly_std_by_type_DD_HgP[t_key] = monthly_std_dd_HgP
        
        print(f"DD_Hg2 - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD_Hg2 - Nessun dato trovato per {t_key}.")

###### DryDep_HgP on the island

for t_key, data_arrays_list_DD_HgP_I in data_arrays_by_type_DD_HgP_I.items():
    
    if data_arrays_list_DD_HgP_I:
        
        combined_data_DD_HgP_I = xr.concat(data_arrays_list_DD_HgP_I, dim='time')
        
        # Calcolo della media annuale (DD_HgP_I)
        
        annual_mean_dd_HgP_i = combined_data_DD_HgP_I.mean(dim='time', skipna=True)
        annual_std_dd_HgP_i = combined_data_DD_HgP_I.std(dim='time', skipna=True)
        
        annual_means_by_type_DD_HgP_I[t_key] = annual_mean_dd_HgP_i
        annual_std_by_type_DD_HgP_I[t_key] = annual_std_dd_HgP_i
        
        print(f"DD_HgP_I - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD_HgP_I)
        
        monthly_mean_dd_HgP_i = combined_data_DD_HgP_I.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd_HgP_i = combined_data_DD_HgP_I.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD_HgP_I[t_key] = monthly_mean_dd_HgP_i
        monthly_std_by_type_DD_HgP_I[t_key] = monthly_std_dd_HgP_i
        
        print(f"DD_HgP_I - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD_HgI - Nessun dato trovato per {t_key}.")


###### DryDep_Hg0

for t_key, data_arrays_list_DD_Hg0 in data_arrays_by_type_DD_Hg0.items():
    
    if data_arrays_list_DD_Hg0:
        
        combined_data_DD_Hg0 = xr.concat(data_arrays_list_DD_Hg0, dim='time')
        
        # Calcolo della media annuale (DD_Hg0)
        
        annual_mean_dd_Hg0 = combined_data_DD_Hg0.mean(dim='time', skipna=True)
        annual_std_dd_Hg0 = combined_data_DD_Hg0.std(dim='time', skipna=True)
        
        annual_means_by_type_DD_Hg0[t_key] = annual_mean_dd_Hg0
        annual_std_by_type_DD_Hg0[t_key] = annual_std_dd_Hg0
        
        print(f"DD_Hg0 - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD_Hg0)
        
        monthly_mean_dd_Hg0 = combined_data_DD_Hg0.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd_Hg0 = combined_data_DD_Hg0.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD_Hg0[t_key] = monthly_mean_dd_Hg0
        monthly_std_by_type_DD_Hg0[t_key] = monthly_std_dd_Hg0
        
        print(f"DD_Hg0 - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD_Hg0 - Nessun dato trovato per {t_key}.")

###### DryDep_Hg0 on the island

for t_key, data_arrays_list_DD_Hg0_I in data_arrays_by_type_DD_Hg0_I.items():
    
    if data_arrays_list_DD_Hg0_I:
        
        combined_data_DD_Hg0_I = xr.concat(data_arrays_list_DD_Hg0_I, dim='time')
        
        # Calcolo della media annuale (DD_Hg0_I)
        
        annual_mean_dd_Hg0_i = combined_data_DD_Hg0_I.mean(dim='time', skipna=True)
        annual_std_dd_Hg0_i = combined_data_DD_Hg0_I.std(dim='time', skipna=True)
        
        annual_means_by_type_DD_Hg0_I[t_key] = annual_mean_dd_Hg0_i
        annual_std_by_type_DD_Hg0_I[t_key] = annual_std_dd_Hg0_i
        
        print(f"DD_Hg0_I - Annual mean and std calculation for {t_key} completed.")

        # Calcolo della media mensile (DD_Hg0_I)
        
        monthly_mean_dd_Hg0_I = combined_data_DD_Hg0_I.groupby('time.month').mean(dim='time', skipna=True)
        monthly_std_dd_Hg0_I = combined_data_DD_Hg0_I.groupby('time.month').std(dim='time', skipna=True)
        
        monthly_means_by_type_DD_Hg0_I[t_key] = monthly_mean_dd_Hg0_I
        monthly_std_by_type_DD_Hg0_I[t_key] = monthly_std_dd_Hg0_I
        
        print(f"DD_Hg0_I - Monthly mean and std calculation for {t_key} completed.")
    else:
        print(f"DD_Hg0_I - Nessun dato trovato per {t_key}.")


##########################################
### 
##########################################

# Create a new dictionary to store the converted data
annual_means_by_type_DD_micro = {}
monthly_means_by_type_DD_micro = {}
annual_std_by_type_DD_micro = {}
monthly_std_by_type_DD_micro = {}

annual_means_by_type_DD_micro_Hg2 = {}
monthly_means_by_type_DD_micro_Hg2 = {}
annual_std_by_type_DD_micro_Hg2 = {}
monthly_std_by_type_DD_micro_Hg2 = {}

annual_means_by_type_DD_micro_HgP = {}
monthly_means_by_type_DD_micro_HgP = {}
annual_std_by_type_DD_micro_HgP = {}
monthly_std_by_type_DD_micro_HgP = {}

annual_means_by_type_DD_micro_Hg0 = {}
monthly_means_by_type_DD_micro_Hg0 = {}
annual_std_by_type_DD_micro_Hg0 = {}
monthly_std_by_type_DD_micro_Hg0 = {}

###ISLAND

annual_means_by_type_DD_micro_I = {}
monthly_means_by_type_DD_micro_I = {}
annual_std_by_type_DD_micro_I = {}
monthly_std_by_type_DD_micro_I = {}

annual_means_by_type_DD_micro_Hg2_I = {}
monthly_means_by_type_DD_micro_Hg2_I = {}
annual_std_by_type_DD_micro_Hg2_I = {}
monthly_std_by_type_DD_micro_Hg2_I = {}

annual_means_by_type_DD_micro_HgP_I = {}
monthly_means_by_type_DD_micro_HgP_I = {}
annual_std_by_type_DD_micro_HgP_I = {}
monthly_std_by_type_DD_micro_HgP_I = {}

annual_means_by_type_DD_micro_Hg0_I = {}
monthly_means_by_type_DD_micro_Hg0_I = {}
annual_std_by_type_DD_micro_Hg0_I = {}
monthly_std_by_type_DD_micro_Hg0_I = {}

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD.items():
    
    annual_means_by_type_DD_micro[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD.items():
    
    monthly_means_by_type_DD_micro[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD.items():
    
    annual_std_by_type_DD_micro[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD.items():
    
    monthly_std_by_type_DD_micro[t_key] = data_array * cf_units_dd_year

print("DD TOTAL - Unit conversion completed.")

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD_Hg2.items():
    
    annual_means_by_type_DD_micro_Hg2[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD_Hg2.items():
    
    monthly_means_by_type_DD_micro_Hg2[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD_Hg2.items():
    
    annual_std_by_type_DD_micro_Hg2[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD_Hg2.items():
    
    monthly_std_by_type_DD_micro_Hg2[t_key] = data_array * cf_units_dd_year

print("DD Hg2 - Unit conversion completed.")

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD_HgP.items():
    
    annual_means_by_type_DD_micro_HgP[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD_HgP.items():
    
    monthly_means_by_type_DD_micro_HgP[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD_HgP.items():
    
    annual_std_by_type_DD_micro_HgP[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD_HgP.items():
    
    monthly_std_by_type_DD_micro_HgP[t_key] = data_array * cf_units_dd_year

print("DD HgP - Unit conversion completed.")

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD_Hg0.items():
    
    annual_means_by_type_DD_micro_Hg0[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD_Hg0.items():
    
    monthly_means_by_type_DD_micro_Hg0[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD_Hg0.items():
    
    annual_std_by_type_DD_micro_Hg0[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD_Hg0.items():
    
    monthly_std_by_type_DD_micro_Hg0[t_key] = data_array * cf_units_dd_year

print("DD Hg0 - Unit conversion completed.")


#####ISLAND

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD_I.items():
    
    annual_means_by_type_DD_micro_I[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD_I.items():
    
    monthly_means_by_type_DD_micro_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD_I.items():
    
    annual_std_by_type_DD_micro_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD_I.items():
    
    monthly_std_by_type_DD_micro_I[t_key] = data_array * cf_units_dd_year

print("DD_I TOTAL - Unit conversion completed.")

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD_Hg2_I.items():
    
    annual_means_by_type_DD_micro_Hg2_I[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD_Hg2_I.items():
    
    monthly_means_by_type_DD_micro_Hg2_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD_Hg2_I.items():
    
    annual_std_by_type_DD_micro_Hg2_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD_Hg2_I.items():
    
    monthly_std_by_type_DD_micro_Hg2_I[t_key] = data_array * cf_units_dd_year

print("DD_I Hg2 - Unit conversion completed.")

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD_HgP_I.items():
    
    annual_means_by_type_DD_micro_HgP_I[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD_HgP_I.items():
    
    monthly_means_by_type_DD_micro_HgP_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD_HgP_I.items():
    
    annual_std_by_type_DD_micro_HgP_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD_HgP_I.items():
    
    monthly_std_by_type_DD_micro_HgP_I[t_key] = data_array * cf_units_dd_year

print("DD_I HgP - Unit conversion completed.")

# Iterate through the dictionary and multiply each DataArray by the conversion factor

for t_key, data_array in annual_means_by_type_DD_Hg0_I.items():
    
    annual_means_by_type_DD_micro_Hg0_I[t_key] = data_array * cf_units_dd_year

for t_key, data_array in monthly_means_by_type_DD_Hg0_I.items():
    
    monthly_means_by_type_DD_micro_Hg0_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in annual_std_by_type_DD_Hg0_I.items():
    
    annual_std_by_type_DD_micro_Hg0_I[t_key] = data_array * cf_units_dd_year
    
for t_key, data_array in monthly_std_by_type_DD_Hg0_I.items():
    
    monthly_std_by_type_DD_micro_Hg0_I[t_key] = data_array * cf_units_dd_year

print("DD_I Hg0 - Unit conversion completed.")


##########################################
### PLOTTING STANDARD things 
##########################################

import pandas as pd
import matplotlib.pyplot as plt

# --- 1. Nuove Variabili di Selezione ---
# Coordinate dell'Isola Macquarie
target_lat = -54.5
target_lon = 158.95

# La chiave specifica che hai scelto
selected_key = 'OX' 

# Una lista dei tuoi dizionari (come definiti in precedenza)
list_of_data_dicts = [
    monthly_means_by_type_DD_micro_Hg0,
    monthly_means_by_type_DD_micro_Hg2,
    monthly_means_by_type_DD_micro_HgP,
    #monthly_means_by_type_DD_micro
]

# Etichette per le serie nel grafico
series_labels = ['DD_Hg0', 'DD_Hg2', 'DD_HgP']#, 'DD_Total']

# Dizionario per contenere le serie temporali finali (month-only)
monthly_data_series = {}

# --- 2. Processo di Estrazione al Punto Specifico ---
for data_dict, label in zip(list_of_data_dicts, series_labels):
    if selected_key in data_dict:
        da = data_dict[selected_key]
        
        # **NUOVO PASSO CRUCIALE:** # Seleziona il punto sulla mappa più vicino alle coordinate fornite.
        # Questo riduce la DataArray da (month, lat, lon) a (month).
        point_series = da.sel(
            lat=target_lat, 
            lon=target_lon, 
            method='nearest' # Trova il punto della griglia più vicino
        )
        
        # 3. Converti la DataArray risultante in pandas.Series e memorizza
        # .to_series() è necessario per inserire i dati nel DataFrame di pandas
        monthly_data_series[label] = point_series.to_series()
    else:
        print(f"Attenzione: La chiave '{selected_key}' non è stata trovata nel dizionario per {label}")

# 4. Combina tutte le serie in un unico DataFrame
df_plot = pd.DataFrame(monthly_data_series)

# Opzionale: Rinomina l'indice 'month' per un'etichetta migliore nel plot
df_plot.index.name = "Mese"


# 5. Crea il Grafico ad Area Stacked
plt.figure(figsize=(12, 6))

ax = df_plot.plot.area(
    stacked=True, # Area Plot Stacked (somma delle aree)
    ax=plt.gca(),
    title=f"Monthly Hg Dry Deposition (2018-2022 period) @ MQI - Perturbation key: {selected_key} (Lat: {target_lat}, Lon: {target_lon})"
)

# 6. Personalizzazione e Visualizzazione
ax.set_xlabel("Month")
ax.set_ylabel("Hg dry deposition (ug * m-2 * year-1)")
ax.legend(title="Perturbation keys", loc='upper right')
ax.grid(axis='y', linestyle='--', alpha=0.6)

# Imposta etichette sull'asse X per i 12 mesi

ax.set_xticks(range(1, 13))
ax.set_xticklabels(month_names)

plt.tight_layout()
plt.show()




##########################################
### Plotting on the island 
##########################################

import pandas as pd
import matplotlib.pyplot as plt

# --- 1. Nuove Variabili di Selezione ---
# Coordinate dell'Isola Macquarie
target_lat = -54.5
target_lon = 158.95

# La chiave specifica che hai scelto
selected_key = 'WM' 

# Una lista dei tuoi dizionari (come definiti in precedenza)
list_of_data_dicts = [
    monthly_means_by_type_DD_micro_Hg0_I,
    monthly_means_by_type_DD_micro_Hg2_I,
    monthly_means_by_type_DD_micro_HgP_I,
    #monthly_means_by_type_DD_micro
]

# Etichette per le serie nel grafico
series_labels = ['DD_Hg0_I', 'DD_Hg2_I', 'DD_HgP_I']#, 'DD_Total']

# Dizionario per contenere le serie temporali finali (month-only)
monthly_data_series = {}

# --- 2. Processo di Estrazione al Punto Specifico ---
for data_dict, label in zip(list_of_data_dicts, series_labels):
    if selected_key in data_dict:
        da = data_dict[selected_key]
        
        # **NUOVO PASSO CRUCIALE:** # Seleziona il punto sulla mappa più vicino alle coordinate fornite.
        # Questo riduce la DataArray da (month, lat, lon) a (month).
        point_series = da.sel(
            lat=target_lat, 
            lon=target_lon, 
            method='nearest' # Trova il punto della griglia più vicino
        )
        
        # 3. Converti la DataArray risultante in pandas.Series e memorizza
        # .to_series() è necessario per inserire i dati nel DataFrame di pandas
        monthly_data_series[label] = point_series.to_series()
    else:
        print(f"Attenzione: La chiave '{selected_key}' non è stata trovata nel dizionario per {label}")

# 4. Combina tutte le serie in un unico DataFrame
df_plot = pd.DataFrame(monthly_data_series)

# Opzionale: Rinomina l'indice 'month' per un'etichetta migliore nel plot
df_plot.index.name = "Mese"


# 5. Crea il Grafico ad Area Stacked
plt.figure(figsize=(12, 6))

ax = df_plot.plot.area(
    stacked=True, # Area Plot Stacked (somma delle aree)
    ax=plt.gca(),
    title=f"Monthly Hg Dry Deposition (only 2018) @ MQI - Perturbation key: {selected_key} (Lat: {target_lat}, Lon: {target_lon})"
)

# 6. Personalizzazione e Visualizzazione
ax.set_xlabel("Month")
ax.set_ylabel("Hg dry deposition (ug * m-2 * year-1)")
ax.legend(title="Perturbation keys", loc='upper right')
ax.grid(axis='y', linestyle='--', alpha=0.6)

# Imposta etichette sull'asse X per i 12 mesi

ax.set_xticks(range(1, 13))
ax.set_xticklabels(month_names)

plt.tight_layout()
plt.show()




##########################################
### Hg2 to SSA 
##########################################

data_types = ['HgChem','Met']

SSAloss_dict_x = {}

for y in years:
    
    for t in perturbation_types:
            
            # Build the key dynamically
           
            path_key_Hg2toSSA= f"{perturbation_type_map[t]}_{y}_{data_types[0]}"
            # print(path_key_Hg2toSSA)
            path_key_Met= f"{perturbation_type_map[t]}_{y}_{data_types[1]}"
            # print(path_key_Met)
        
            # Hg2GasToSSA is in molecules/cm3*s, must be converted in ug/m2*year
        
            dataset_Hg2toSSA = datasets_by_year[path_key_Hg2toSSA]['Hg2GasToSSA']
            dataset_met1 = datasets_by_year[path_key_Met]['Met_AIRVOL']
            dataset_met2 = datasets_by_year[path_key_Met]['AREA']
            dataset = (( dataset_met1 * dataset_Hg2toSSA * conv ) / dataset_met2) # << Met_AIRVOL*AREA*Hg2GasToSSA*conv

            result_key = f"{perturbation_type_map[t]}_{y}"
            # Save the dataset in the main dictionary only if it's valid
            if dataset is not None:
                SSAloss_dict_x[result_key] = dataset

# Access the dataset for a specific run
sample_dataset = SSAloss_dict_x['NP_2018']

# Then, you would typically look at its dimensions/coordinates
# If it's an xarray object:
print(sample_dataset.coords)
print(sample_dataset.dims)
print(sample_dataset['lat'].values)

import xarray as xr
import pandas as pd # pandas is often imported for xarray operations

# --- Assumed Input Variables ---
# SSAloss_dict_x: Your dictionary containing DataArrays, e.g., {'NP_2018': <DataArray>, 'OL_2018': <DataArray>, ...}
# perturbation_types: The unique prefixes (scenarios) deduced from your keys.
perturbation_types = ['NP', 'OL', 'OM', 'SSAC', 'SSA0', 'WL', 'WM', 'OX']
# -------------------------------

# Dictionary to store the final 5-year averages for each scenario
SSAloss_5year_avg = {}

for p_type in perturbation_types:
    # 1. Gather keys for the current perturbation type across all years
    keys_for_p_type = [k for k in SSAloss_dict_x.keys() if k.startswith(p_type)]
    
    # 2. Extract DataArrays and Concatenate along the 'time' dimension
    datasets_to_combine = [SSAloss_dict_x[k] for k in keys_for_p_type]
    
    
    # This creates a single DataArray (e.g., 60 months for 5 years)
    combined_dataset = xr.concat(datasets_to_combine, dim='time')
    
    # 3. Calculate Time-Weighted Average (Correct for month length)

    # 3a. Get the length of each month (e.g., 31, 28, 31, ...)
    # xarray's .dt accessor pulls this from the 'time' coordinate
    month_length = combined_dataset.time.dt.days_in_month

    # 3b. Normalize the month lengths to create weights that sum to 1
    # We use .where() to avoid weights for any months that might be missing (NaN)
    weights = month_length.where(~combined_dataset.isnull())
    weights /= weights.sum(dim='time')
    
    # 3c. Apply the weights and sum along the 'time' dimension
    # (Data * Weights) summed is mathematically equivalent to the weighted mean.
    five_year_avg = (combined_dataset * weights).sum(dim='time')
    
    # 4. Store the result with the perturbation type as the key
    SSAloss_5year_avg[p_type] = five_year_avg

print("Calculation Complete. SSAloss_5year_avg now holds the 5-year weighted average for each scenario.")
# Example of how to inspect a result:
# print(SSAloss_5year_avg['NP'])

import xarray as xr
import pandas as pd
import numpy as np

# --- ASSUMED INPUT VARIABLES ---
# SSAloss_dict_x: Your dictionary containing DataArrays, e.g., 
# {'NP_2018': <DataArray>, 'OL_2018': <DataArray>, ...}
# Each DataArray is assumed to have a 'time' dimension.
# -------------------------------

# The unique prefixes (scenarios)
perturbation_types = ['NP', 'OL', 'OM', 'SSAC', 'SSA0', 'WL', 'WM', 'OX']

# Dictionary to store the final 12-month climatology for each scenario
SSAloss_monthly_climatology = {}

print("Starting calculation of 5-year monthly climatology...")

for p_type in perturbation_types:
    # 1. Gather keys for the current perturbation type across all years
    keys_for_p_type = [k for k in SSAloss_dict_x.keys() if k.startswith(p_type)]
    
    # 2. Extract DataArrays and Concatenate along the 'time' dimension
    datasets_to_combine = [SSAloss_dict_x[k] for k in keys_for_p_type]
    
    # This creates a single, long DataArray (e.g., 60 months for 5 years)
    combined_dataset = xr.concat(datasets_to_combine, dim='time')
    
    # 3. Calculate the Climatological Monthly Average
    
    # We use .groupby('time.month') to group all 60 time points into 12 groups (Jan, Feb, ..., Dec).
    # Then we use .mean() to average the data within those 5 instances of each month.
    # The result has a new dimension called 'month' (size 12).
    monthly_climatology = combined_dataset.groupby('time.month').mean(dim='time')
    
    # 4. Store the result
    SSAloss_monthly_climatology[p_type] = monthly_climatology
    print(f"Climatology calculated for scenario: {p_type}")

print("\nCalculation Complete.")
print("SSAloss_monthly_climatology now holds the 8 DataArrays needed for time series plotting.")
# Example of how to inspect the dimensions of a result:
# print(SSAloss_monthly_climatology['NP'])

# --- MOCKUP FOR DEMONSTRATION PURPOSES ---
# If your SSAloss_dict_x is not defined, this block runs a small example
if 'SSAloss_dict_x' not in globals():
    print("\n--- Running Mockup Example ---")
    
    # Create mock time coordinates for 5 years (60 months)
    time_coords = pd.date_range('2018-01-15', periods=60, freq='M')
    
    # Mock spatial coordinates
    mock_lev = np.array([1000, 500])
    mock_lat = np.array([45, 55])
    mock_lon = np.array([0, 10])
    
    # Create a simple, synthetic SSA loss value that varies seasonally
    seasonal_data = np.sin(np.linspace(0, 2*np.pi * 5, 60))[:, np.newaxis, np.newaxis, np.newaxis]
    mock_data = (seasonal_data + 1) * 1e-4
    
    # Create mock dictionary structure
    mock_dict = {}
    for year in range(2018, 2023):
        start_idx = (year - 2018) * 12
        end_idx = start_idx + 12
        
        # Mock DataArray for a single year (e.g., 'NP_2018')
        mock_da = xr.DataArray(
            mock_data[start_idx:end_idx, :, :, :].squeeze(), # squeeze removes the unnecessary dimension
            coords={
                'time': time_coords[start_idx:end_idx],
                'lev': mock_lev, 
                'lat': mock_lat, 
                'lon': mock_lon
            },
            dims=['time', 'lev', 'lat', 'lon'],
            name='SSAloss'
        )
        mock_dict[f'NP_{year}'] = mock_da
    
    # Run the calculation with mock data
    SSAloss_dict_x = mock_dict
    
    # Re-run the core logic to calculate climatology on the mock data
    mock_climatology = {}
    
    # Simplified loop for mock example
    datasets_to_combine = [SSAloss_dict_x[k] for k in SSAloss_dict_x.keys()]
    combined_dataset = xr.concat(datasets_to_combine, dim='time')
    mock_climatology['NP'] = combined_dataset.groupby('time.month').mean(dim='time')
    
    print("\n--- Mock Climatology Result (NP) ---")
    print(mock_climatology['NP'])
    print("Dimensions:", mock_climatology['NP'].dims) # Should show ('month', 'lev', 'lat', 'lon')


# Access the dataset for a specific run
sample_dataset = SSAloss_5year_avg['OX']

# Then, you would typically look at its dimensions/coordinates
# If it's an xarray object:
print(sample_dataset.coords)
print(sample_dataset.dims)
print(sample_dataset['lat'].values)

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

# --- NOTE TO USER: SETUP ASSUMPTION ---
# This code assumes you have already defined the 'SSAloss_5year_avg' 
# dictionary in your Python session, containing the 5-year averages 
# for each scenario, including 'NP'.
# -------------------------------------

def plot_ssa_loss_map(ssa_loss_dict, case_name='NP', dim_to_sum='lev', title_suffix="5-Year Average (2018-2022)"):
    """
    Generates a global map for a specified case, calculated as the SUM (Total) 
    across the vertical dimension.
    
    Args:
        ssa_loss_dict (dict): Dictionary containing 3D xarray DataArrays (lev, lat, lon).
        case_name (str): The key in the dictionary to plot (e.g., 'NP').
        dim_to_sum (str): The dimension to sum over to get a 2D map (e.g., 'lev').
        title_suffix (str): Text to append to the plot title.
    """
    
    if case_name not in ssa_loss_dict:
        print(f"Error: Case '{case_name}' not found in the dictionary.")
        return

    # 1. Select the 3D data (lev, lat, lon)
    data_3d = ssa_loss_dict[case_name]
    
    # 2. Calculate the SUM across the vertical dimension ('lev') to get a 2D map (lat, lon)
    # The result represents the TOTAL integrated SSA loss throughout the entire atmospheric column.
    data_2d = data_3d.sum(dim=dim_to_sum) # CHANGED FROM .mean() TO .sum()

    # Determine plot properties
    plot_title = f"{case_name} Sea Salt Uptake - Column Total (Summed Over Level) {title_suffix}"
    
    # Set up the figure and the GeoAxes
    fig = plt.figure(figsize=(12, 8))
    
    # Define the projection for the map display
    # PlateCarree is a simple, rectangular global map projection
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    # Add map features for context
    ax.coastlines(resolution='50m', color='gray')
    ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    
    # 3. Plot the data using xarray's built-in plotting method (pcolormesh)
    # The 'transform' argument tells cartopy what projection the data coordinates are in.
    # The data is already in standard lat/lon, so we use ccrs.PlateCarree().
    plot = data_2d.plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='viridis',  # You can choose any matplotlib colormap
        cbar_kwargs={
            'label': 'ug/m2*year', # Updated label
            'shrink': 0.8
        },
        extend='both', # Extend the colorbar for values outside the normal range
        # You might need to adjust vmin/vmax based on your actual data range
        # vmin=1e-5,
        # vmax=1e-3, 
        
    )

    ax.set_title(plot_title)
    
    # Optional: Set the extent of the map (full globe)
    ax.set_global() 
    
    plt.show()

# --- Placeholder for the execution call ---
# You need to run this function with your actual dictionary.
# Since I cannot access your environment, I will show the call here.
# Assuming SSAloss_5year_avg exists:
# plot_ssa_loss_map(SSAloss_5year_avg, case_name='NP')

# --- MOCKUP DATA FOR RUNNABILITY (Replace with your actual execution) ---
# Since this code must be runnable, I will generate mock data that matches your structure
# If you run this file directly, it will use the mock data.
if 'SSAloss_5year_avg' not in globals():
    print("\n--- Generating MOCK Data for Demo ---")
    
    # Create mock coordinates
    mock_lev = np.linspace(0, 1, 10)
    mock_lat = np.arange(-89, 90, 4)
    mock_lon = np.arange(-180, 180, 5)
    
    # Create mock data array (3D: lev, lat, lon)
    # This mock data is designed to show a higher loss in the Arctic (high lat)
    lon_grid, lat_grid = np.meshgrid(mock_lon, mock_lat)
    mock_values = (1 + np.sin(np.deg2rad(lat_grid) * 3)) * np.cos(np.deg2rad(lon_grid) / 5)
    
    # Add a level dimension to the mock data and scale it
    mock_data = mock_values[np.newaxis, :, :] * np.exp(-mock_lev[:, np.newaxis, np.newaxis] * 2) * 1e-4
    
    # Create the mock DataArray
    mock_da = xr.DataArray(
        mock_data,
        coords={'lev': mock_lev, 'lat': mock_lat, 'lon': mock_lon},
        dims=['lev', 'lat', 'lon'],
        name='SSAloss'
    )
    
    # Create the mock SSAloss_5year_avg dictionary
    mock_SSAloss_5year_avg = {'NP': mock_da}
    
    # Execute the plotting function with the mock data
    plot_ssa_loss_map(mock_SSAloss_5year_avg, case_name='NP')
else:
    # Execute the plotting function with the real data
    plot_ssa_loss_map(SSAloss_5year_avg, case_name='OX')

##### Other plot 

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import calendar

def plot_climatological_timeseries(climatology_dict, case_name='NP', plot_title="Climatological SSA Loss Cycle (Global Mean, Column Summed)"):
    """
    Plots the 12-month climatological cycle for a specified scenario.
    
    The data is spatially averaged (lat/lon) and summed vertically (lev) 
    to represent the total global atmospheric impact over the year.
    
    Args:
        climatology_dict (dict): Dictionary containing the 4D monthly climatology DataArrays.
        case_name (str): The key in the dictionary to plot (e.g., 'NP').
        plot_title (str): The title for the resulting plot.
    """
    
    if case_name not in climatology_dict:
        print(f"Error: Case '{case_name}' not found in the climatology dictionary.")
        return

    # 1. Select the 4D climatology data (month, lev, lat, lon)
    data_4d = climatology_dict[case_name]
    
    # 2. Reduce the data to a 1D time series (size 12)
    # a. Sum across the vertical dimension ('lev') to get the Column Total.
    data_column_sum = data_4d.sum(dim='lev')
    
    # b. Calculate the Global Mean (mean across 'lat' and 'lon').
    # The final result is a 1D array indexed by 'month'.
    time_series_data = data_column_sum.mean(dim=['lat', 'lon'])
    
    # 3. Prepare plot coordinates
    months = time_series_data['month'].values
    month_names = [calendar.month_abbr[m] for m in months]
    
    # 4. Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot the time series
    ax.plot(months, time_series_data.values, marker='o', linestyle='-', color='indigo', linewidth=2)
    
    # Set X-axis to display month names
    ax.set_xticks(months)
    ax.set_xticklabels(month_names)
    
    # Labeling and Titling
    ax.set_title(f"{case_name} {plot_title}", fontsize=14)
    ax.set_xlabel("Month of the Year", fontsize=12)
    ax.set_ylabel("Total Global SSA Loss (Column Summed)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Improve aesthetics
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.tight_layout()
    plt.show()

# --- MOCKUP DATA FOR RUNNABILITY (Replace with your actual execution) ---
# If your SSAloss_monthly_climatology is not defined, this block runs a small example
if 'SSAloss_monthly_climatology' not in globals():
    print("\n--- Running Mockup Example for Plotting ---")
    
    # Create mock 4D climatology data (12 months, 2 lev, 2 lat, 2 lon)
    mock_months = np.arange(1, 13)
    mock_lev = np.array([1000, 500])
    mock_lat = np.array([45, 55])
    mock_lon = np.array([0, 10])
    
    # Create synthetic seasonal data (peaks in summer, dips in winter)
    seasonal_data = (np.cos(np.linspace(0, 2 * np.pi, 12) + np.pi/2) * 0.5 + 1.5) * 1e-4
    
    # Expand data to 4D structure
    mock_4d_data = seasonal_data[:, np.newaxis, np.newaxis, np.newaxis] * np.ones((12, 2, 2, 2))
    
    mock_da = xr.DataArray(
        mock_4d_data,
        coords={'month': mock_months, 'lev': mock_lev, 'lat': mock_lat, 'lon': mock_lon},
        dims=['month', 'lev', 'lat', 'lon'],
        name='SSAloss'
    )
    
    # Create the mock climatology dictionary
    mock_climatology = {'NP': mock_da}
    
    # Execute the plotting function with the mock data
    plot_climatological_timeseries(mock_climatology, case_name='NP')

else:
    # Execute the plotting function with the real data
    plot_climatological_timeseries(SSAloss_monthly_climatology, case_name='OX')


##########################################
### SEA SALT UPTAKE COMPARISON - GLOBAL 
##########################################

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import calendar

# --- NOTE TO USER: SETUP ASSUMPTION ---
# This code assumes the 'SSAloss_monthly_climatology' dictionary has 
# already been created using the 'calculate_climatology.py' script 
# (4D DataArrays: month, lev, lat, lon).
# -------------------------------------

# Define all scenarios
ALL_PERTURBATION_TYPES = ['NP', 'OL', 'OM', 'SSAC', 'SSA0', 'WL', 'WM', 'OX']



def plot_climatological_timeseries_comparison(climatology_dict, plot_title="Sea Salt Uptake - Comparison (Global Mean (2018-2022), Column Summed)"):
    """
    Plots the 12-month climatological cycle for ALL specified scenarios 
    on a single graph for comparison, using custom colors and linestyles.
    
    The data is spatially averaged (lat/lon) and summed vertically (lev) 
    to represent the total global atmospheric impact over the year.
    
    Args:
        climatology_dict (dict): Dictionary containing the 4D monthly climatology DataArrays.
        plot_title (str): The title for the resulting plot.
    """
    
    if not all(p in climatology_dict for p in ALL_PERTURBATION_TYPES):
        missing = [p for p in ALL_PERTURBATION_TYPES if p not in climatology_dict]
        print(f"Error: Missing cases in dictionary: {missing}. Cannot plot all series.")
        return

    # 1. Setup the plot
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Get month coordinates (same for all cases)
    first_case = ALL_PERTURBATION_TYPES[0]
    months = climatology_dict[first_case]['month'].values
    
    # Use the global 1-character month names as requested
   # month_names = month_names

    # 2. Loop through all perturbation types and plot the series
    for p_type in ALL_PERTURBATION_TYPES:
        data_4d = climatology_dict[p_type]
        
        # a. Sum across the vertical dimension ('lev') to get the Column Total.
        data_column_sum = data_4d.sum(dim='lev')
        
        # b. Calculate the Global Mean (mean across 'lat' and 'lon').
        time_series_data = data_column_sum.mean(dim=['lat', 'lon'])
        
        # c. Plot the series using custom colors, linestyles, and legend labels
        ax.plot(
            months, 
            time_series_data.values, 
            marker='o', 
            linewidth=2,
            color=colors.get(p_type, 'gray'),             # Use new colors dictionary
            linestyle=linestyles.get(p_type, '-'),        # Use new linestyles dictionary
            label=legend_map.get(p_type, p_type)          # Use new legend_map dictionary
        )
    
    # 3. Final Plot Customization
    
    # Set X-axis to display 1-character month names
    ax.set_xticks(months)
    ax.set_xticklabels(month_names)
    
    # Labeling and Titling
    ax.set_title(plot_title, fontsize=16, fontweight='bold')
    ax.set_xlabel("Month of the Year", fontsize=12)
    ax.set_ylabel("Total Global Sea Salt Uptake (ug/m2*year)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Add Legend to distinguish the lines
    ax.legend(title="Scenario", frameon=True, shadow=True, loc='best')
    
    # Improve aesthetics
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.tight_layout()
    plt.show()

# --- MOCKUP DATA FOR RUNNABILITY (Replace with your actual execution) ---
# If your SSAloss_monthly_climatology is not defined, this block runs a small example
if 'SSAloss_monthly_climatology' not in globals():
    print("\n--- Running Mockup Example for Plotting All Series ---")
    
    # Create mock 4D climatology data (12 months, 2 lev, 2 lat, 2 lon)
    mock_months = np.arange(1, 13)
    mock_lev = np.array([1000, 500])
    mock_lat = np.array([45, 55])
    mock_lon = np.array([0, 10])
    
    # Create synthetic seasonal data (peaks in summer, dips in winter)
    base_seasonal = (np.cos(np.linspace(0, 2 * np.pi, 12) + np.pi/2) * 0.5 + 1.5)
    
    # Create the mock climatology dictionary with varied scales
    mock_climatology = {}
    
    for i, p_type in enumerate(ALL_PERTURBATION_TYPES):
        # Scale each type slightly differently
        scale = 1e-4 * (1 + i * 0.1) 
        
        # Expand data to 4D structure
        mock_4d_data = (base_seasonal * scale)[:, np.newaxis, np.newaxis, np.newaxis] * np.ones((12, 2, 2, 2))
        
        mock_da = xr.DataArray(
            mock_4d_data,
            coords={'month': mock_months, 'lev': mock_lev, 'lat': mock_lat, 'lon': mock_lon},
            dims=['month', 'lev', 'lat', 'lon'],
            name='SSAloss'
        )
        mock_climatology[p_type] = mock_da
    
    # Execute the plotting function with the mock data
    plot_climatological_timeseries_comparison(mock_climatology)

else:
    # Execute the plotting function with the real data
    plot_climatological_timeseries_comparison(SSAloss_monthly_climatology)


##########################################
### SEA SALT UPTAKE COMPARISON - ISLAND 
##########################################

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import calendar

# --- GEOGRAPHIC FOCUS ---
# Coordinates for Macquarie Island (as requested)
TARGET_LAT = -54.5
TARGET_LON = 158.95
# ------------------------

# --- NOTE TO USER: SETUP ASSUMPTION ---
# This code assumes the 'SSAloss_monthly_climatology' dictionary has 
# already been created using the 'calculate_climatology.py' script 
# (4D DataArrays: month, lev, lat, lon).
# -------------------------------------

# Define all scenarios
ALL_PERTURBATION_TYPES = ['NP', 'OL', 'OM', 'SSAC', 'SSA0', 'WL', 'WM', 'OX']


def plot_climatological_timeseries_comparison(climatology_dict):
    """
    Plots the 12-month climatological cycle for ALL specified scenarios 
    on a single graph, focused on the location defined by TARGET_LAT/LON.
    
    The data is summed vertically (lev) to represent the total atmospheric 
    impact at the single selected grid point.
    
    Args:
        climatology_dict (dict): Dictionary containing the 4D monthly climatology DataArrays.
    """
    
    if not all(p in climatology_dict for p in ALL_PERTURBATION_TYPES):
        missing = [p for p in ALL_PERTURBATION_TYPES if p not in climatology_dict]
        print(f"Error: Missing cases in dictionary: {missing}. Cannot plot all series.")
        return

    # 1. Setup the plot
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Get month coordinates (same for all cases)
    first_case = ALL_PERTURBATION_TYPES[0]
    months = climatology_dict[first_case]['month'].values
   # month_names = MONTH_ABBREVIATIONS
    
    plot_title = f"Sea Salt Uptake - Comparison at Lat: {TARGET_LAT}, Lon: {TARGET_LON}"

    # 2. Loop through all perturbation types and plot the series
    for p_type in ALL_PERTURBATION_TYPES:
        data_4d = climatology_dict[p_type]
        
        # a. Select the single grid cell closest to the target coordinates
        # Use method='nearest' to ensure we pick the closest grid point
        data_point = data_4d.sel(
            lat=TARGET_LAT, 
            lon=TARGET_LON, 
            method='nearest'
        )
        
        # b. Sum across the vertical dimension ('lev') to get the Column Total at that point
        # The result is a 1D array indexed by 'month'.
        time_series_data = data_point.sum(dim='lev')
        
        # c. Plot the series using custom colors, linestyles, and legend labels
        ax.plot(
            months, 
            time_series_data.values, 
            marker='o', 
            linewidth=2,
            color=colors.get(p_type, 'gray'),             
            linestyle=linestyles.get(p_type, '-'),        
            label=legend_map.get(p_type, p_type)          
        )
    
    # 3. Final Plot Customization
    
    # Set X-axis to display 1-character month names
    ax.set_xticks(months)
    ax.set_xticklabels(month_names)
    
    # Labeling and Titling
    ax.set_title(plot_title, fontsize=16, fontweight='bold')
    ax.set_xlabel("Month of the Year", fontsize=12)
    ax.set_ylabel("Sea Salt Uptake (ug/m2*year)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Add Legend to distinguish the lines
    ax.legend(title="Scenario", frameon=True, shadow=True, loc='best')
    
    # Improve aesthetics
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.tight_layout()
    plt.show()

# --- MOCKUP DATA FOR RUNNABILITY (Replace with your actual execution) ---
# If your SSAloss_monthly_climatology is not defined, this block runs a small example
if 'SSAloss_monthly_climatology' not in globals():
    print(f"\n--- Running Mockup Example Focused on Lat: {TARGET_LAT}, Lon: {TARGET_LON} ---")
    
    # Create mock 4D climatology data (12 months, 2 lev, 2 lat, 2 lon)
    mock_months = np.arange(1, 13)
    mock_lev = np.array([1000, 500])
    
    # Use coordinates near the target for the mock data
    mock_lat = np.array([-54, -55])
    mock_lon = np.array([158, 159])
    
    # Create synthetic seasonal data (peaks in Southern Hemisphere summer, Dec-Feb)
    # Shifted cosine function to peak around month 1 (Jan)
    base_seasonal = (np.cos(np.linspace(0, 2 * np.pi, 12) + 2.6) * 0.5 + 1.5)
    
    # Create the mock climatology dictionary with varied scales
    mock_climatology = {}
    
    for i, p_type in enumerate(ALL_PERTURBATION_TYPES):
        # Scale each type slightly differently
        scale = 1e-4 * (1 + i * 0.1) 
        
        # Expand data to 4D structure
        # Add a slight spatial gradient so the selection is meaningful
        spatial_gradient = np.array([[[1.1, 1.0], [0.9, 0.8]]]) * scale
        mock_4d_data = base_seasonal[:, np.newaxis, np.newaxis, np.newaxis] * spatial_gradient
        
        mock_da = xr.DataArray(
            mock_4d_data,
            coords={'month': mock_months, 'lev': mock_lev, 'lat': mock_lat, 'lon': mock_lon},
            dims=['month', 'lev', 'lat', 'lon'],
            name='SSAloss'
        )
        mock_climatology[p_type] = mock_da
    
    # Execute the plotting function with the mock data
    plot_climatological_timeseries_comparison(mock_climatology)

else:
    # Execute the plotting function with the real data
    plot_climatological_timeseries_comparison(SSAloss_monthly_climatology)




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
