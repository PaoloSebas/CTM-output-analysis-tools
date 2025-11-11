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
