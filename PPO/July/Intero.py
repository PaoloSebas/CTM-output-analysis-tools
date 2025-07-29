import xarray as xr
print('XARRAY version:', xr.__version__)
import netCDF4
import numpy as np
print('NUMPY version:', np.__version__)
import datetime # For handling dates
import pandas as pd
print('PANDA version:', pd.__version__)
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import seaborn as sns
import cartopy.feature as cfeature
import os

# --- Configuration ---
# Define base directories for different scenarios
base_directory_s1 = 'D:/PPO/NEW_SET_2025/SCENARIO_1/'
base_directory_s2 = 'D:/PPO/NEW_SET_2025/SCENARIO_2/'
base_directory_s3 = 'D:/PPO/NEW_SET_2025/SCENARIO_3/'

# Define base filenames for species concentration files for each scenario
species_conc_base_filename_s1 = 'SCENARIO1_GEOSChem.SpeciesConc.'
species_conc_base_filename_s2 = 'SCENARIO2_GEOSChem.SpeciesConc.'
species_conc_base_filename_s3 = 'SCENARIO3_GEOSChem.SpeciesConc.'

# Define base filenames for state met files for each scenario
state_met_base_filename_s1 = 'SCENARIO1_GEOSChem.StateMet.'
state_met_base_filename_s2 = 'SCENARIO2_GEOSChem.StateMet.'
state_met_base_filename_s3 = 'SCENARIO3_GEOSChem.StateMet.'

# Define the start and end dates for the analysis
start_date_str = '20190801'
end_date_str = '20200801'

# --- Processing Each Scenario and File Type ---

# Dictionary to store the concatenated datasets for easier access
concatenated_datasets = {}

###############################################################################

# --- Scenario 1 ---
print("Processing SCENARIO 1...")
# Species Concentration files S1
species_conc_pattern_s1 = os.path.join(base_directory_s1, species_conc_base_filename_s1 + '*.nc4')
multi_ds_species_s1 = xr.open_mfdataset(
    species_conc_pattern_s1,
    combine='by_coords',
    compat='override',
    data_vars='minimal',
    coords='minimal',
    parallel=True # Keep parallel=True if dask is installed and working, otherwise remove
)
concatenated_datasets['s1_species'] = multi_ds_species_s1
print(f"  SCENARIO 1 SpeciesConc loaded. Dimensions: {multi_ds_species_s1.dims}")

#multi_ds_species_s1_ppb = multi_ds_species_s1 * 1.0e9
#multi_ds_species_s1_ppb.attrs['units'] = 'ppb'
#print(multi_ds_species_s1_ppb.attrs['units'])

# StatMet files S1
StateMet_pattern_s1 = os.path.join(base_directory_s1, state_met_base_filename_s1 + '*.nc4')
multi_ds_StateMet_s1 = xr.open_mfdataset(
    StateMet_pattern_s1,
    combine='by_coords',
    compat='override',
    data_vars='minimal',
    coords='minimal',
    parallel=True # Keep parallel=True if dask is installed and working, otherwise remove
)
concatenated_datasets['s1_StateMet'] = multi_ds_StateMet_s1
print(f"  SCENARIO 1 StateMet loaded. Dimensions: {multi_ds_StateMet_s1.dims}")

##############################################################################################

# --- Scenario 2 ---
print("Processing SCENARIO 2...")
# Species Concentration files
species_conc_pattern_s2 = os.path.join(base_directory_s2, species_conc_base_filename_s2 + '*.nc4')
multi_ds_species_s2 = xr.open_mfdataset(
    species_conc_pattern_s2,
    combine='by_coords',
    compat='override',
    data_vars='minimal',
    coords='minimal',
    parallel=True # Keep parallel=True if dask is installed and working, otherwise remove
)
concatenated_datasets['s2_species'] = multi_ds_species_s2
print(f"  SCENARIO 2 SpeciesConc loaded. Dimensions: {multi_ds_species_s2.dims}")

#multi_ds_species_s2_ppb = multi_ds_species_s2 * 1.0e9
#multi_ds_species_s2_ppb.attrs['units'] = 'ppbv'
#print(multi_ds_species_s2_ppb.attrs['units'])

# StatMet files S2
StateMet_pattern_s2 = os.path.join(base_directory_s2, state_met_base_filename_s2 + '*.nc4')
multi_ds_StateMet_s2 = xr.open_mfdataset(
    StateMet_pattern_s2,
    combine='by_coords',
    compat='override',
    data_vars='minimal',
    coords='minimal',
    parallel=True # Keep parallel=True if dask is installed and working, otherwise remove
)
concatenated_datasets['s2_StateMet'] = multi_ds_StateMet_s2
print(f"  SCENARIO 2 StateMet loaded. Dimensions: {multi_ds_StateMet_s2.dims}")

##############################################################################################

# --- Scenario 3 ---
print("Processing SCENARIO 3...")
# Species Concentration files S3
species_conc_pattern_s3 = os.path.join(base_directory_s3, species_conc_base_filename_s3 + '*.nc4')
multi_ds_species_s3 = xr.open_mfdataset(
    species_conc_pattern_s3,
    combine='by_coords',
    compat='override',
    data_vars='minimal',
    coords='minimal',
    parallel=True # Keep parallel=True if dask is installed and working, otherwise remove
)
concatenated_datasets['s3_species'] = multi_ds_species_s3
print(f"  SCENARIO 3 SpeciesConc loaded. Dimensions: {multi_ds_species_s3.dims}")

#multi_ds_species_s3_ppb = multi_ds_species_s3 * 1.0e9
#multi_ds_species_s3_ppb.attrs['units'] = 'ppbv'
#print(multi_ds_species_s3_ppb.attrs['units'])

# StatMet files S3nci
StateMet_pattern_s3 = os.path.join(base_directory_s3, state_met_base_filename_s3 + '*.nc4')
multi_ds_StateMet_s3 = xr.open_mfdataset(
    StateMet_pattern_s3,
    combine='by_coords',
    compat='override',
    data_vars='minimal',
    coords='minimal',
    parallel=True # Keep parallel=True if dask is installed and working, otherwise remove
)
concatenated_datasets['s3_StateMet'] = multi_ds_StateMet_s3
print(f"  SCENARIO 3 StateMet loaded. Dimensions: {multi_ds_StateMet_s3.dims}")

print(concatenated_datasets)
print(multi_ds_species_s1)

#Define variable and datasets name
H2CO = 'SpeciesConcVV_CH2O'
HO2 = 'SpeciesConcVV_HO2'
AMG = 'Met_AD'

### ON ALL LEVELS 
H2CO_X_S1= multi_ds_species_s1[H2CO]
HO2_X_S1 = multi_ds_species_s1[HO2]
Air_mass_grid_S1 = multi_ds_StateMet_s1[AMG]

H2CO_X_S2= multi_ds_species_s2[H2CO]
HO2_X_S2 = multi_ds_species_s2[HO2]
Air_mass_grid_S2 = multi_ds_StateMet_s2[AMG]

H2CO_X_S3= multi_ds_species_s3[H2CO]
HO2_X_S3 = multi_ds_species_s3[HO2]
Air_mass_grid_S3 = multi_ds_StateMet_s3[AMG]

#SCENARIO 1
H2CO_multids_S1_surface = H2CO_X_S1.isel(lev=0)
HO2_multids_S1_surface  = HO2_X_S1.isel(lev=0)

H2CO_multids_S1_middle = H2CO_X_S1.isel(lev=23) ### MIDDLE TROPOSPHERE
HO2_multids_S1_middle  = HO2_X_S1.isel(lev=23)

#SCENARIO 2
H2CO_multids_S2_surface = H2CO_X_S2.isel(lev=0)
HO2_multids_S2_surface = HO2_X_S2.isel(lev=0)

H2CO_multids_S2_middle = H2CO_X_S2.isel(lev=23) ### MIDDLE TROPOSPHERE
HO2_multids_S2_middle = HO2_X_S2.isel(lev=23)

#SCENARIO 3
H2CO_multids_S3_surface = H2CO_X_S3.isel(lev=0)
HO2_multids_S3_surface = HO2_X_S3.isel(lev=0)

H2CO_multids_S3_middle = H2CO_X_S3.isel(lev=23) ### MIDDLE TROPOSPHERE
HO2_multids_S3_middle = HO2_X_S3.isel(lev=23)

#All levels 
H2CO_S1_ppb = H2CO_X_S1 * 1e09
H2CO_S2_ppb = H2CO_X_S2 * 1e09
H2CO_S3_ppb = H2CO_X_S3 * 1e09

H2CO_S1_ppb.attrs['units'] = 'ppb'
H2CO_S2_ppb.attrs['units'] = 'ppb'
H2CO_S3_ppb.attrs['units'] = 'ppb'

HO2_S1_ppt = HO2_X_S1 * 1e12
HO2_S2_ppt = HO2_X_S2 * 1e12
HO2_S3_ppt = HO2_X_S3 * 1e12

HO2_S1_ppt.attrs['units'] = 'ppt'
HO2_S2_ppt.attrs['units'] = 'ppt'
HO2_S3_ppt.attrs['units'] = 'ppt'

# Surface_ppb_ppt
H2CO_S1_ppb_surf = H2CO_multids_S1_surface * 1e09   #ppb
H2CO_S2_ppb_surf = H2CO_multids_S2_surface * 1e09   #ppb
H2CO_S3_ppb_surf = H2CO_multids_S3_surface * 1e09   #ppb

HO2_S1_ppt_surf = HO2_multids_S1_surface * 1e12   #ppt
HO2_S2_ppt_surf = HO2_multids_S2_surface * 1e12   #ppt
HO2_S3_ppt_surf = HO2_multids_S3_surface * 1e12   #ppt

H2CO_S1_ppb_surf.attrs['units'] = 'ppb'
H2CO_S2_ppb_surf.attrs['units'] = 'ppb'
H2CO_S3_ppb_surf.attrs['units'] = 'ppb'

HO2_S1_ppt_surf.attrs['units'] = 'ppt'
HO2_S2_ppt_surf.attrs['units'] = 'ppt'
HO2_S3_ppt_surf.attrs['units'] = 'ppt'

# Middle troposphere_ppb_ppt
H2CO_S1_ppb_middle = H2CO_multids_S1_middle * 1e09   #ppb
H2CO_S2_ppb_middle = H2CO_multids_S2_middle * 1e09   #ppb
H2CO_S3_ppb_middle = H2CO_multids_S3_middle * 1e09   #ppb

HO2_S1_ppt_middle = HO2_multids_S1_middle * 1e12   #ppt
HO2_S2_ppt_middle = HO2_multids_S2_middle * 1e12   #ppt
HO2_S3_ppt_middle = HO2_multids_S3_middle * 1e12   #ppt

H2CO_S1_ppb_middle.attrs['units'] = 'ppb'
H2CO_S2_ppb_middle.attrs['units'] = 'ppb'
H2CO_S3_ppb_middle.attrs['units'] = 'ppb'

HO2_S1_ppt_middle.attrs['units'] = 'ppt'
HO2_S2_ppt_middle.attrs['units'] = 'ppt'
HO2_S3_ppt_middle.attrs['units'] = 'ppt'

# Surf MEAN
H2CO_S1_mean_surf = H2CO_S1_ppb_surf.mean(dim='time')
HO2_S1_mean_surf =  HO2_S1_ppt_surf.mean(dim='time')

H2CO_S2_mean_surf = H2CO_S2_ppb_surf.mean(dim='time')
HO2_S2_mean_surf =  HO2_S2_ppt_surf.mean(dim='time')

H2CO_S3_mean_surf = H2CO_S3_ppb_surf.mean(dim='time')
HO2_S3_mean_surf =  HO2_S3_ppt_surf.mean(dim='time')

H2CO_S1_mean_surf.attrs['units'] = 'ppb'
HO2_S1_mean_surf.attrs['units'] = 'ppt'
H2CO_S2_mean_surf.attrs['units'] = 'ppb'
HO2_S2_mean_surf.attrs['units'] = 'ppt'
H2CO_S3_mean_surf.attrs['units'] = 'ppb'
HO2_S3_mean_surf.attrs['units'] = 'ppt'

# Middle troposphere MEAN

H2CO_S1_mean_middle = H2CO_S1_ppb_middle.mean(dim='time')
HO2_S1_mean_middle =  HO2_S1_ppt_middle.mean(dim='time')

H2CO_S2_mean_middle = H2CO_S2_ppb_middle.mean(dim='time')
HO2_S2_mean_middle =  HO2_S2_ppt_middle.mean(dim='time')

H2CO_S3_mean_middle = H2CO_S3_ppb_middle.mean(dim='time')
HO2_S3_mean_middle =  HO2_S3_ppt_middle.mean(dim='time')

H2CO_S1_mean_middle.attrs['units'] = 'ppb'
HO2_S1_mean_middle.attrs['units'] = 'ppt'
H2CO_S2_mean_middle.attrs['units'] = 'ppb'
HO2_S2_mean_middle.attrs['units'] = 'ppt'
H2CO_S3_mean_middle.attrs['units'] = 'ppb'
HO2_S3_mean_middle.attrs['units'] = 'ppt'

#ABSOLUTE DIFFERENCES
# S2 S1
H2CO_S2_S1_diff_surf = H2CO_S2_mean_surf - H2CO_S1_mean_surf
HO2_S2_S1_diff_surf  = HO2_S2_mean_surf - HO2_S1_mean_surf

H2CO_S2_S1_diff_middle = H2CO_S2_mean_middle - H2CO_S1_mean_middle
HO2_S2_S1_diff_middle  = HO2_S2_mean_middle - HO2_S1_mean_middle

#S3 S1
H2CO_S3_S1_diff_surf = H2CO_S3_mean_surf - H2CO_S1_mean_surf
HO2_S3_S1_diff_surf  = HO2_S3_mean_surf - HO2_S1_mean_surf

H2CO_S3_S1_diff_middle = H2CO_S3_mean_middle - H2CO_S1_mean_middle
HO2_S3_S1_diff_middle  = HO2_S3_mean_middle - HO2_S1_mean_middle

#PERCENT

H2CO_S2_S1_perc_surf   = (H2CO_S2_S1_diff_surf/H2CO_S1_mean_surf) * 100
H2CO_S2_S1_perc_middle = (H2CO_S2_S1_diff_middle/H2CO_S1_mean_middle) * 100

H2CO_S3_S1_perc_surf   = (H2CO_S3_S1_diff_surf/H2CO_S1_mean_surf) * 100
H2CO_S3_S1_perc_middle = (H2CO_S3_S1_diff_middle/H2CO_S1_mean_middle) * 100

HO2_S2_S1_perc_surf   = (HO2_S2_S1_diff_surf/HO2_S1_mean_surf) * 100
HO2_S2_S1_perc_middle = (HO2_S2_S1_diff_middle/HO2_S1_mean_middle) * 100

HO2_S3_S1_perc_surf   = (HO2_S3_S1_diff_surf/HO2_S1_mean_surf) * 100
HO2_S3_S1_perc_middle = (HO2_S3_S1_diff_middle/HO2_S1_mean_middle) * 100



def calcola_burden(
    mr_species: xr.DataArray,
    met_ad: xr.DataArray,
    mw_species_g_per_mol: float,
    mw_dry_air_g_per_mol: float = 28.97
) -> float:
    """
    Calcola il burden totale (massa totale in kg) di una specie atmosferica,
    sommando su ( lev, lat, lon) o su tutte le dimensioni.
    Questa funzione assume che i DataArray di input siano già stati caricati
    e allineati (se necessario) e che abbiano le dimensioni corrette.

    Args:
        mr_species_da (xr.DataArray): DataArray di xarray contenente il mixing ratio della specie.
                                      Dovrebbe includere gli attributi 'units' (es. 'mol mol-1', 'ppbv').
        met_ad_da (xr.DataArray): DataArray di xarray contenente la massa dell'aria secca nel box grid (in kg).
                                  Deve avere le stesse dimensioni di mr_species_da.
        mw_species_g_per_mol (float): Peso molecolare della specie in grammi/mole (g/mol).
                                      Es: Per HO2 (radical idrossile) = 33.007 g/mol
        mw_dry_air_g_per_mol (float, optional): Peso molecolare dell'aria secca in grammi/mole (g/mol).
                                                Default: 28.97 g/mol (valore standard).

    Returns:
        float: Il burden totale della specie in kilogrammi (kg) per l'intero periodo
               e volume simulato.

    Raises:
        TypeError: Se i tipi di input non sono corretti (non xr.DataArray per i dati).
        ValueError: Se gli attributi 'units' sono mancanti o le dimensioni non corrispondono.
    """
    if not isinstance(mr_species, xr.DataArray) or not isinstance(met_ad, xr.DataArray):
        raise TypeError("Gli input 'mr_species' e 'met_ad' devono essere oggetti xarray.DataArray.")
    if not isinstance(mw_species_g_per_mol, (int, float)) or not isinstance(mw_dry_air_g_per_mol, (int, float)):
        raise TypeError("I pesi molecolari devono essere numeri (int o float).")

    # Controlla e converti le unità del mixing ratio se necessario
    mr_units = mr_species.attrs.get('units', 'mol mol-1').lower().strip()
    conversion_factor_to_mol_mol = 1.0 # Fattore per convertire a mol/mol
    
    if 'ppb' in mr_units or 'ppb' in mr_units:
        conversion_factor_to_mol_mol = 1e9 # 1 ppbv = 1e-9 mol/mol
        #print(f"DEBUG: Unità rilevate per il mixing ratio: '{mr_units}'. Convertendo da ppbv a mol/mol (dividendo per 1e9).")
        print(f"  ")
    elif 'ppt' in mr_units or 'ppt' in mr_units:
        conversion_factor_to_mol_mol = 1e6 # 1 ppmv = 1e-6 mol/mol
        print(f"DEBUG: Unità rilevate per il mixing ratio: '{mr_units}'. Convertendo da ppmv a mol/mol (dividendo per 1e6).")
        print(f"   ")
    elif 'mol mol-1' in mr_units or 'v/v' in mr_units or 'mol/mol' in mr_units:
        #print(f"DEBUG: Unità rilevate per il mixing ratio: '{mr_units}'. Già in mol/mol, nessuna conversione necessaria.")
        print(f"  ")
    else:
        print(f"ATTENZIONE: Unità di mixing ratio '{mr_units}' non riconosciute. Assumo 'mol mol-1'.")

    # Converti i pesi molecolari da g/mol a kg/mol per coerenza con Met_AD (che è in kg)
    mw_species_kg_per_mol = mw_species_g_per_mol / 1000.0
    mw_dry_air_kg_per_mol = mw_dry_air_g_per_mol / 1000.0

    # Calcola la massa della specie in ogni box grid (in kg)
    # xarray gestisce automaticamente l'allineamento delle dimensioni e le operazioni elemento per elemento.
    # Verifica preventiva che le dimensioni siano compatibili
    if mr_species.shape != met_ad.shape:
         # È possibile che le dimensioni siano compatibili ma l'ordine sia diverso,
         # xarray gestirebbe l'allineamento per nome delle dimensioni.
         # Una verifica più rigorosa potrebbe essere: set(mr_species.dims) != set(met_ad.dims)
         # Per GEOS-Chem, solitamente le dimensioni sono le stesse.
        print("ATTENZIONE: Le forme dei DataArray di mixing ratio e Met_AD sono diverse.")
        print(f"Shape MR: {mr_species.shape}, Shape Met_AD: {met_ad.shape}")
        # xarray dovrebbe comunque allinearli per nome delle coordinate se esistono.

    massa_specie_per_box = (mr_species / conversion_factor_to_mol_mol) * \
                           (met_ad / mw_dry_air_kg_per_mol) * \
                           mw_species_kg_per_mol

    # Somma su tutte le dimensioni per ottenere il burden totale
    total_annual_average_burden_kg = (massa_specie_per_box.sum().compute().item())/12
    monthly_burden_kg = (massa_specie_per_box.sum(dim=['lat', 'lon', 'lev'], skipna=True).compute())

    return total_annual_average_burden_kg,  monthly_burden_kg

# --- Esempio di Utilizzo della Funzione Semplificata ---

# Pesi molecolari
mw_H2CO = 30.031 # g/mol
mw_HO2 = 33.0058   # g/mol
mw_dry_air = 28.97 # g/mol (puoi ometterlo se usi il default nella funzione)

print("\n--- SCENARIO 1 ---")
# S1 Chiama la Funzione
burden_H2CO_S1_ann, burden_H2CO_S1_month = calcola_burden(
    mr_species=H2CO_X_S1,
    met_ad=Air_mass_grid_S1,
    mw_species_g_per_mol=mw_H2CO,
    mw_dry_air_g_per_mol=mw_dry_air
    )
print(f"(H2CO): global annual average burden: {burden_H2CO_S1_ann/1e9} Tg")

burden_HO2_S1_ann, burden_HO2_S1_month  = calcola_burden(
    mr_species=HO2_X_S1,
    met_ad=Air_mass_grid_S1,
    mw_species_g_per_mol=mw_HO2,
    mw_dry_air_g_per_mol=mw_dry_air
    )
print(f"(HO2): global annual average burden: {burden_HO2_S1_ann/1e9} Tg")

# S2 Chiama la Funzione
print("\n--- SCENARIO 2 ---")
burden_H2CO_S2_ann, burden_H2CO_S2_month = calcola_burden(
    mr_species=H2CO_X_S2,
    met_ad=Air_mass_grid_S2,
    mw_species_g_per_mol=mw_H2CO,
    mw_dry_air_g_per_mol=mw_dry_air
    )
print(f"(H2CO): global annual average burden: {burden_H2CO_S2_ann/1e9} Tg")

burden_HO2_S2_ann, burden_HO2_S2_month = calcola_burden(
    mr_species=HO2_X_S2,
    met_ad=Air_mass_grid_S2,
    mw_species_g_per_mol=mw_HO2,
    mw_dry_air_g_per_mol=mw_dry_air
    )
print(f"(HO2): global annual average burden: {burden_HO2_S2_ann/1e9} Tg")

# 3. S3 Chiama la Funzione
print("\n--- SCENARIO 3 ---")
burden_H2CO_S3_ann, burden_H2CO_S3_month = calcola_burden(
    mr_species=H2CO_X_S3,
    met_ad=Air_mass_grid_S3,
    mw_species_g_per_mol=mw_H2CO,
    mw_dry_air_g_per_mol=mw_dry_air
    )
print(f"(H2CO): global annual average burden: {burden_H2CO_S3_ann/1e9} Tg")

burden_HO2_S3_ann, burden_HO2_S3_month = calcola_burden(
    mr_species=HO2_X_S3,
    met_ad=Air_mass_grid_S3,
    mw_species_g_per_mol=mw_HO2,
    mw_dry_air_g_per_mol=mw_dry_air
    )
print(f"(HO2): global annual average burden: {burden_HO2_S3_ann/1e9} Tg")

print(f"    ")
print(f"########## (H2CO)")

Absolute_difference_S2_S1_H2CO = burden_H2CO_S2_ann - burden_H2CO_S1_ann
print(f"Absolute difference S2-S1 (H2CO): {Absolute_difference_S2_S1_H2CO/1e9} Tg")

Absolute_difference_S3_S1_H2CO = burden_H2CO_S3_ann - burden_H2CO_S1_ann
print(f"Absolute difference S3-S1 (H2CO): {Absolute_difference_S3_S1_H2CO/1e9} Tg")

Percentage_S2_S1_H2CO = ((burden_H2CO_S2_ann - burden_H2CO_S1_ann)/(burden_H2CO_S1_ann))*100
print(f"% S2 vs S1 (H2CHO): {Percentage_S2_S1_H2CO}")

Percentage_S3_S1_H2CO = ((burden_H2CO_S3_ann - burden_H2CO_S1_ann)/(burden_H2CO_S1_ann))*100
print(f"% S3 vs S1 (H2CHO): {Percentage_S3_S1_H2CO}")

print(f"     ")
print(f"########## (HO2)")

Absolute_difference_S2_S1 = burden_HO2_S2_ann - burden_HO2_S1_ann
print(f"Absolute difference S2-S1 (HO2): {Absolute_difference_S2_S1/1e9} Tg")

Absolute_difference_S3_S1 = burden_HO2_S3_ann - burden_HO2_S1_ann
print(f"Absolute difference S3-S1 (HO2): {Absolute_difference_S3_S1/1e9} Tg")

Percentage_S2_S1 = ((burden_HO2_S2_ann - burden_HO2_S1_ann)/(burden_HO2_S1_ann))*100
print(f"% S2 vs S1 (HO2): {Percentage_S2_S1}")

Percentage_S3_S1 = ((burden_HO2_S3_ann - burden_HO2_S1_ann)/(burden_HO2_S1_ann))*100
print(f"% S3 vs S1 (HO2): {Percentage_S3_S1}")

data = {
    'Metric': [
        'H2CO Burden (S1)',
        'H2CO Burden (S2)',
        'H2CO Burden (S3)',
        'Abs. Diff. S2-S1 (H2CO)',
        'Abs. Diff. S3-S1 (H2CO)',
        'Perc. S2 vs S1 (H2CO)',
        'Perc. S3 vs S1 (H2CO)',
        'HO2 Burden (S1)',
        'HO2 Burden (S2)',
        'HO2 Burden (S3)',
        'Abs. Diff. S2-S1 (HO2)',
        'Abs. Diff. S3-S1 (HO2)',
        'Perc. S2 vs S1 (HO2)',
        'Perc. S3 vs S1 (HO2)'
    ],
    'Value': [
        burden_H2CO_S1_ann / 1e9,
        burden_H2CO_S2_ann / 1e9,
        burden_H2CO_S3_ann / 1e9,
        Absolute_difference_S2_S1_H2CO / 1e9,
        Absolute_difference_S3_S1_H2CO / 1e9,
        Percentage_S2_S1_H2CO,
        Percentage_S3_S1_H2CO,
        burden_HO2_S1_ann / 1e9,
        burden_HO2_S2_ann / 1e9,
        burden_HO2_S3_ann / 1e9,
        Absolute_difference_S2_S1 / 1e9,
        Absolute_difference_S3_S1 / 1e9,
        Percentage_S2_S1,
        Percentage_S3_S1
    ],
    'Unit': [
        'Tg', 'Tg', 'Tg', 'Tg', 'Tg', '%', '%',
        'Tg', 'Tg', 'Tg', 'Tg', 'Tg', '%', '%'
    ]
}

# Crea il DataFrame di Pandas
df_results = pd.DataFrame(data)

# Formatta i valori numerici per una migliore leggibilità
# Ad esempio, a 4 cifre decimali, o notazione scientifica se preferisci
df_results['Value'] = df_results['Value'].map(lambda x: f"{x:.4f}" if abs(x) >= 0.0001 else f"{x:.2e}") # Formatta in base alla grandezza

print("\n--- TABLE: Summary of the results ---")
print(df_results.to_string(index=False)) # .to_string() per una stampa completa senza troncamenti

ho2_time_series_S1 = burden_HO2_S1_month/1e9
ho2_time_series_S2 = burden_HO2_S2_month/1e9
ho2_time_series_S3 = burden_HO2_S3_month/1e9

ho2_time_series_S1.attrs['units'] = 'Tg'
ho2_time_series_S2.attrs['units'] = 'Tg'
ho2_time_series_S3.attrs['units'] = 'Tg'

h2co_time_series_S1 = burden_H2CO_S1_month/1e9
h2co_time_series_S2 = burden_H2CO_S2_month/1e9
h2co_time_series_S3 = burden_H2CO_S3_month/1e9

h2co_time_series_S1.attrs['units'] = 'Tg'
h2co_time_series_S2.attrs['units'] = 'Tg'
h2co_time_series_S3.attrs['units'] = 'Tg'

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# 1. Create the figure and a GridSpec object
# The GridSpec is defined with 2 rows and 2 columns.
# This conceptually divides our plotting area into a 2x2 grid.
fig = plt.figure(figsize=(10, 8)) # Adjust figure size as needed
gs = gridspec.GridSpec(2, 2, figure=fig)

#Obtaining units
try:
    units = H2CO_S1_mean_surf.attrs['units']
except KeyError:
    units = 'Units not known' # Fallback se le unità non sono negli attributi

print(H2CO_S1_mean_surf.attrs['units'])

# 2. Define the subplots using the GridSpec
# Top-left plot: Spans row 0, column 0
ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())

# L'argomento transform=ccrs.PlateCarree() è cruciale e dice a Cartopy che i dati
# (latitudine e longitudine) sono nella proiezione Plate Carree (coordinate geografiche standard).
H2CO_S1_mean_surf.plot(ax=ax1,
                    transform=ccrs.PlateCarree(),
                    cmap='Greys', # Scegli una colormap (es. 'viridis', 'plasma', 'jet', 'RdBu_r')
                    cbar_kwargs={'label': f' H2CO ({units})', 'orientation': 'horizontal'}) # Etichetta per la barra colori

# Aggiungi elementi della mappa
ax1.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax1.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.6)
ax1.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray') # Opzionale, per distinguere terra
ax1.add_feature(cfeature.OCEAN, facecolor='lightblue') # Opzionale, per distinguere mare

# Imposta l'estensione globale
ax1.set_global()

# Aggiungi le griglie e le etichette delle coordinate
ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Aggiungi un titolo
#time_val = H2CO_S1_mean_surf.time.values
title_str = f'Global annual average of H2CO at surface level' # a {time_val}'
#if 'lev' in H2CO_ds_ppb.dims:
 #   # Se il livello è stato selezionato, aggiungilo al titolo
  #  lev_val = H2CO_ds_ppb.lev.values
   # title_str += f' (Surface level: {lev_val})'
    
plt.title(title_str)

ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())

#Obtaining units
try:
    units2 = HO2_S1_mean_surf.attrs['units']
except KeyError:
    units = 'Units not known' # Fallback se le unità non sono negli attributi
# Plot dei dati
# L'argomento transform=ccrs.PlateCarree() è cruciale e dice a Cartopy che i dati
# (latitudine e longitudine) sono nella proiezione Plate Carree (coordinate geografiche standard).
HO2_S1_mean_surf.plot(ax=ax2,
                    transform=ccrs.PlateCarree(),
                    cmap='cool', # Scegli una colormap (es. 'viridis', 'plasma', 'jet', 'RdBu_r')
                    cbar_kwargs={'label': f' HO2  ({units2})','orientation': 'horizontal'}) # Etichetta per la barra colori

# Aggiungi elementi della mappa
ax2.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax2.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.6)
ax2.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray') # Opzionale, per distinguere terra
ax2.add_feature(cfeature.OCEAN, facecolor='lightblue') # Opzionale, per distinguere mare

# Imposta l'estensione globale
ax2.set_global()

# Aggiungi le griglie e le etichette delle coordinate
ax2.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Aggiungi un titolo
title_str_HO2 = f' Global annual average HO2 at surface level'
#if 'lev' in HO2_S1_mean_surf.dims:
    # Se il livello è stato selezionato, aggiungilo al titolo
 #   lev_val = HO2_S1_mean_surf.lev.values
  #  title_str += f' (Surface level: {lev_val})'
# plt.title(title_str)
ax2.set_title(title_str_HO2)

# Create the time series plot
ax3 = fig.add_subplot(gs[1, :])

# Define colors for clarity - NEW: one color per scenario
scenario_colors = {
    'S1': 'orange',
    'S2': 'black',
    'S3': 'red'
}

# Plot HO2 on the first axis (left Y-axis)
# Removed markers, assigned colors by scenario
ho2_time_series_S1.plot(ax=ax3, marker='', linestyle='-', color=scenario_colors['S1'], label='HO2_S1')
ho2_time_series_S2.plot(ax=ax3, marker='', linestyle='-', color=scenario_colors['S2'], label='HO2_S2')
ho2_time_series_S3.plot(ax=ax3, marker='', linestyle='-', color=scenario_colors['S3'], label='HO2_S3')

# Set labels and title for the first axis
ax3.set_xlabel('Time')
ax3.set_ylabel(f'HO2 Burden ({ho2_time_series_S1.attrs["units"]})') # Changed color here too #, color=scenario_colors['S1']
ax3.tick_params(axis='y') # And here # , labelcolor=scenario_colors['S1']
ax3.set_title('HO2 and H2CO burden over Time')
ax3.grid(True)

# Create a second Y-axis that shares the same X-axis
ax4 = ax3.twinx()

# Plot H2CO on the second axis (right Y-axis)
# Removed markers, assigned colors by scenario
h2co_time_series_S1.plot(ax=ax4, marker='', linestyle='--', color=scenario_colors['S1'], label='H2CO_S1')
h2co_time_series_S2.plot(ax=ax4, marker='', linestyle='--', color=scenario_colors['S2'], label='H2CO_S2')
h2co_time_series_S3.plot(ax=ax4, marker='', linestyle='--', color=scenario_colors['S3'], label='H2CO_S3')

ax4.set_title('')

# Set labels for the second axis
ax4.set_ylabel(f'H2CO Burden ({h2co_time_series_S1.attrs["units"]})') 
ax4.tick_params(axis='y') # And here

# --- 4. Create a Combined Legend ---
lines, labels = ax3.get_legend_handles_labels()
lines2, labels2 = ax4.get_legend_handles_labels()
ax4.legend(lines + lines2, labels + labels2, loc='lower left', bbox_to_anchor=(0, 0))

plt.tight_layout()
plt.show()


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs # Assicurati di importare cartopy
import cartopy.feature as cfeature # Assicurati di importare cartopy
import numpy as np # Assicurati di importare numpy
import matplotlib.colors as colors
#%matplotlib qt


# --- Define the desired min, max, and center for the colormap for the ABSOLUTE DIFFERENCE ---
vmin_diff = -0.3
vmax_diff = 0.3
vcenter_diff = 0.0
num_levels_diff = 21 # Number of discrete levels

# --- Define the desired min, max, and center for the colormap for the PERCENTAGE DIFFERENCE ---
vmin_perc = -5.0 # Example: -50%
vmax_perc = 5.0  # Example: +50%
vcenter_perc = 0.0
num_levels_perc = 21

# Colormaps
cmap_diff = 'coolwarm' # Diverging colormap for absolute difference (or RdBu_r)
cmap_perc = 'coolwarm' # Another good diverging colormap for percentage difference

# Create the figure and a GridSpec object
# The GridSpec is defined with 5 rows and 2 columns.
# This conceptually divides our plotting area into a 5x2 grid.
fig = plt.figure(figsize=(14, 20)) # Adjust figure size as needed
gs = gridspec.GridSpec(5, 2, figure=fig)

# [0:0] Define the subplots using the GridSpec - SURFACE LEVEL

#Obtaining units
try:
    units = H2CO_S1_mean_surf.attrs['units']
except KeyError:
    units = 'Units not known' # Fallback se le unità non sono negli attributi

ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())

# L'argomento transform=ccrs.PlateCarree() è cruciale e dice a Cartopy che i dati
# (latitudine e longitudine) sono nella proiezione Plate Carree (coordinate geografiche standard).
H2CO_S1_mean_surf.plot(ax=ax1,
                    transform=ccrs.PlateCarree(),
                    cmap='Greys', # Scegli una colormap (es. 'viridis', 'plasma', 'jet', 'RdBu_r')
                    cbar_kwargs={'label': f' H2CO ({units})', 'pad': 0.05}) #, 'orientation': 'horizontal'}) # Etichetta per la barra colori

# Aggiungi elementi della mappa
ax1.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax1.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.6)
ax1.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray') # Opzionale, per distinguere terra
ax1.add_feature(cfeature.OCEAN, facecolor='lightblue') # Opzionale, per distinguere mare

# Imposta l'estensione globale
ax1.set_global()

# Aggiungi le griglie e le etichette delle coordinate
ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Aggiungi un titolo
#time_val = H2CO_S1_mean_surf.time.values
title_str = f'S1 - H2CO - Annual average - surface level' # a {time_val}'
#if 'lev' in H2CO_ds_ppb.dims:
 #   # Se il livello è stato selezionato, aggiungilo al titolo
  #  lev_val = H2CO_ds_ppb.lev.values
   # title_str += f' (Surface level: {lev_val})'
    
plt.title(title_str)

# [0:1] Define the subplots using the GridSpec - MIDDLE LEVEL

#Obtaining units
try:
    units = H2CO_S1_mean_middle.attrs['units']
except KeyError:
    units = 'Units not known' # Fallback se le unità non sono negli attributi

ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())

# L'argomento transform=ccrs.PlateCarree() è cruciale e dice a Cartopy che i dati
# (latitudine e longitudine) sono nella proiezione Plate Carree (coordinate geografiche standard).
H2CO_S1_mean_middle.plot(ax=ax2,
                    transform=ccrs.PlateCarree(),
                    cmap='Greys', # Scegli una colormap (es. 'viridis', 'plasma', 'jet', 'RdBu_r')
                    cbar_kwargs={'label': f' H2CO ({units})', 'pad': 0.05}) #, 'orientation': 'horizontal'}) # Etichetta per la barra colori

# Aggiungi elementi della mappa
ax2.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax2.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.6)
ax2.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray') # Opzionale, per distinguere terra
ax2.add_feature(cfeature.OCEAN, facecolor='lightblue') # Opzionale, per distinguere mare

# Imposta l'estensione globale
ax2.set_global()

# Aggiungi le griglie e le etichette delle coordinate
ax2.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Aggiungi un titolo
#time_val = H2CO_S1_mean_middle.time.values
title_str2 = f'S1 - H2CO - Annual average - Middle troposphere' # a {time_val}'
#if 'lev' in H2CO_ds_ppb.dims:
 #   # Se il livello è stato selezionato, aggiungilo al titolo
  #  lev_val = H2CO_ds_ppb.lev.values
   # title_str += f' (Surface level: {lev_val})'
    
plt.title(title_str2)


##############################################################################
###                                  S2-S1                             #######
##############################################################################

# [1:0] Define the subplots using the GridSpec - Abs diff at SURFACE LEVEL
# S1-S1

ax3 = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 1: Absolute Difference S2-S1 (Left Plot) ---

H2CO_S2_S1_diff_surf.plot(ax=ax3,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_diff,
                       levels=np.linspace(vmin_diff, vmax_diff, num_levels_diff),
                       norm=colors.TwoSlopeNorm(vmin=vmin_diff, vcenter=vcenter_diff, vmax=vmax_diff),
                       cbar_kwargs={'label': f'Abs. Diff. ({units})'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax3.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax3.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax3.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax3.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax3.set_global()
ax3.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_diff1 = f'S2-S1 Abs. Diff. at surface level'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax3.set_title(title_str_diff1)


# [1:1] Define the subplots using the GridSpec - Abs diff at MIDDLE LEVEL
# S2-S1

ax4 = fig.add_subplot(gs[1, 1], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 2: Absolute Difference S2-S1 (Right Plot) MIDDLE TROPO ---

H2CO_S2_S1_diff_middle.plot(ax=ax4,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_diff,
                       levels=np.linspace(vmin_diff, vmax_diff, num_levels_diff),
                       norm=colors.TwoSlopeNorm(vmin=vmin_diff, vcenter=vcenter_diff, vmax=vmax_diff),
                       cbar_kwargs={'label': f'Abs. Diff. ({units})'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax4.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax4.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax4.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax4.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax4.set_global()
ax4.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_diff2 = f'S2-S1 Abs. Diff. at middle troposphere'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax4.set_title(title_str_diff2)


######################################################
######################################################

# [2:0] Define the subplots using the GridSpec - % at SURFACE LEVEL
# S2-S1 %

ax5 = fig.add_subplot(gs[2, 0], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 1: % S2-S1 (Left Plot) ---

H2CO_S2_S1_perc_surf.plot(ax=ax5,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_perc,
                       levels=np.linspace(vmin_perc, vmax_perc, num_levels_perc),
                       norm=colors.TwoSlopeNorm(vmin=vmin_perc, vcenter=vcenter_perc, vmax=vmax_perc),
                       cbar_kwargs={'label': f'%'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax5.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax5.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax5.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax5.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax5.set_global()
ax5.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_perc1 = f'S2-S1 % Difference at surface level'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax5.set_title(title_str_perc1)

# [2:1] Define the subplots using the GridSpec - % at MIDDLE LEVEL
# S2-S1 %

ax6 = fig.add_subplot(gs[2, 1], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 2: % Difference S2-S1 (Right Plot) MIDDLE TROPO ---

H2CO_S2_S1_perc_middle.plot(ax=ax6,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_perc,
                       levels=np.linspace(vmin_perc, vmax_perc, num_levels_perc),
                       norm=colors.TwoSlopeNorm(vmin=vmin_perc, vcenter=vcenter_perc, vmax=vmax_perc),
                       cbar_kwargs={'label': f'%'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax6.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax6.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax6.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax6.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax6.set_global()
ax6.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_perc2 = f'S2-S1 % Difference at middle troposphere'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax6.set_title(title_str_perc2)


##############################################################################
###                                  S3-S1                             #######
##############################################################################

# [3:0] Define the subplots using the GridSpec - Abs diff at SURFACE LEVEL
# S3-S1

ax7 = fig.add_subplot(gs[3, 0], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 1: Absolute Difference S3-S1 (Left Plot) ---

H2CO_S3_S1_diff_surf.plot(ax=ax7,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_diff,
                       levels=np.linspace(vmin_diff, vmax_diff, num_levels_diff),
                       norm=colors.TwoSlopeNorm(vmin=vmin_diff, vcenter=vcenter_diff, vmax=vmax_diff),
                       cbar_kwargs={'label': f'Abs. Diff. ({units})'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax7.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax7.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax7.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax7.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax7.set_global()
ax7.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_diff5 = f'S3-S1 Abs. Diff. at surface level'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax7.set_title(title_str_diff5)

# [3:1] Define the subplots using the GridSpec - Abs diff at MIDDLE LEVEL
# S3-S1

ax8 = fig.add_subplot(gs[3, 1], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 2: Absolute Difference S2-S1 (Right Plot) MIDDLE TROPO ---

H2CO_S3_S1_diff_middle.plot(ax=ax8,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_diff,
                       levels=np.linspace(vmin_diff, vmax_diff, num_levels_diff),
                       norm=colors.TwoSlopeNorm(vmin=vmin_diff, vcenter=vcenter_diff, vmax=vmax_diff),
                       cbar_kwargs={'label': f'Abs. Diff. ({units})'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax8.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax8.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax8.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax8.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax8.set_global()
ax8.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_diff6 = f'S3-S1 Abs. Diff. at middle troposphere'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax8.set_title(title_str_diff6)

######################################################
######################################################

# [4:0] Define the subplots using the GridSpec - % at SURFACE LEVEL
# S3-S1 %

ax9 = fig.add_subplot(gs[4, 0], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 1: % S3-S1 (Left Plot) ---

H2CO_S3_S1_perc_surf.plot(ax=ax9,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_perc,
                       levels=np.linspace(vmin_perc, vmax_perc, num_levels_perc),
                       norm=colors.TwoSlopeNorm(vmin=vmin_perc, vcenter=vcenter_perc, vmax=vmax_perc),
                       cbar_kwargs={'label': f'%'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax9.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax9.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax9.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax9.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax9.set_global()
ax9.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_perc5 = f'S3-S1 % Difference at surface level'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax9.set_title(title_str_perc5)

# [4:1] Define the subplots using the GridSpec - % at MIDDLE LEVEL
# S3-S1 %

ax10 = fig.add_subplot(gs[4, 1], projection=ccrs.PlateCarree()) # gs[row, column]

# --- Plot 2: % Difference S3-S1 (Right Plot) MIDDLE TROPO ---

H2CO_S3_S1_perc_middle.plot(ax=ax10,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_perc,
                       levels=np.linspace(vmin_perc, vmax_perc, num_levels_perc),
                       norm=colors.TwoSlopeNorm(vmin=vmin_perc, vcenter=vcenter_perc, vmax=vmax_perc),
                       cbar_kwargs={'label': f'%'}) #, 'orientation': 'horizontal'})

# Add map elements for the first plot
ax10.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax10.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax10.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax10.add_feature(cfeature.OCEAN, facecolor='lightblue')
ax10.set_global()
ax10.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)

# Add title for the first plot
title_str_perc7 = f'S3-S1 % Difference at middle troposphere'
#if 'lev' in H2CO_X_S2_S1_diff.dims:
 #   lev_val_diff = H2CO_X_S2_S1_diff.lev.values
  #  title_str_diff += f' (Surface level: {lev_val_diff})'
ax10.set_title(title_str_perc7)

plt.tight_layout() # Adatta i margini per evitare sovrapposizioni
