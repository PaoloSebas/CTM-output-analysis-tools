# IMPORT necessary Libraries
# xarray for handling labeled multi-dimensional arrays (NetCDF files)
import xarray as xr
# matplotlib for plotting functionalities
import matplotlib.pyplot as plt
# numpy for numerical operations
import numpy as np

# gcpy is a set of utilities for GEOS-Chem output analysis and plotting.
# We explicitly import the necessary plotting functions.
from gcpy.plot.single_panel import single_panel

# --- Unused imports removed for cleaner code ---
# netCDF4 is often used by xarray internally, no need to import explicitly unless direct use.
# from gcpy.plot.compare_single_level import compare_single_level # Not used in this script
# from gcpy.plot.compare_zonal_mean import compare_zonal_mean   # Not used in this script
# from matplotlib.backends.backend_pdf import PdfPages           # Not used in this script

# --- Configuration Section ---
# Define the species to analyze
SPECIES_NAME = 'HO2'
# Define the variable name in the NetCDF file for this species
# GEOS-Chem typically names these 'SpeciesConcVV_SPECIESNAME'
SPECIES_VAR_NAME = f'SpeciesConcVV_{SPECIES_NAME}'

# Define input file paths for initial and final scenarios
# Using f-strings for clarity in paths
# NOTE: The original script used 2020-01-01 files but plot title mentioned 08/2019.
# Assuming the files provided are correct for the analysis period.
# If these files represent a specific day *within* a longer August 2019 period,
# the title should reflect that carefully.
INITIAL_FILE_PATH = '/mnt/d/PPO/NEW_SET_2025/SCENARIO_1/SCENARIO1_GEOSChem.SpeciesConc.20200101_0000z.nc4'
FINAL_FILE_PATH = '/mnt/d/PPO/NEW_SET_2025/SCENARIO_3/SCENARIO3_GEOSChem.SpeciesConc.20200101_0000z.nc4'

# Define plot properties
PLOT_TITLE_DATE = '01/2020' # Adjusted to match filenames, change if 08/2019 is truly intended
PLOT_LEVEL = 0 # Surface level for the plot
COLORMAP = plt.get_cmap("seismic") # Good for diverging data (positive/negative changes)
COLOR_SCALE_MIN = -6 # Hardcoded vmin for plot color scale
COLOR_SCALE_MAX = 6  # Hardcoded vmax for plot color scale
GEOGRAPHICAL_EXTENT = [-180, 180, -90, 90] # Global extent


# --- READING COMPLETE DATASETS ---

# Initial Scenario Dataset
print(f"Reading initial dataset from: {INITIAL_FILE_PATH}")
ds_initial = xr.open_dataset(INITIAL_FILE_PATH)
print('\n--- Initial Dataset Info ---')
print(ds_initial)
print('------------------------------------------------------------------------------------')

# Final Scenario Dataset
print(f"\nReading final dataset from: {FINAL_FILE_PATH}")
ds_final = xr.open_dataset(FINAL_FILE_PATH)
print('\n--- Final Dataset Info ---')
print(ds_final)
print('------------------------------------------------------------------------------------')


# --- SUBSETTING on the defined species ---
# Extract the DataArray for the specified species from both datasets
print(f'\n########################## SUBSETTING ON SPECIES: {SPECIES_NAME} ##################################')
species_data_initial = ds_initial[SPECIES_VAR_NAME]
species_data_final = ds_final[SPECIES_VAR_NAME]

print('\n--- Initial Dataset after SUBSETTING ---')
print(species_data_initial)
print('------------------------------------------------------------------------------------')

print('\n--- Final Dataset after SUBSETTING ---')
print(species_data_final)
print('------------------------------------------------------------------------------------')


# --- CALCULATIONS ---

# Convert concentration from mol mol-1 (dry air) to ppb (parts per billion)
# GEOS-Chem SpeciesConcVV is typically mol mol-1.
# This conversion is crucial for interpreting values as ppb.
CONVERSION_FACTOR_TO_PPB = 1e+9 # 1e9 mol/mol = 1 ppb

species_data_initial_ppb = species_data_initial * CONVERSION_FACTOR_TO_PPB
species_data_final_ppb = species_data_final * CONVERSION_FACTOR_TO_PPB

# Summing over time (if multiple time steps exist in the file)
# This will average/sum the concentrations across time for each lat/lon/lev grid cell.
# If only one time step, it effectively keeps the original data.
species_sum_time_initial_ppb = species_data_initial_ppb.sum(dim='time', keep_attrs=True)
species_sum_time_final_ppb = species_data_final_ppb.sum(dim='time', keep_attrs=True)

# --- Burden Calculation (Important Note) ---
# The original "burden" calculation was a sum of mixing ratios.
# A true atmospheric burden (in moles or mass) requires multiplying by the air mass of each grid cell.
# For simplicity, here we retain the original logic of summing mixing ratios across all dimensions,
# but it's important to be aware this is NOT a true burden in terms of total moles or mass.
# It's an index proportional to burden for comparison.

initial_total_mixing_ratio_sum = np.sum(species_data_initial_ppb.values)
final_total_mixing_ratio_sum = np.sum(species_data_final_ppb.values)

# Percentage Variation of Total Mixing Ratio Sum
burden_variation_percentage = (
    (final_total_mixing_ratio_sum - initial_total_mixing_ratio_sum) / initial_total_mixing_ratio_sum
) * 100

# Burden-like calculation for the surface level (lev=0)
initial_surface_mixing_ratio_sum = np.sum(species_data_initial_ppb.isel(lev=PLOT_LEVEL).values)
final_surface_mixing_ratio_sum = np.sum(species_data_final_ppb.isel(lev=PLOT_LEVEL).values)

# Percentage Variation of Surface Mixing Ratio Sum
burden_variation_surface_percentage = (
    (final_surface_mixing_ratio_sum - initial_surface_mixing_ratio_sum) / initial_surface_mixing_ratio_sum
) * 100

print(f'\n--- Burden-like Calculations for {SPECIES_NAME} ---')
print(f'Sum of Initial mixing ratios (ppb): {initial_total_mixing_ratio_sum:.2f}')
print(f'Sum of Final mixing ratios (ppb): {final_total_mixing_ratio_sum:.2f}')
print(f'Total Mixing Ratio Sum % Variation: {burden_variation_percentage:.2f}%')
print(f'Surface Mixing Ratio Sum % Variation: {burden_variation_surface_percentage:.2f}%')
print('------------------------------------------------------------------------------------')


# --- Calculate Differences for Plotting ---

# Percentage change at each grid point (after summing over time)
# Note: Handle division by zero or very small initial values if they are possible.
# Adding .where() to avoid NaN from 0/0 and inf from X/0 where X is non-zero
diff_plot_percentage = (
    ((species_sum_time_final_ppb - species_sum_time_initial_ppb) / species_sum_time_initial_ppb) * 100
).where(species_sum_time_initial_ppb != 0) # Set to NaN where initial is zero

# Absolute difference in ppb at each grid point (after summing over time)
diff_plot_absolute = species_sum_time_final_ppb - species_sum_time_initial_ppb

# Find min/max values for potential dynamic color scale (for percentage plot)
max_val_percentage = np.max(diff_plot_percentage.values[~np.isnan(diff_plot_percentage.values)]) # Ignore NaNs
min_val_percentage = np.min(diff_plot_percentage.values[~np.isnan(diff_plot_percentage.values)]) # Ignore NaNs

# Helper function to find the value with the largest absolute magnitude
def get_max_modulus(val1, val2):
    """Returns the value (either val1 or val2) that has the largest absolute magnitude."""
    return val1 if abs(val1) > abs(val2) else val2

# Calculate the maximum absolute range for the percentage plot, if dynamic limits were used
# (Currently, hardcoded limits COLOR_SCALE_MIN/MAX are used below)
dynamic_abs_max_percentage = abs(get_max_modulus(max_val_percentage, min_val_percentage))

print(f'\n--- Difference Statistics for {SPECIES_NAME} (after summing over time) ---')
print(f'Max Percentage Change: {max_val_percentage:.2f}%')
print(f'Min Percentage Change: {min_val_percentage:.2f}%')
print(f'Max Absolute Difference (ppb): {np.max(diff_plot_absolute.values):.2e}') # Use scientific notation
print(f'Min Absolute Difference (ppb): {np.min(diff_plot_absolute.values):.2e}') # Use scientific notation
print(f'Calculated dynamic max scale for percentage plot (abs): {dynamic_abs_max_percentage:.2f}')
print('------------------------------------------------------------------------------------')


# --- PLOTTING ---
print(f'\nGenerating plot for {SPECIES_NAME} % Change at Level {PLOT_LEVEL}...')

single_panel(
    # Data to plot: Percentage change at the specified level
    diff_plot_percentage.isel(lev=PLOT_LEVEL),
    # Title for the plot, using constants for consistency
    title=f'{SPECIES_NAME} % CHANGE SCENARIO3 vs SCENARIO1 {PLOT_TITLE_DATE} - LEV={PLOT_LEVEL}',
    # Colormap for diverging data
    comap=COLORMAP,
    # Whether to use a logarithmic color scale (False for percentage change)
    log_color_scale=False,
    # Define color bar limits (vmin/vmax).
    # Using hardcoded values for consistent comparison as per original script.
    # To use dynamically calculated values, uncomment the lines below and comment out the hardcoded ones.
    vmin=COLOR_SCALE_MIN,
    vmax=COLOR_SCALE_MAX,
    # vmin= -dynamic_abs_max_percentage, # Alternative: dynamically set limits based on data range
    # vmax= dynamic_abs_max_percentage,   # Alternative: dynamically set limits based on data range
    # Geographical extent of the plot
    extent=GEOGRAPHICAL_EXTENT
)

# Display the plot
plt.show()

print('\nScript execution complete.')
