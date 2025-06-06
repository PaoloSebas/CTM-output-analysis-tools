# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 22:08:27 2025

@author: paolos
"""

# --- Part 1: Explanation of the Logic and Formulas ---
#
# This script calculates the monthly global burden (total mass) of HO2 at two
# different vertical levels (0 and 15) based on GEOS-Chem model output.
# It processes data from three different scenarios, reading monthly species
# concentration and state met files.  It also calculates and plots the absolute
# difference in HO2 concentration (ppb) at Level 0 between Scenario 2 and 1,
# and Scenario 3 and 1.
#
# The core logic involves the following steps for each scenario, each month,
# and each specified vertical level:
#
# 1. Read Data:
#   - Open the GEOS-Chem SpeciesConc file to get the volume mixing ratio of HO2
#     ('SpeciesConcVV_HO2') at the specified vertical level (lev) and the first
#     time step.
#   - Open the corresponding GEOS-Chem StateMet file to get the air density
#     ('Met_AIRDEN') at the same vertical level and the grid cell area ('AREA')
#     (Note: AREA is 2D and does not depend on the vertical level).
#
# 2. Calculate Mass Concentration:
#   - The volume mixing ratio of HO2 is used as a proxy for its mole fraction.
#   - The number of moles of HO2 per cubic meter is estimated using the air density
#     at the specific level and the molar mass of air (M_AIR), assuming ideal
#     gas behavior locally.
#     Formula: moles_HO2/m^3 ≈ (mole_fraction_HO2 * air_density_level) / M_AIR
#   - The mass concentration of HO2 (kg/m³) is then calculated by multiplying the
#     moles per cubic meter by the molar mass of HO2 (M_HO2).
#     Formula: mass_ho2/m^3 = (moles_HO2/m^3) * M_HO2
#
# 3. Calculate Mass per Grid Cell:
#   - The mass of HO2 in each grid cell at the specific vertical level is
#     calculated by multiplying the mass concentration at that level by the
#     grid cell area.
#     Formula: mass_ho2_per_cell_level (kg) = mass_ho2/m^3_level * AREA (m^2)
#
# 4. Calculate Global Burden per Level:
#   - The global burden of HO2 for each month at the specific vertical level
#     is obtained by summing the mass of HO2 in all grid cells at that level.
#     Formula: Global Burden_level (kg) = Σ (mass_HO2_per_cell_level) over all lat/lon
#
# 5. Plotting (Burden):
#   - Finally, the script plots the monthly global burden of HO2 for each scenario
#     and each of the two specified vertical levels over the specified time period,
#     with colors and legend grouped by vertical level.
#
# 6. Plotting (Difference Maps):
#   -  Calculate the absolute difference in HO2 concentration (ppb) at Level 0
#       between Scenario 2 and Scenario 1, and Scenario 3 and Scenario 1 for each month.
#   -  Create global maps of these differences for each month.
#
# --- Part 2: Import Libraries ---
import xarray as xr
import numpy as np
import os
print("Current Working Directory:", os.getcwd())
from datetime import date
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import traceback  # Import traceback for detailed error logging

# --- Part 3: Configuration ---
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

# Define the vertical levels to consider
vertical_levels = [0, 15]
level_names = {0: 'Level 0 (Surface)', 15: 'Level 15'}
scenario_names = {'Scenario 1', 'Scenario 2', 'Scenario 3'}

# Define color scales for each level
level_colors = {
    0: [(1.0, 0.5, 0.0), (0.8, 0.4, 0.0), (0.6, 0.3, 0.0)],  # Orange scale for level 0 (RGB tuples)
    15: [(0.8, 0, 0), (0.6, 0, 0), (0.4, 0, 0)]                        # Red scale (RGB tuples)
}
scenario_order = {'Scenario 1': 0, 'Scenario 2': 1, 'Scenario 3': 2}

# Create dictionaries to easily access directories and filenames for each scenario
scenario_directories = {
    'Scenario 1': base_directory_s1,
    'Scenario 2': base_directory_s2,
    'Scenario 3': base_directory_s3,
}

species_conc_base_filenames = {
    'Scenario 1': species_conc_base_filename_s1,
    'Scenario 2': species_conc_base_filename_s2,
    'Scenario 3': species_conc_base_filename_s3,
}

state_met_base_filenames = {
    'Scenario 1': state_met_base_filename_s1,
    'Scenario 2': state_met_base_filename_s2,
    'Scenario 3': state_met_base_filename_s3,
}

# --- Part 4: Constants and Helper Function ---
# Define physical constants used in the calculations
AVOGADRO = 6.02214076e23  # molecules/mol
M_AIR = 28.97e-3          # kg/mol (molar mass of air)
M_HO2 = 33.007e-3         # kg/mol (molar mass of HO2)
PPM_TO_PPB = 1e3

# --- Helper Function to Generate Monthly Dates ---
def generate_monthly_dates(start_str, end_str):
    # Convert start and end date strings to datetime.date objects
    start = date(int(start_str[:4]), int(start_str[4:6]), int(start_str[6:]))
    end = date(int(end_date_str[:4]), int(end_date_str[4:6]), int(end_date_str[6:]))
    dates = []
    current = start
    # Loop through months from start to end date
    while current <= end:
        # Format the current date as YYYYMMDD string and append to the list
        dates.append(current.strftime('%Y%m%d'))
        # Increment to the next month
        if current.month == 12:
            current = date(current.year + 1, 1, 1)
        else:
            current = date(current.year, current.month + 1, 1)
    return dates

# --- Part 5: Main Data Processing and Storage ---
# --- Main Script for Plotting Monthly Global Burden of HO2 at Multiple Levels (kg) ---
print("\n--- Standalone Script: Plotting Monthly Global Burden of HO2 at Multiple Levels (kg) ---")

# Generate the list of monthly dates based on the start and end dates
monthly_dates = generate_monthly_dates(start_date_str, end_date_str)
# Initialize a dictionary to store the monthly global burden of HO2 for each scenario and level
ho2_burden_monthly_kg = {scenario: {level: [] for level in vertical_levels} for scenario in scenario_directories}
# Initialize a dictionary to store monthly HO2 mixing ratios at level 0 for each scenario
ho2_mf_monthly_ppb_level0 = {scenario: [] for scenario in scenario_directories}
# Initialize latitude and longitude (assuming they are consistent across files)
latitude = None
longitude = None
error_occurred = False  # Flag to track if any error occurred during processing

# --- Loop through each scenario ---
for scenario_name, base_dir in scenario_directories.items():
    # Get the base filenames for species concentration and state met for the current scenario
    species_conc_base = species_conc_base_filenames[scenario_name]
    state_met_base = state_met_base_filenames[scenario_name]

    # --- Loop through each monthly date ---
    for date_str in monthly_dates:
        try:
            # Construct the full file paths for the species concentration and state met files
            species_conc_filepath = os.path.join(base_dir, f'{species_conc_base}{date_str}_0000z.nc4')
            state_met_filepath = os.path.join(base_dir, f'{state_met_base}{date_str}_0000z.nc4')

            # --- Check if files exist before attempting to open them ---
            if not os.path.exists(species_conc_filepath):
                print(f"Error: Species concentration file not found: {species_conc_filepath}")
                for level in vertical_levels:
                    ho2_burden_monthly_kg[scenario_name][level].append(np.nan)
                if 0 in vertical_levels:
                    ho2_mf_monthly_ppb_level0[scenario_name].append(None)
                error_occurred = True
                continue  # Skip to the next iteration (next month)
            if not os.path.exists(state_met_filepath):
                print(f"Error: State met file not found: {state_met_filepath}")
                for level in vertical_levels:
                    ho2_burden_monthly_kg[scenario_name][level].append(np.nan)
                if 0 in vertical_levels:
                    ho2_mf_monthly_ppb_level0[scenario_name].append(None)
                error_occurred = True
                continue  # Skip to the next iteration (next month)
            # --- Open the datasets using xarray ---
            try:
                with xr.open_dataset(species_conc_filepath) as ds_species, \
                        xr.open_dataset(state_met_filepath) as ds_state_met:

                    # --- Loop through each vertical level ---
                    for level in vertical_levels:
                        # --- Extract relevant variables at the current level and first time step ---
                        try:
                            ho2_mf = ds_species['SpeciesConcVV_HO2'].isel(time=0, lev=level)  # HO2 volume mixing ratio at current level
                            air_density = ds_state_met['Met_AIRDEN'].isel(time=0, lev=level)  # Air density at current level
                            area = ds_state_met['AREA']  # Grid cell area (does not depend on level)
                        except KeyError as e:
                            print(f"Error: Variable not found in dataset for {scenario_name}, {date_str}, level {level}: {e}")
                            traceback.print_exc()
                            ho2_burden_monthly_kg[scenario_name][level].append(np.nan)
                            if level == 0:
                                ho2_mf_monthly_ppb_level0[scenario_name].append(None)
                            error_occurred = True
                            continue  # Skip to the next level
                        except Exception as e:
                            print(f"Error accessing data for {scenario_name}, {date_str}, level {level}: {e}")
                            traceback.print_exc()
                            ho2_burden_monthly_kg[scenario_name][level].append(np.nan)
                            if level == 0:
                                ho2_mf_monthly_ppb_level0[scenario_name].append(None)
                            error_occurred = True
                            continue  # Skip to the next level

                        # --- Calculate mass concentration of HO2 (kg/m3) ---
                        ho2_mol_fraction = ho2_mf
                        # Assuming partial pressure of HO2 is approximately equal to its mole fraction (simplified for trace gas)
                        partial_pressure_ho2 = ho2_mol_fraction
                        # Calculate moles of HO2 per cubic meter using the ideal gas law approximation (related to air density)
                        moles_ho2_per_m3 = (partial_pressure_ho2 * air_density) / M_AIR
                        # Convert moles per cubic meter to mass per cubic meter using the molar mass of HO2
                        mass_ho2_per_m3 = moles_ho2_per_m3 * M_HO2

                        # --- Calculate the mass of HO2 in each grid cell at the current level (kg) ---
                        mass_ho2_per_cell = mass_ho2_per_m3 * area

                        # --- Calculate the global burden of HO2 at the current level (kg) by summing over all grid cells ---
                        global_burden_ho2_kg = mass_ho2_per_cell.sum(dim=['lat', 'lon']).item()
                        # Append the calculated global burden to the dictionary for the current scenario and level
                        ho2_burden_monthly_kg[scenario_name][level].append(global_burden_ho2_kg)

                        # --- Store HO2 mixing ratio (ppb) at level 0 ---
                        if level == 0:
                            ho2_mf_ppb_level0 = ho2_mf * PPM_TO_PPB * 1e3  # Multiply by 1e9 to get ppb
                            ho2_mf_monthly_ppb_level0[scenario_name].append(ho2_mf_ppb_level0.values)
                            if latitude is None:
                                latitude = ds_species['lat'].values
                                longitude = ds_species['lon'].values
            except Exception as e:
                print(f"Error opening/processing dataset for {scenario_name}, {date_str}: {e}")
                traceback.print_exc()  # Print the full traceback
                for level in vertical_levels:
                    ho2_burden_monthly_kg[scenario_name][level].append(np.nan)
                if 0 in vertical_levels:
                    ho2_mf_monthly_ppb_level0[scenario_name].append(None)
                error_occurred = True

        # --- Handle potential errors ---
        except Exception as e:
            print(f"An unexpected error occurred while processing {scenario_name}, {date_str}: {e}")
            traceback.print_exc()
            for level in vertical_levels:
                ho2_burden_monthly_kg[scenario_name][level].append(np.nan)
            if 0 in vertical_levels:
                ho2_mf_monthly_ppb_level0[scenario_name].append(None)
            error_occurred = True

# --- Part 6: Plotting Global Burden ---
# --- Create the plot of global burden ---
if not error_occurred:  # Only plot if no errors occurred during data processing
    plt.figure(figsize=(15, 8))
    # Convert monthly date strings to datetime objects for proper plotting
    date_objects = [date(int(d[:4]), int(d[4:6]), 1) for d in monthly_dates]

    # --- Plot the results for each scenario and each level with specified colors ---
    for level in vertical_levels:
        for scenario_name in scenario_names:
            data_to_plot = ho2_burden_monthly_kg[scenario_name][level]
            color_index = scenario_order[scenario_name]
            color = level_colors[level][color_index]
            label = f'{scenario_name}, {level_names[level]}'
            print(f"Plotting Burden: Level={level}, Scenario={scenario_name}, Color={color}")  # Debug print
            plt.plot(date_objects, data_to_plot, marker='o', linestyle='-', color=color, label=label)

    # --- Add plot labels and title for burden ---
    plt.xlabel('Month')
    plt.ylabel('Global Burden of HO$_2$ (kg)')
    plt.title('Monthly Global Burden of HO$_2$ at Different Vertical Levels')
    plt.grid(True)

    # --- Create custom legend order for burden ---
    handles, labels = plt.gca().get_legend_handles_labels()
    legend_order = []
    ordered_scenarios = ['Scenario 1', 'Scenario 2', 'Scenario 3']

    # Add Level 0 entries in the desired scenario order
    for scenario_name in ordered_scenarios:
        label_to_find = f'{scenario_name}, {level_names[0]}'
        if label_to_find in labels:
            legend_order.append(labels.index(label_to_find))

    # Add Level 15 entries in the desired scenario order
    for scenario_name in ordered_scenarios:
        label_to_find = f'{scenario_name}, {level_names[15]}'
        if label_to_find in labels:
            legend_order.append(labels.index(label_to_find))

    plt.legend([handles[i] for i in legend_order], [labels[i] for i in legend_order])

    # Format the x-axis ticks to show Month-Year
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y'))
    plt.xticks(rotation=45, ha='right')  # Rotate labels for better readability
    plt.tight_layout()

    # --- Save the burden plot as a PNG file ---
    plt.savefig('ho2_global_burden.png')  # Specify the filename
    plt.close()  # Close the figure to free up memory

else:
    print("Skipping burden plot due to errors during data processing.")


# --- Part 7: Plotting Difference Maps ---
# --- Plotting Global Maps of Absolute HO2 Difference at Level 0 (ppb) ---
print("\n--- Plotting Global Maps of Absolute HO2 Difference at Level 0 (ppb) ---")

if not error_occurred and latitude is not None and longitude is not None:  # only plot if no errors and lat/lon are loaded
    n_months = len(monthly_dates)
    for i in range(n_months):
        if (ho2_mf_monthly_ppb_level0['Scenario 1'][i] is not None and
                ho2_mf_monthly_ppb_level0['Scenario 2'][i] is not None and
                ho2_mf_monthly_ppb_level0['Scenario 3'][i] is not None):

            ho2_diff_s2_s1 = np.abs(ho2_mf_monthly_ppb_level0['Scenario 2'][i] - ho2_mf_monthly_ppb_level0['Scenario 1'][i])
            ho2_diff_s3_s1 = np.abs(ho2_mf_monthly_ppb_level0['Scenario 3'][i] - ho2_mf_monthly_ppb_level0['Scenario 1'][i])

            fig, axes = plt.subplots(1, 2, figsize=(15, 6), subplot_kw={'projection': ccrs.PlateCarree()})
            fig.suptitle(f'Absolute HO$_2$ Difference at Level 0 - {date_objects[i].strftime("%b-%Y")}', fontsize=16)

            # Plot Scenario 2 - Scenario 1 difference
            ax1 = axes[0]
            try:
                im1 = ax1.contourf(longitude, latitude, ho2_diff_s2_s1, cmap='viridis', levels=np.linspace(0, np.nanmax(ho2_diff_s2_s1), 21), extend='both')
                ax1.coastlines()
                ax1.add_feature(cfeature.BORDERS, linestyle=':')
                ax1.set_title('Scenario 2 - Scenario 1 (ppb)')
                fig.colorbar(im1, ax=ax1, label='ppb', orientation='vertical', shrink=0.7)
            except Exception as e:
                print(f"Error plotting Scenario 2 - Scenario 1 difference for {date_objects[i]}: {e}")
                traceback.print_exc()
                error_occurred = True

            # Plot Scenario 3 - Scenario 1 difference
            ax2 = axes[1]
            try:
                im2 = ax2.contourf(longitude, latitude, ho2_diff_s3_s1, cmap='viridis', levels=np.linspace(0, np.nanmax(ho2_diff_s3_s1), 21), extend='both')
                ax2.coastlines()
                ax2.add_feature(cfeature.BORDERS, linestyle=':')
                ax2.set_title('Scenario 3 - Scenario 1 (ppb)')
                fig.colorbar(im2, ax=ax2, label='ppb', orientation='vertical', shrink=0.7)
            except Exception as e:
                print(f"Error plotting Scenario 3 - Scenario 1 difference for {date_objects[i]}: {e}")
                traceback.print_exc()
                error_occurred = True

            plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to prevent overlap with suptitle

            # --- Save the difference map as a PNG file ---
            filename = f"ho2_difference_level0_{date_objects[i].strftime('%Y%m')}.png"
            plt.savefig(filename)  # Save with a filename indicating the month
            plt.close(fig)  # Close the figure

        else:
            print(f"Skipping difference map for {date_objects[i]} due to missing data.")
else:
    print("Skipping difference maps plot due to errors during data processing or missing lat/lon.")

# --- Part 8: Plotting Percentage Burden Variation at Level 0 ---
# --- Create a bar plot of the percentage variation in burden at Level 0 ---
print("\n--- Part 8: Plotting Percentage Burden Variation at Level 0 ---")

if not error_occurred:
    level_0 = 0
    scenario_1_burden_level0 = np.array(ho2_burden_monthly_kg['Scenario 1'][level_0])
    scenario_2_burden_level0 = np.array(ho2_burden_monthly_kg['Scenario 2'][level_0])
    scenario_3_burden_level0 = np.array(ho2_burden_monthly_kg['Scenario 3'][level_0])

    # Convert monthly date strings to datetime objects for x-axis labels
    date_objects = [date(int(d[:4]), int(d[4:6]), 1) for d in monthly_dates]
    month_labels = [d.strftime('%b-%Y') for d in date_objects]
    x = np.arange(len(month_labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(figsize=(16, 8))
    rects1 = ax.bar(x - width / 2,
                    (scenario_2_burden_level0 - scenario_1_burden_level0) / scenario_1_burden_level0 * 100,
                    width, label='Scenario 2 vs Scenario 1')
    rects2 = ax.bar(x + width / 2,
                    (scenario_3_burden_level0 - scenario_1_burden_level0) / scenario_1_burden_level0 * 100,
                    width, label='Scenario 3 vs Scenario 1')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Percentage Variation (%)')
    ax.set_xlabel('Month')
    ax.set_title('Percentage Variation of HO$_2$ Burden at Level 0 Relative to Scenario 1')
    ax.set_xticks(x)
    ax.set_xticklabels(month_labels, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', linestyle='--')

    fig.tight_layout()
    plt.savefig('ho2_burden_percentage_variation_level0.png')
    plt.close()

else:
    print("Skipping percentage burden variation plot due to errors during data processing.")
    

# --- Part 9: Calculate and Plot Overlapped Vertical Profiles for All Scenarios and Cities (Monthly) ---
print("\n--- Part 9: Calculate and Plot Overlapped Vertical Profiles for All Scenarios and Cities (Monthly) ---")

# Configuration for vertical profiles
profile_locations = {
    "Sydney": (-33.87, 151.21),
    "New York City": (40.71, -74.01),
    "Tokyo": (35.68, 139.76),
}
city_names = list(profile_locations.keys())
scenario_names_list = list(scenario_directories.keys())
scenario_colors = {'Scenario 1': 'blue', 'Scenario 2': 'red', 'Scenario 3': 'green'}
line_styles = {'Sydney': '-', 'New York City': '--', 'Tokyo': ':'}

# Constants
AVOGADRO = 6.02214076e23  # molecules/mol
M_AIR = 28.97e-3          # kg/mol (molar mass of air)
M3_TO_CM3 = 1e6           # Conversion factor from m^3 to cm^3

# Get unique months from the monthly_dates list
unique_months = sorted(list(set([d[:6] for d in monthly_dates])))

# Loop through each unique month
for month_str in unique_months:
    plt.figure(figsize=(12, 10))
    all_profiles_plotted_this_month = False
    handles = []
    labels = []

    # Loop through each scenario
    for scenario_name, base_dir in scenario_directories.items():
        species_conc_base = species_conc_base_filenames[scenario_name]
        state_met_base = state_met_base_filenames[scenario_name]

        # Construct the date string for the current month (using the first day)
        current_month_dates = [d for d in monthly_dates if d[:6] == month_str]
        if not current_month_dates:
            print(f"Warning: No data found for {scenario_name} in {month_str}.")
            continue
        date_str = current_month_dates[0] # Use the first day of the month

        try:
            species_conc_filepath = os.path.join(base_dir, f'{species_conc_base}{date_str}_0000z.nc4')
            state_met_filepath = os.path.join(base_dir, f'{state_met_base}{date_str}_0000z.nc4')

            if not os.path.exists(species_conc_filepath) or not os.path.exists(state_met_filepath):
                print(f"Warning: Skipping profiles for {scenario_name} - {month_str} due to missing files.")
                continue

            with xr.open_dataset(species_conc_filepath) as ds_species, \
                 xr.open_dataset(state_met_filepath) as ds_state_met:

                ho2_vvmr = ds_species['SpeciesConcVV_HO2'].isel(time=0)
                air_density = ds_state_met['Met_AIRDEN'].isel(time=0)
                pressure = ds_state_met['Met_PMID'].isel(time=0)
                lat = ds_species['lat']
                lon = ds_species['lon']
                lev = ds_species['lev']

                # Loop through the specified profile locations
                for city, (profile_lat, profile_lon) in profile_locations.items():
                    # Find the closest grid cell
                    lat_diff = np.abs(lat - profile_lat)
                    lon_diff = np.abs(lon - profile_lon)
                    lat_idx = lat_diff.argmin()
                    lon_idx = lon_diff.argmin()

                    ho2_vvmr_profile = ho2_vvmr[:, lat_idx, lon_idx].values
                    air_density_profile = air_density[:, lat_idx, lon_idx].values
                    pressure_profile = pressure[:, lat_idx, lon_idx].values

                    # Calculate HO2 number density (molecules/cm³)
                    ho2_num_density = ho2_vvmr_profile * (air_density_profile / M_AIR) * AVOGADRO / M3_TO_CM3

                    # Plot the profile for the current city and scenario
                    label = f'{city} - {scenario_name}'
                    line, = plt.plot(ho2_num_density, pressure_profile, color=scenario_colors[scenario_name], linestyle=line_styles[city], label=label)
                    handles.append(line)
                    labels.append(label)
                    all_profiles_plotted_this_month = True

        except Exception as e:
            print(f"Error processing overlapped vertical profiles for {scenario_name} - {month_str}: {e}")
            traceback.print_exc()

    # Add labels and title for the monthly plot
    if all_profiles_plotted_this_month:
        try:
            month_date = date.fromisoformat(f'{month_str[:4]}-{month_str[4:]}-01')
            plt.xlabel('HO$_2$ Concentration (molecules/cm³)')
            plt.ylabel('Pressure (hPa)')
            plt.title(f'Vertical Profiles of HO$_2$ - {month_date.strftime("%B %Y")}\n(All Scenarios and Cities)')
            plt.gca().invert_yaxis()  # Invert y-axis for pressure (higher at bottom)
            plt.grid(True)

            # Sort handles and labels by city name
            sorted_indices = sorted(range(len(labels)), key=lambda k: labels[k].split(' - ')[0])
            sorted_handles = [handles[i] for i in sorted_indices]
            sorted_labels = [labels[i] for i in sorted_indices]

            plt.legend(sorted_handles, sorted_labels, loc='best')  # Add the ordered legend
            plt.tight_layout()

            # Save the combined monthly plot
            filename = f'ho2_vertical_profiles_all_scenarios_ordered_legend_{month_str}.png'
            plt.savefig(filename)
            print(f"Combined vertical profiles for {month_date.strftime('%B %Y')} (all scenarios, ordered legend) saved successfully as {filename}")
        except ValueError as ve:
            print(f"Error formatting date for title or filename for month {month_str}: {ve}")
        except Exception as e:
            print(f"Error saving combined vertical profiles for {month_date.strftime('%B %Y')} (all scenarios, ordered legend): {e}")
            traceback.print_exc()
    else:
        print(f"No vertical profile data to plot for {month_str}.")
    plt.close()

print("\n--- End of Part 9: Overlapped Vertical Profile Calculation for All Scenarios and Cities (Monthly) with Ordered Legend ---")

# --- Part 10: Latitude-Pressure Cross-Section at Sydney's Longitude (Monthly) with Fixed Y-axis (Corrected Level Interpretation) ---
print("\n--- Part 10: Latitude-Pressure Cross-Section at Sydney's Longitude (Monthly) with Fixed Y-axis (Corrected Level Interpretation) ---")

# Configuration for cross-section
sydney_longitude = 151.21
selected_scenario = 'Scenario 1'  # Choose the scenario for the cross-section
top_level_index = 0
bottom_level_index = 71

# Constants
AVOGADRO = 6.02214076e23  # molecules/mol
M_AIR = 28.97e-3          # kg/mol (molar mass of air)
M3_TO_CM3 = 1e6           # Conversion factor from m^3 to cm^3

# Get unique months from the monthly_dates list
unique_months = sorted(list(set([d[:6] for d in monthly_dates])))

# Determine pressure at level 0 and 71 from the first available state met file
pressure_at_top = np.nan
pressure_at_bottom = np.nan

first_month = unique_months[0]
base_dir = scenario_directories.get(selected_scenario)
state_met_base = state_met_base_filenames.get(selected_scenario)

if base_dir and state_met_base:
    current_month_dates = [d for d in monthly_dates if d[:6] == first_month]
    if current_month_dates:
        date_str = current_month_dates[0]
        state_met_filepath = os.path.join(base_dir, f'{state_met_base}{date_str}_0000z.nc4')
        if os.path.exists(state_met_filepath):
            try:
                with xr.open_dataset(state_met_filepath) as ds_state_met:
                    pressure = ds_state_met['Met_PMID']
                    lon = ds_state_met['lon']
                    lon_diff = np.abs(lon - sydney_longitude)
                    lon_idx = lon_diff.argmin()

                    if pressure.shape[1] > top_level_index:
                        pressure_at_top = pressure[0, top_level_index, :, lon_idx].values.max()
                    if pressure.shape[1] > bottom_level_index:
                        pressure_at_bottom = pressure[0, bottom_level_index, :, lon_idx].values.min()

                print(f"Pressure at level {bottom_level_index}: {pressure_at_bottom:.2f} hPa")
                print(f"Pressure at level {top_level_index}: {pressure_at_top:.2f} hPa")

            except Exception as e:
                print(f"Error reading pressure data to determine fixed y-axis limits: {e}")

# Loop through each unique month for plotting
for month_str in unique_months:
    plt.figure(figsize=(10, 8))
    cross_section_plotted_this_month = False

    base_dir = scenario_directories.get(selected_scenario)
    species_conc_base = species_conc_base_filenames.get(selected_scenario)
    state_met_base = state_met_base_filenames.get(selected_scenario)

    if not base_dir or not species_conc_base or not state_met_base:
        print(f"Warning: Scenario '{selected_scenario}' not found in configuration for cross-section.")
        continue

    # Construct the date string for the current month (using the first day)
    current_month_dates = [d for d in monthly_dates if d[:6] == month_str]
    if not current_month_dates:
        print(f"Warning: No data found for {selected_scenario} in {month_str}.")
        continue
    date_str = current_month_dates[0] # Use the first day of the month

    try:
        species_conc_filepath = os.path.join(base_dir, f'{species_conc_base}{date_str}_0000z.nc4')
        state_met_filepath = os.path.join(base_dir, f'{state_met_base}{date_str}_0000z.nc4')

        if not os.path.exists(species_conc_filepath) or not os.path.exists(state_met_filepath):
            print(f"Warning: Skipping cross-section for {selected_scenario} - {month_str} due to missing files.")
            continue

        with xr.open_dataset(species_conc_filepath) as ds_species, \
             xr.open_dataset(state_met_filepath) as ds_state_met:

            ho2_vvmr = ds_species['SpeciesConcVV_HO2'].isel(time=0)
            air_density = ds_state_met['Met_AIRDEN'].isel(time=0)
            pressure = ds_state_met['Met_PMID'].isel(time=0)
            lat = ds_species['lat']
            lon = ds_species['lon']
            lev = ds_species['lev']

            # Find the longitude index closest to Sydney's longitude
            lon_diff = np.abs(lon - sydney_longitude)
            lon_idx = lon_diff.argmin()

            ho2_vvmr_slice = ho2_vvmr[:, :, lon_idx].values  # (lev, lat)
            air_density_slice = air_density[:, :, lon_idx].values  # (lev, lat)
            pressure_slice = pressure[:, :, lon_idx].values  # (lev, lat)
            latitude = lat.values

            # Calculate HO2 number density (molecules/cm³)
            ho2_num_density = ho2_vvmr_slice * (air_density_slice / M_AIR) * AVOGADRO / M3_TO_CM3

            # Create the plot
            lat_mesh, p_mesh = np.meshgrid(latitude, pressure_slice[:, 0]) # Use pressure values for y-axis

            contour = plt.contourf(lat_mesh, p_mesh, ho2_num_density, cmap='viridis', levels=np.linspace(np.nanmin(ho2_num_density), np.nanmax(ho2_num_density), 20))
            cbar = plt.colorbar(contour, label='HO$_2$ Concentration (molecules/cm³)')
            plt.xlabel('Latitude (°N)')
            plt.ylabel('Pressure (hPa)')
            month_date = date.fromisoformat(f'{month_str[:4]}-{month_str[4:]}-01')
            plt.title(f'HO$_2$ Cross-Section at Longitude {lon.values[lon_idx]:.2f}°E (near Sydney)\n{month_date.strftime("%B %Y")} - {selected_scenario}')
            plt.gca().invert_yaxis()  # Invert y-axis for pressure (higher at bottom)
            if np.isfinite(pressure_at_bottom) and np.isfinite(pressure_at_top):
                plt.ylim(pressure_at_bottom, pressure_at_top) # Set fixed y-axis limits (higher pressure at bottom)
            plt.tight_layout()
            cross_section_plotted_this_month = True

            # Save the plot
            filename = f'ho2_cross_section_sydney_lon_fixed_yaxis_levels_corrected_{selected_scenario.replace(" ", "_")}_{month_str}.png'
            plt.savefig(filename)
            print(f"Cross-section (fixed y-axis levels 71-0) for {month_date.strftime('%B %Y')} ({selected_scenario}) saved as {filename}")

    except Exception as e:
        print(f"Error processing cross-section (fixed y-axis levels 71-0) for {selected_scenario} - {month_str}: {e}")
        traceback.print_exc()
    finally:
        plt.close()

print("\n--- End of Part 10: Latitude-Pressure Cross-Section at Sydney's Longitude (Monthly) with Fixed Y-axis (Corrected Level Interpretation) ---")
