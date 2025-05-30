# main_3.py
import os
from datetime import datetime
import xarray as xr
import numpy as np # Imported for potential numerical operations on results in main.py
import matplotlib.pyplot as plt # Kept for general plotting needs if any, though actual plots done by plotting.py

# Import core file I/O and variable analysis modules
from file_io_manager import find_geoschem_nc4_files, display_file_summary_table, open_geoschem_datasets
from variable_analyzer import analyze_and_list_variables, display_variable_summary_table, export_variables_to_csv

# Import mercury analysis module
from mercury_analysis import quantify_total_deposition, quantify_hg2_to_ssa_transfer

# Import the refactored plotting module
from plotting import plot_variable_global, plot_variable_antarctic

# Note: gif_creator.py functions are not yet refactored to accept xarray.DataArray
# If they were, you would import them here:
# from gif_creator import create_antarctic_gif, create_macquarie_gif


def display_main_menu():
    """Displays the main menu options to the user."""
    print("\n--- GEOS-Chem Data Analysis Tool ---")
    print("1. Load & Process GEOS-Chem Data Folder")
    print("2. Perform Mercury Analysis")
    print("3. Perform Time-Series Analysis (Future)")
    print("4. Plotting a selected variable")
    print("5. Conduct Spatial Aggregations (Future)")
    print("6. Compare GEOS-Chem Runs (Future)")
    print("7. Exit")
    print("------------------------------------")

def process_data_folder():
    """
    Handles the sequence for prompting a path, analyzing files, opening them,
    and listing variables.
    """
    print("\n--- Load & Process GEOS-Chem Data Folder ---")
    folder_path = input("Please provide the path that contains the nc4 files to be processed: ").strip()

    if not folder_path:
        print("No path provided. Returning to main menu.")
        return None, None, None # Return None for datasets, files_info, and folder_path

    # Validate if folder_path is a valid directory
    if not os.path.isdir(folder_path):
        print(f"Error: The provided path '{folder_path}' is not a valid directory or is inaccessible.")
        return None, None, None

    # Step 1: Find and display summary of files
    print(f"Scanning folder: {folder_path}")
    files_info = find_geoschem_nc4_files(folder_path)
    if not files_info:
        print("No GEOS-Chem .nc4 files found in the specified folder.")
        return None, None, None

    display_file_summary_table(files_info)

    # Step 2: Open all files with Xarray
    opened_datasets = open_geoschem_datasets(folder_path)
    if not opened_datasets:
        print("Failed to open any datasets. Please check file integrity and folder path.")
        # It's safer to close only if opened_datasets is a dict and not None
        if isinstance(opened_datasets, dict):
            for ds in opened_datasets.values():
                ds.close()
        return None, None, None

    # Step 3: Analyze and list all variables
    all_variables_info = analyze_and_list_variables(opened_datasets)
    if not all_variables_info:
        print("No variables found in the opened datasets.")
        # Ensure all datasets are closed if no variables found
        for ds in opened_datasets.values():
            ds.close()
        return None, None, None

    display_variable_summary_table(all_variables_info)

    # Step 4: Ask user to export to CSV
    if all_variables_info:
        save_csv_choice = input("Do you want to export this variable list to a CSV file? (yes/no): ").strip().lower()
        if save_csv_choice == 'yes':
            export_variables_to_csv(all_variables_info)
        else:
            print("CSV export skipped.")

    print("\nFiles processed and variables listed. Datasets remain open for further analysis.")
    print("You can now proceed with other operations using the loaded data.")

    return opened_datasets, files_info, folder_path

def perform_mercury_analysis(opened_datasets: dict[str, xr.Dataset], files_info: list[dict], folder_path: str):
    """
    Orchestrates the mercury-specific analysis using the loaded datasets and file paths.
    """
    print("\n--- Perform Mercury Analysis ---")

    if opened_datasets is None or not opened_datasets:
        print("No datasets are currently loaded. Please use '1. Load & Process Data Folder' first.")
        return

    # Map descriptions to their full file paths for mercury_analysis_new.py
    file_paths_by_description = {}
    if files_info and folder_path:
        for file_info in files_info:
            desc = file_info['description']
            filename = file_info['filename']
            full_path = os.path.join(folder_path, filename)
            # Take the first path found for each description.
            if desc not in file_paths_by_description:
                file_paths_by_description[desc] = full_path
    else:
        print("Error: File information or folder path is missing. Cannot proceed with mercury analysis.")
        return

    required_mercury_files = {
        'WetLossConv': None,
        'WetLossLS': None,
        'DryDep': None,
        'MercuryChem': None,
        'StateMet': None
    }

    # Check if all required files are present and get their paths
    all_files_present = True
    for desc in required_mercury_files.keys():
        if desc in file_paths_by_description:
            required_mercury_files[desc] = file_paths_by_description[desc]
        else:
            print(f"Missing required file type for mercury analysis: '{desc}'.")
            all_files_present = False

    if not all_files_present:
        print("Not all required files for mercury analysis were found. Returning.")
        return

    print("\n--- Quantifying Total Wet and Dry Deposition ---")
    # --- TEMPORARY: REMOVE TRY-EXCEPT TO SEE FULL TRACEBACK ---
    # try:
    # ... (file path assignments - keep these) ...
    total_wet_dep, total_dry_dep = quantify_total_deposition(
        required_mercury_files['WetLossConv'],
        required_mercury_files['WetLossLS'],
        required_mercury_files['DryDep']
    )
    if total_wet_dep is not None and total_dry_dep is not None:
        # Assuming results are xarray DataArrays, sum the first time step
        if 'time' in total_wet_dep.dims and total_wet_dep['time'].size > 0:
            print(f"Total Wet Deposition calculated. First time step sum (kg_Hg s-1): {total_wet_dep.isel(time=0).sum().item():.2e}")
        else:
            print(f"Total Wet Deposition calculated. Overall sum (kg_Hg s-1): {total_wet_dep.sum().item():.2e}")

        if 'time' in total_dry_dep.dims and total_dry_dep['time'].size > 0:
            print(f"Total Dry Deposition calculated. First time step sum (kg_Hg s-1): {total_dry_dep.isel(time=0).sum().item():.2e}")
        else:
            print(f"Total Dry Deposition calculated. Overall sum (kg_Hg s-1): {total_dry_dep.sum().item():.2e}")

        print("Results are xarray.DataArray objects.")
    else:
        print("Total deposition calculation failed (quantify_total_deposition returned None).") # Updated message for clarity
    # except Exception as e:
    #     print(f"An error occurred during total deposition calculation: {e}")
    # --- END TEMPORARY CHANGE ---

    print("\n--- Quantifying Hg2 (RGM) Transfer to Sea Salt Aerosol (SSA) ---")
    try:
        # Calling mercury_analysis_new functions with file paths.
        # IMPORTANT: The functions in 'mercury_analysis_new.py' will re-open these files.
        # For efficiency, 'mercury_analysis_new.py' should ideally be refactored
        # to accept the 'xarray.Dataset' objects directly from 'opened_datasets'.
        hg2_to_ssa_transfer = quantify_hg2_to_ssa_transfer(
            required_mercury_files['MercuryChem'],
            required_mercury_files['StateMet']
        )
        if hg2_to_ssa_transfer is not None:
            if 'time' in hg2_to_ssa_transfer.dims and hg2_to_ssa_transfer['time'].size > 0:
                print(f"Hg2 to SSA Transfer calculated. First time step sum (kg_Hg s-1): {hg2_to_ssa_transfer.isel(time=0).sum().item():.2e}")
            else:
                print(f"Hg2 to SSA Transfer calculated. Overall sum (kg_Hg s-1): {hg2_to_ssa_transfer.sum().item():.2e}")
            print("Results are an xarray.DataArray object.")
        else:
            print("Hg2 to SSA transfer calculation failed.")
    except Exception as e:
        print(f"An error occurred during Hg2 to SSA transfer calculation: {e}")

    print("\nMercury analysis complete. Remember to check the console for calculation summaries.")


def plot_selected_variable(opened_datasets: dict[str, xr.Dataset]):
    """
    Allows the user to select a variable from loaded datasets and choose a plotting option.
    Inspired by the plotting submenu from the old main code.
    """
    print("\n--- Plotting a selected variable ---")
    if opened_datasets is None or not opened_datasets:
        print("No datasets are currently loaded. Please use '1. Load & Process Data Folder' first.")
        return

    # Create a mapping of display index to (dataset_key, variable_name, xarray.DataArray)
    available_vars_map = {}
    display_idx = 1
    print("\nAvailable variables across all loaded datasets:")
    for ds_key, ds_obj in opened_datasets.items():
        for var_name in ds_obj.data_vars:
            available_vars_map[display_idx] = (ds_key, var_name, ds_obj[var_name])
            # Display variable name and units if available
            units = ds_obj[var_name].attrs.get('units', 'unitless')
            print(f"{display_idx}. {var_name} (from {ds_key}, units: {units})")
            display_idx += 1

    if not available_vars_map:
        print("No variables found in the loaded datasets to plot.")
        return

    selected_data_array: xr.DataArray = None
    selected_var_name: str = ""

    while selected_data_array is None:
        try:
            var_choice = input("Select a variable number to plot (or 0 to go back): ").strip()
            if var_choice == '0':
                return # Go back to main menu
            var_index = int(var_choice)

            if 1 <= var_index <= len(available_vars_map):
                _ds_key, selected_var_name, selected_data_array = available_vars_map[var_index]
                print(f"Selected variable: {selected_var_name} from dataset '{_ds_key}'")
            else:
                print("Invalid variable number. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

    # data_to_plot_slice will be the 2D (lat, lon) slice used for static plots
    # selected_data_array (original) will be used for GIF creation across time
    data_to_plot_slice = selected_data_array

    # Handle time dimension if present
    if 'time' in data_to_plot_slice.dims and data_to_plot_slice['time'].size > 1:
        print(f"\nSelected variable '{selected_var_name}' has a time dimension with {data_to_plot_slice['time'].size} steps.")
        print("Please choose a time step for the static plot:")
        print("1. First time step")
        print("2. Last time step")
        print("3. Specific time index (0 to N-1)")
        # For GIFs, the full time series (selected_data_array) will be used later
        
        time_choice_made = False
        while not time_choice_made:
            time_choice = input("Enter your choice (1-3): ").strip()
            if time_choice == '1':
                data_to_plot_slice = data_to_plot_slice.isel(time=0)
                time_choice_made = True
            elif time_choice == '2':
                data_to_plot_slice = data_to_plot_slice.isel(time=-1)
                time_choice_made = True
            elif time_choice == '3':
                try:
                    idx = int(input(f"Enter time index (0 to {data_to_plot_slice['time'].size - 1}): ").strip())
                    if 0 <= idx < data_to_plot_slice['time'].size:
                        data_to_plot_slice = data_to_plot_slice.isel(time=idx)
                        time_choice_made = True
                    else:
                        print("Index out of range. Please try again.")
                except ValueError:
                    print("Invalid input. Please enter a number.")
            else:
                print("Invalid choice. Please try again.")
    elif 'time' in data_to_plot_slice.dims and data_to_plot_slice['time'].size == 1:
        print(f"Variable '{selected_var_name}' has only one time step. Using it directly for static plot.")
        data_to_plot_slice = data_to_plot_slice.isel(time=0) # Ensure it's a 2D slice if possible

    # Handle level dimension if present (for 3D variables)
    if 'lev' in data_to_plot_slice.dims and data_to_plot_slice['lev'].size > 1:
        print(f"\nSelected variable '{selected_var_name}' has a 'lev' (level) dimension with {data_to_plot_slice['lev'].size} levels.")
        print("Please choose a level for the static plot:")
        # Display available levels if reasonable
        if data_to_plot_slice['lev'].size < 10:
            print(f"Available levels (indices 0 to {data_to_plot_slice['lev'].size - 1}): {data_to_plot_slice['lev'].values}")
        else:
            print(f"Levels range from {data_to_plot_slice['lev'].values[0]} to {data_to_plot_slice['lev'].values[-1]} (indices 0 to {data_to_plot_slice['lev'].size - 1}).")

        level_choice_made = False
        while not level_choice_made:
            try:
                lev_choice = input("Enter level index (0 to N-1): ").strip()
                lev_idx = int(lev_choice)
                if 0 <= lev_idx < data_to_plot_slice['lev'].size:
                    data_to_plot_slice = data_to_plot_slice.isel(lev=lev_idx)
                    level_choice_made = True
                else:
                    print("Level index out of range. Please try again.")
            except ValueError:
                print("Invalid input. Please enter a number.")
    elif 'lev' in data_to_plot_slice.dims and data_to_plot_slice['lev'].size == 1:
        print(f"Variable '{selected_var_name}' has only one level. Using it directly for static plot.")
        data_to_plot_slice = data_to_plot_slice.isel(lev=0) # Ensure it's a 2D slice if possible


    # Plotting options sub-menu
    while True:
        print("\n--- Plotting options for selected variable ---")
        print("1. Global Distribution Plot")
        print("2. Antarctic Plot")
        print("3. Antarctic GIF Animation (Requires gif_creator.py refactor)")
        print("4. Macquarie Island GIF Animation (Requires gif_creator.py refactor)")
        print("0. Back to variable selection")

        plot_choice = input("Enter your choice: ").strip()

        if plot_choice == '0':
            break # Go back to variable selection
        elif plot_choice == '1':
            # Call the actual refactored function from plotting.py
            plot_variable_global(data_to_plot_slice, selected_var_name)
        elif plot_choice == '2':
            # Call the actual refactored function from plotting.py
            plot_variable_antarctic(data_to_plot_slice, selected_var_name)
        elif plot_choice == '3':
            # This would call the refactored GIF function, which needs the original (unsliced) DataArray
            # because it will iterate through time steps itself.
            # E.g.: create_antarctic_gif(selected_data_array, selected_var_name)
            print("\n--- GIF Creation Placeholder ---")
            print("Antarctic GIF creation requires 'gif_creator.py' to be refactored to accept xarray.DataArray.")
            print("Ensure the variable has a 'time' dimension for animation.")
        elif plot_choice == '4':
            # E.g.: create_macquarie_gif(selected_data_array, selected_var_name)
            print("\n--- GIF Creation Placeholder ---")
            print("Macquarie Island GIF creation requires 'gif_creator.py' to be refactored to accept xarray.DataArray.")
            print("Ensure the variable has a 'time' dimension for animation.")
        else:
            print("Invalid choice. Please try again.")


def handle_future_feature(feature_name: str):
    """Generic handler for features not yet implemented."""
    print(f"\n[Feature coming soon]: {feature_name}...")

def main():
    """Main function to run the GEOS-Chem data analysis tool."""
    print("Welcome to the new GEOS-Chem Data Analysis Tool!")

    current_opened_datasets = None
    current_files_info = None
    current_folder_path = None

    while True:
        display_main_menu()
        choice = input("Enter your choice (1-7): ").strip()

        if choice == '1':
            current_opened_datasets, current_files_info, current_folder_path = process_data_folder()
        elif choice == '2':
            perform_mercury_analysis(current_opened_datasets, current_files_info, current_folder_path)
        elif choice == '3':
            handle_future_feature("Perform Time-Series Analysis")
        elif choice == '4': # NEW CHOICE: Plotting
            plot_selected_variable(current_opened_datasets)
        elif choice == '5':
            handle_future_feature("Conduct Spatial Aggregations")
        elif choice == '6':
            handle_future_feature("Compare GEOS-Chem Runs")
        elif choice == '7': # Updated exit number
            print("Exiting the application. Goodbye!")
            # Ensure datasets are closed upon exit
            if current_opened_datasets:
                print("Closing opened datasets...")
                for ds in current_opened_datasets.values():
                    ds.close()
            break
        else:
            print("Invalid choice. Please enter a number between 1 and 7.")

if __name__ == "__main__":
    main()
