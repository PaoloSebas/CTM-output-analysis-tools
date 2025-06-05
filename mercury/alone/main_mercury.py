# main_mercury.py

import xarray as xr
import os
import sys
from file_io_manager import find_geoschem_nc4_files, open_geoschem_datasets
from mercury_analysis_2 import quantify_total_deposition, quantify_hg2_to_ssa_transfer
from showing_results import display_mercury_results
from plotting import plot_variable_global # Import the plotting function

def main():
    print("================================================================================")
    print("      --- GEOS-Chem Mercury Deposition & Transfer Analysis ---      ")
    print("================================================================================")
    print("Welcome to the GEOS-Chem Mercury Analysis Tool. This program will guide you through:")
    print("1. Specifying a folder to automatically load GEOS-Chem NetCDF files.")
    print("2. Performing a combined mercury-related calculation (Deposition and Transfer).")
    print("3. Displaying numerical results.")
    print("4. Generating global maps of averaged rates (new!).\n") # Update welcome message

    # --- Step 1: User Input for Data Folder ---
    print("--- Step 1: Specify Data Folder ---")
    print("Please provide the path to the folder containing your GEOS-Chem .nc4 output files.")
    print("Example: '/path/to/my/geoschem_output_data/' or 'data/' (for a local subfolder)")

    folder_path = input("Enter the folder path: ").strip()

    # Validate folder path
    if not os.path.isdir(folder_path):
        print(f"\nError: Folder not found at '{folder_path}'. Please ensure the path is correct and accessible.")
        print("Exiting application.")
        return

    # --- Step 2: Automatic File Discovery and Loading ---
    print(f"\n--- Step 2: Discovering and Loading GEOS-Chem Files from '{folder_path}' ---")
    print("The tool will now automatically scan the specified folder and attempt to open")
    print("NetCDF files recognized as GEOS-Chem output using xarray.")

    opened_datasets = open_geoschem_datasets(folder_path)

    if not opened_datasets:
        print("\nNo GEOS-Chem .nc4 datasets were successfully loaded from the specified folder.")
        print("Please ensure your files are present and follow the 'GEOSChem.<description>.<date>.nc4' naming convention.")
        print("Exiting application.")
        return

    print("\nSuccessfully loaded the following GEOS-Chem datasets:")
    for desc, ds in opened_datasets.items():
        print(f"  - '{desc}' (Variables: {list(ds.data_vars.keys())})")
    print("-" * 80)

    # Initialize results variables to None
    total_dry_deposition_kg_s = None
    total_wet_deposition_kg_s = None
    hg2_to_ssa_transfer_rate_kg_s = None

    # Retrieve all potentially needed datasets upfront
    ds_mercurychem = opened_datasets.get('MercuryChem')
    ds_statemet = opened_datasets.get('StateMet')
    ds_drydep = opened_datasets.get('DryDep')      # Get DryDep dataset
    ds_wetconv = opened_datasets.get('WetLossConv') # Get WetLossConv dataset
    ds_wetls = opened_datasets.get('WetLossLS')     # Get WetLossLS dataset

    # --- Step 3: Run the combined analysis directly ---
    print("\n--- Running: Combined Wet/Dry Deposition AND Sea Salt Transfer Analysis ---")

    # Run deposition if datasets are available
    required_dep_ds = ['MercuryChem', 'StateMet', 'DryDep', 'WetLossConv', 'WetLossLS']
    if not all(opened_datasets.get(ds_name) is not None for ds_name in required_dep_ds):
        print(f"\nWARNING: Not all required datasets for total deposition analysis were found ({', '.join(required_dep_ds)}).")
        print("Skipping Total Wet and Dry Deposition calculation.")
        total_dry_deposition_kg_s = None
        total_wet_deposition_kg_s = None
    else:
        total_dry_deposition_kg_s, total_wet_deposition_kg_s = quantify_total_deposition(
            ds_mercurychem, ds_statemet, ds_drydep, ds_wetconv, ds_wetls
        )

    # Run transfer if datasets are available
    required_transfer_ds = ['MercuryChem', 'StateMet']
    if not all(opened_datasets.get(ds_name) is not None for ds_name in required_transfer_ds):
        print(f"\nWARNING: Required datasets for Hg(II) to SSA transfer analysis not fully loaded ({', '.join(required_transfer_ds)}).")
        print("Skipping Hg(II) Gas to Sea Salt Aerosol Transfer calculation.")
        hg2_to_ssa_transfer_rate_kg_s = None
    else:
        hg2_to_ssa_transfer_rate_kg_s = quantify_hg2_to_ssa_transfer(ds_mercurychem, ds_statemet)
            
    # Display combined results
    display_mercury_results(
        total_dry_deposition_kg_s,
        total_wet_deposition_kg_s,
        hg2_to_ssa_transfer_rate_kg_s
    )

    # --- Step 4: Generate Global Plots of Averaged Rates ---
    print("\n--- Step 4: Generating Global Maps of Averaged Rates ---")

    # Helper function to prepare data for plotting
    def prepare_for_global_plot(data_array: xr.DataArray, process_name: str) -> xr.DataArray | None:
        if data_array is None:
            print(f"  - Skipping plot for {process_name}: Data not available.")
            return None

        # Average over time
        averaged_da = data_array.mean(dim='time', keep_attrs=True) # Keep attributes for units/long_name

        # If there's a 'lev' dimension, sum over it (to get total column rate for 2D map)
        if 'lev' in averaged_da.dims:
            print(f"  - Summing '{process_name}' over 'lev' dimension for plotting.")
            averaged_da = averaged_da.sum(dim='lev', keep_attrs=True)

        print(f"  - Prepared '{process_name}' data for plotting. Shape: {averaged_da.shape}")
        return averaged_da

    # Prepare and plot Dry Deposition
    avg_dry_dep_for_plot = prepare_for_global_plot(total_dry_deposition_kg_s, "Dry Deposition")
    if avg_dry_dep_for_plot is not None:
        plot_variable_global(avg_dry_dep_for_plot, var_name="Average Dry Deposition Rate")

    # Prepare and plot Wet Deposition
    avg_wet_dep_for_plot = prepare_for_global_plot(total_wet_deposition_kg_s, "Wet Deposition")
    if avg_wet_dep_for_plot is not None:
        plot_variable_global(avg_wet_dep_for_plot, var_name="Average Wet Deposition Rate")

    # Prepare and plot Hg(II) to SSA Transfer
    avg_hg2_ssa_for_plot = prepare_for_global_plot(hg2_to_ssa_transfer_rate_kg_s, "Hg(II) Gas to Sea Salt Aerosol Transfer")
    if avg_hg2_ssa_for_plot is not None:
        plot_variable_global(avg_hg2_ssa_for_plot, var_name="Average Hg(II) to Sea Salt Aerosol Transfer Rate")

    print("\n--- Global Map Generation Complete ---")
    print("-" * 80)


    # --- Final Cleanup: Close all opened datasets ---
    print("\n--- Final Step: Cleaning Up Resources ---")
    print("Closing all opened xarray datasets to free up memory.")
    for desc, ds in opened_datasets.items():
        if ds is not None:
            ds.close()
            print(f"  - Dataset '{desc}' closed.")
    print("Resource cleanup complete.")
    print("\nExiting the application. Thank you for using the tool!")


if __name__ == "__main__":
    main()
