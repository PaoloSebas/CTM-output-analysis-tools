# main_mercury.py
import os
import xarray as xr
from datetime import datetime

# Import modules for file management, analysis, and result display
from file_io_manager import find_geoschem_nc4_files, open_geoschem_datasets
from mercury_analysis_2 import quantify_total_deposition, quantify_hg2_to_ssa_transfer
from showing_results import display_mercury_results, offer_and_execute_plotting

def run_mercury_analysis_tool():
    """
    Orchestrates the GEOS-Chem mercury analysis workflow.
    This didactic script guides the user through file loading,
    analysis selection, result display, and optional plotting.
    """
    print("\n" + "="*80)
    print("      GEOS-Chem Mercury Analysis Tool (Didactic Version)".center(80))
    print("="*80)
    print("Welcome to the GEOS-Chem Mercury Analysis Tool. This program will guide you through:")
    print("1. Automatically identifying and loading relevant GEOS-Chem NetCDF files.")
    print("2. Performing specific mercury-related calculations.")
    print("3. Displaying numerical results.")
    print("4. Optionally visualizing results through plots.\n")

    # --- Step 1: User Input for Data Folder ---
    print("--- Step 1: Specify Data Folder ---")
    print("Please provide the path to the folder containing your GEOS-Chem .nc4 output files.")
    print("Example: '/path/to/my/geoschem_output_data/' or 'data/' (for a local folder)")

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

    # Find and open datasets using the file_io_manager
    # Note: open_geoschem_datasets groups by description, so we get a dict of datasets.
    opened_datasets = open_geoschem_datasets(folder_path)

    if not opened_datasets:
        print("\nNo GEOS-Chem .nc4 datasets were successfully loaded. Please ensure your files")
        print("are present and follow the 'GEOSChem.<description>.<date>.nc4' naming convention.")
        print("Exiting application.")
        return

    print("\nSuccessfully loaded the following GEOS-Chem datasets:")
    for desc, ds in opened_datasets.items():
        print(f"  - '{desc}' (Variables: {list(ds.data_vars.keys())})")
    print("-" * 80)

    # --- Step 3: Mercury Analysis Menu ---
    print("\n--- Step 3: Select Mercury Analysis ---")
    print("This tool focuses on specific mercury-related calculations.")
    print("Select an option below to proceed with the analysis.")

    menu_options = {
        "1": "Quantify Total Wet and Dry Deposition of Mercury",
        "2": "Quantify Gaseous Hg(II) to Sea Salt Aerosol (SSA) Transfer Rate",
        "0": "Exit Program"
    }

    while True:
        print("\nMercury Analysis Options:")
        for key, value in menu_options.items():
            print(f"  [{key}] {value}")

        choice = input("Enter your choice (e.g., 1 or 2): ").strip()

        wet_dep_result = None
        dry_dep_result = None
        hg2_ssa_transfer_result = None
        analysis_performed = False

        if choice == '1':
            print(f"\nInitiating: {menu_options['1']}...")
            # We need to pass the specific datasets required by quantify_total_deposition
            # These are identified by their 'description' keys from open_geoschem_datasets
            # Ensure the correct datasets are passed.
            ds_wetconv = opened_datasets.get('WetLossConv')
            ds_wetls = opened_datasets.get('WetLossLS')
            ds_drydep = opened_datasets.get('DryDep')

            if ds_wetconv is None:
                print("Error: 'WetLossConv' dataset not found for total deposition calculation.")
                print("Please ensure you have files like 'GEOSChem.WetLossConv....nc4'.")
                continue
            if ds_wetls is None:
                print("Error: 'WetLossLS' dataset not found for total deposition calculation.")
                print("Please ensure you have files like 'GEOSChem.WetLossLS....nc4'.")
                continue
            if ds_drydep is None:
                print("Error: 'DryDep' dataset not found for total deposition calculation.")
                print("Please ensure you have files like 'GEOSChem.DryDep....nc4'.")
                continue

            # Call the analysis function, passing the opened xarray datasets
            wet_dep_result, dry_dep_result = quantify_total_deposition(
                ds_wetconv, ds_wetls, ds_drydep
            )
            analysis_performed = True
            break # Exit menu loop after analysis
        elif choice == '2':
            print(f"\nInitiating: {menu_options['2']}...")
            # We need to pass the specific datasets required by quantify_hg2_to_ssa_transfer
            ds_mercurychem = opened_datasets.get('MercuryChem')
            ds_statemet = opened_datasets.get('StateMet')

            if ds_mercurychem is None:
                print("Error: 'MercuryChem' dataset not found for Hg(II) to SSA transfer calculation.")
                print("Please ensure you have files like 'GEOSChem.MercuryChem....nc4'.")
                continue
            if ds_statemet is None:
                print("Error: 'StateMet' dataset not found for Hg(II) to SSA transfer calculation.")
                print("Please ensure you have files like 'GEOSChem.StateMet....nc4'.")
                continue

            # Call the analysis function
            hg2_ssa_transfer_result = quantify_hg2_to_ssa_transfer(
                ds_mercurychem, ds_statemet
            )
            analysis_performed = True
            break # Exit menu loop after analysis
        elif choice == '0':
            print("Exiting the GEOS-Chem Mercury Analysis Tool. Goodbye!")
            analysis_performed = False # Indicate no analysis was performed
            break
        else:
            print("Invalid choice. Please enter a number corresponding to an option.")

    # --- Step 4: Display Results (Textual Summary) ---
    if analysis_performed:
        print("\n--- Step 4: Displaying Analysis Results ---")
        if wet_dep_result is not None or dry_dep_result is not None:
            display_mercury_results(wet_dep_result=wet_dep_result, dry_dep_result=dry_dep_result)
        elif hg2_ssa_transfer_result is not None:
            display_mercury_results(hg2_ssa_transfer_result=hg2_ssa_transfer_result)
        else:
            print("No results to display. Analysis may have failed or produced no valid data.")

        # --- Step 5: Offer Plotting ---
        print("\n--- Step 5: Plotting Results (Optional) ---")
        plot_choice = input("Do you want to plot the results? (yes/no): ").strip().lower()
        if plot_choice == 'yes':
            print("Initiating plotting process...")
            # Pass the results of the analysis to the plotting function
            offer_and_execute_plotting(
                wet_dep_result=wet_dep_result,
                dry_dep_result=dry_dep_result,
                hg2_ssa_transfer_result=hg2_ssa_transfer_result
            )
        else:
            print("Skipping plotting. Analysis complete.")

    # --- Final Cleanup: Close all opened datasets ---
    print("\n--- Final Step: Cleaning Up Resources ---")
    print("Closing all opened xarray datasets to free up memory.")
    for desc, ds in opened_datasets.items():
        if ds is not None: # Ensure the dataset was actually opened
            ds.close()
            print(f"  - Dataset '{desc}' closed.")
    print("Resource cleanup complete. Thank you for using the tool!")


if __name__ == "__main__":
    run_mercury_analysis_tool()
