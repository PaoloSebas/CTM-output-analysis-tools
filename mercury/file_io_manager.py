# file_io_manager.py
import os
import re
import xarray as xr
from tabulate import tabulate

def find_geoschem_nc4_files(folder_path: str) -> list[dict]:
    """
    Scans a given folder path for GEOS-Chem .nc4 files and extracts
    their names and a basic description. This version is more flexible
    regarding where 'GEOSChem.' appears in the filename.

    Args:
        folder_path (str): The path to the directory containing .nc4 files.

    Returns:
        list[dict]: A list of dictionaries, where each dictionary contains
                    'filename' and 'description' for each found .nc4 file.
                    Returns an empty list if the folder is not found or no files.
    """
    if not os.path.isdir(folder_path):
        print(f"Error: Folder not found at '{folder_path}'.")
        return []

    nc4_files_info = []
    # Updated Regex:
    # Looks for 'GEOSChem.' anywhere in the filename, then captures the next alphanumeric part.
    # `.*?` is a non-greedy match for any characters before 'GEOSChem.'.
    # `GEOSChem\.` matches the literal string "GEOSChem."
    # `([a-zA-Z0-9]+)` is the capturing group for the description (e.g., DryDep, MercuryChem).
    # `\..*\.nc4$` matches the dot after the description, any subsequent characters,
    # and ensures it ends with '.nc4'.
    geoschem_pattern = re.compile(r'.*?GEOSChem\.([a-zA-Z0-9]+)\..*\.nc4$')

    for filename in os.listdir(folder_path):
        if filename.endswith('.nc4'): # First, filter by the .nc4 extension
            match = geoschem_pattern.search(filename)
            description = "Unknown"
            if match:
                # If a match is found by the regex, it's considered a GEOS-Chem related file
                description = match.group(1) # Extract the captured description part
                nc4_files_info.append({'filename': filename, 'description': description})
    return nc4_files_info

def display_file_summary_table(files_info: list[dict]):
    """
    Displays a formatted table of file names and descriptions.

    Args:
        files_info (list[dict]): A list of dictionaries, each with 'filename' and 'description'.
    """
    if not files_info:
        print("No .nc4 files found or recognized in the specified folder.")
        return

    headers = ["ID", "Filename", "Description"]
    table_data = []
    for i, file_info in enumerate(files_info):
        table_data.append([i + 1, file_info['filename'], file_info['description']])

    print("\n--- Folder Content Summary ---")
    print(tabulate(table_data, headers=headers, tablefmt="grid"))
    print("------------------------------")

def open_geoschem_datasets(folder_path: str) -> dict[str, xr.Dataset]:
    """
    Opens all recognized GEOS-Chem .nc4 files in a folder using xarray.
    It returns a dictionary where keys are simplified descriptions
    (e.g., 'DryDep', 'MercuryChem') and values are the opened xarray Datasets.
    If multiple files of the same description exist (e.g., time series),
    they are opened as a single dataset using open_mfdataset.

    Args:
        folder_path (str): The path to the directory containing .nc4 files.

    Returns:
        dict[str, xr.Dataset]: A dictionary of opened xarray Datasets.
                               Returns an empty dictionary if no files are opened.
    """
    files_info = find_geoschem_nc4_files(folder_path)
    if not files_info:
        print("No GEOS-Chem .nc4 files found to open.")
        return {}

    # Group files by their description
    grouped_files = {}
    for file_info in files_info:
        desc = file_info['description']
        full_path = os.path.join(folder_path, file_info['filename'])
        if desc not in grouped_files:
            grouped_files[desc] = []
        grouped_files[desc].append(full_path)

    opened_datasets = {}
    print("\nAttempting to open files with Xarray...")
    for desc, paths in grouped_files.items():
        try:
            # Use open_mfdataset for potentially multiple files per description
            # combine="by_coords" is often good for time series where files are split by time
            ds = xr.open_mfdataset(paths, combine="by_coords", engine="netcdf4")
            opened_datasets[desc] = ds
            print(f"Successfully opened '{desc}' dataset(s).")
        except FileNotFoundError as e:
            print(f"Error: File not found for '{desc}'. {e}")
        except Exception as e:
            print(f"Error opening dataset for '{desc}': {e}")

    return opened_datasets


if __name__ == "__main__":
    # --- Example Usage for Testing ---
    print("--- Testing file_io_manager.py (Updated Logic) ---")

    # Create a dummy folder and some dummy .nc4 files for testing
    test_folder = "test_geoschem_data_flexible"
    os.makedirs(test_folder, exist_ok=True)

    # Create empty dummy files with a prefix before "GEOSChem."
    dummy_files = [
        "MyRun_ANU_GEOSChem.DryDep.20190101_0000z.nc4",
        "MyRun_ANU_GEOSChem.DryDep.20190102_0000z.nc4", # Simulate multiple files for mfdataset
        "Simulation1_GEOSChem.MercuryChem.20190101_0000z.nc4",
        "Base_GEOSChem.StateMet.20190101_0000z.nc4",
        "Experiment_GEOSChem.WetLossConv.20190101_0000z.nc4",
        "Just_A_Test_GEOSChem.SomeOtherVar.20190101_0000z.nc4",
        "non_geoschem_file.txt", # Should be ignored
        "GEOSChem.DirectStart.20190101_0000z.nc4" # Still works for files starting directly
    ]
    for df in dummy_files:
        with open(os.path.join(test_folder, df), 'w') as f:
            f.write("dummy content") # Real nc4 files are binary, this is just for file existence

    print(f"Created dummy files in: {test_folder}\n")

    # Test find_geoschem_nc4_files
    found_files = find_geoschem_nc4_files(test_folder)
    print("Files found by find_geoschem_nc4_files:")
    for f_info in found_files:
        print(f"  - {f_info['filename']} (Description: {f_info['description']})")

    # Test display_file_summary_table
    display_file_summary_table(found_files)

    # Test open_geoschem_datasets (will likely fail for empty files, but shows logic)
    print("\nAttempting to open datasets (expect errors for dummy files as they are not valid NetCDF):")
    opened_datasets = open_geoschem_datasets(test_folder)
    print(f"Number of datasets supposedly opened: {len(opened_datasets)}")

    # Clean up dummy files and folder
    for df in dummy_files:
        os.remove(os.path.join(test_folder, df))
    os.rmdir(test_folder)
    print(f"\nCleaned up '{test_folder}' folder and dummy files.")
    print("\n--- Testing Complete ---")
