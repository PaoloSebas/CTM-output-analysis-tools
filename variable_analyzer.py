# variable_analyzer.py
import xarray as xr
from tabulate import tabulate
import csv # Import the csv module
from datetime import datetime # For unique filenames
import os # <--- ADD THIS LINE


def analyze_and_list_variables(opened_datasets: dict[str, xr.Dataset]) -> list[dict]:
    """
    Analyzes a dictionary of opened xarray Datasets to extract all data variables
    along with their originating dataset (file description) and other metadata
    like dimensions, units, and long_name (description).

    Args:
        opened_datasets (dict[str, xr.Dataset]): A dictionary where keys are
                                                 dataset descriptions (e.g., 'DryDep')
                                                 and values are xarray.Dataset objects.

    Returns:
        list[dict]: A list of dictionaries, where each dictionary contains:
                    - 'id': A sequential ID number.
                    - 'variable_name': The name of the variable.
                    - 'source_dataset': The description of the dataset it came from.
                    - 'dimensions': A tuple of the variable's dimensions.
                    - 'units': The units of the variable (if available).
                    - 'long_name': The long description of the variable (if available).
    """
    if not opened_datasets:
        print("No datasets were provided for variable analysis.")
        return []

    variable_list = []
    current_id = 1

    for dataset_desc, ds in opened_datasets.items():
        for var_name, var_data in ds.data_vars.items():
            units = var_data.attrs.get('units', 'N/A')
            long_name = var_data.attrs.get('long_name', 'No description available')
            variable_info = {
                'id': current_id,
                'variable_name': var_name,
                'source_dataset': dataset_desc,
                'dimensions': tuple(var_data.dims),
                'units': units,
                'long_name': long_name
            }
            variable_list.append(variable_info)
            current_id += 1
    return variable_list

def display_variable_summary_table(variable_list: list[dict]):
    """
    Displays a formatted table of variables, their source datasets, dimensions,
    units, and long description.

    Args:
        variable_list (list[dict]): A list of dictionaries, as returned by analyze_and_list_variables.
    """
    if not variable_list:
        print("No variables found or processed for display.")
        return

    headers = ["ID", "Variable Name", "Source File Type", "Dimensions", "Units", "Description"]
    table_data = []
    for var_info in variable_list:
        table_data.append([
            var_info['id'],
            var_info['variable_name'],
            var_info['source_dataset'],
            ", ".join(var_info['dimensions']),
            var_info['units'],
            var_info['long_name']
        ])

    print("\n--- Available Variables Summary ---")
    print(tabulate(table_data, headers=headers, tablefmt="grid"))
    print("-----------------------------------")


def export_variables_to_csv(variable_list: list[dict], output_folder: str = ".") -> str or None:
    """
    Exports the list of variables to a CSV file.

    Args:
        variable_list (list[dict]): A list of dictionaries, as returned by analyze_and_list_variables.
        output_folder (str): The folder where the CSV file should be saved.
                             Defaults to the current working directory.

    Returns:
        str or None: The absolute path to the created CSV file if successful,
                     None otherwise.
    """
    if not variable_list:
        print("No variables to export to CSV.")
        return None

    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Generate a unique filename with a timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_filename = os.path.join(output_folder, f"geoschem_variables_summary_{timestamp}.csv")

    try:
        # Determine headers from the keys of the first dictionary
        # Assuming all dictionaries have the same keys
        fieldnames = list(variable_list[0].keys())

        with open(csv_filename, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader() # Write the header row
            writer.writerows(variable_list) # Write all data rows

        print(f"Variable summary successfully exported to: {os.path.abspath(csv_filename)}")
        return os.path.abspath(csv_filename)
    except IOError as e:
        print(f"Error writing CSV file '{csv_filename}': {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred during CSV export: {e}")
        return None


if __name__ == "__main__":
    # --- Example Usage for Testing (Updated with export) ---
    print("--- Testing variable_analyzer.py (with CSV export) ---")

    # Create dummy xarray Datasets (same as previous test)
    ds_drydep = xr.Dataset(
        {
            "DryDep_Hg0": (("time", "lat", "lon"), [[[1, 2], [3, 4]]]),
            "DryDep_Hg2": (("time", "lat", "lon"), [[[5, 6], [7, 8]]]),
            "DryDep_ALL": (("time", "lat", "lon"), [[[9, 10], [11, 12]]]),
        },
        coords={
            "time": [0], "lat": [10, 20], "lon": [30, 40]
        },
        attrs={'Description': 'DryDep data'}
    )
    ds_drydep["DryDep_Hg0"].attrs['units'] = 'kg/s'
    ds_drydep["DryDep_Hg0"].attrs['long_name'] = 'Dry Deposition Flux of Elemental Mercury'
    ds_drydep["DryDep_Hg2"].attrs['units'] = 'kg/s'
    ds_drydep["DryDep_Hg2"].attrs['long_name'] = 'Dry Deposition Flux of Gaseous Oxidized Mercury'

    ds_mercurychem = xr.Dataset(
        {
            "Hg2GasToSSA": (("time", "lev", "lat", "lon"), [[[[1.1]]], [[[2.2]]]]),
            "Conc_Hg0": (("time", "lev", "lat", "lon"), [[[[3.3]]], [[[4.4]]]]),
        },
        coords={
            "time": [0, 1], "lev": [1,2], "lat": [50], "lon": [60]
        },
        attrs={'Description': 'MercuryChem data'}
    )
    ds_mercurychem["Hg2GasToSSA"].attrs['units'] = 'molec cm-3 s-1'
    ds_mercurychem["Hg2GasToSSA"].attrs['long_name'] = 'Conversion Rate of Hg2 to Soluble Species'
    ds_mercurychem["Conc_Hg0"].attrs['units'] = 'ppbv'
    ds_mercurychem["Conc_Hg0"].attrs['long_name'] = 'Concentration of Gaseous Elemental Mercury'

    dummy_opened_datasets = {
        "DryDep": ds_drydep,
        "MercuryChem": ds_mercurychem
    }

    analyzed_vars = analyze_and_list_variables(dummy_opened_datasets)
    display_variable_summary_table(analyzed_vars)

    # Test the new export function
    print("\nAttempting to export variables to CSV:")
    csv_path = export_variables_to_csv(analyzed_vars)
    if csv_path:
        print(f"Check your current directory for: {os.path.basename(csv_path)}")

    print("\n--- Testing Complete ---")
