import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# --- Mercury Analysis Functions (from mercury_analysis.py) ---
MOLAR_MASS_HG_G_PER_MOL = 200.59
MOLAR_MASS_HG_KG_PER_MOL = MOLAR_MASS_HG_G_PER_MOL / 1000
AVOGADRO_NUMBER = 6.022e23
CM2_TO_M2 = 1e-4
CM3_TO_M3 = 1e-6

def quantify_total_deposition(wetlossconv_file, wetlossls_file, drydep_file):
    ds_wetconv = ds_wetls = ds_drydep = None
    total_wet_deposition = total_dry_deposition = None
    try:
        ds_wetconv = xr.open_dataset(wetlossconv_file)
        ds_wetls = xr.open_dataset(wetlossls_file)
        ds_drydep = xr.open_dataset(drydep_file)
        dryd_vars_to_extract = [
            'DryDep_Hg0', 'DryDep_HgCl2', 'DryDep_HgOHOH', 'DryDep_HgClOH',
            'DryDep_HgBr2', 'DryDep_HgBrOH', 'DryDep_Hg2ORGP', 'DryDep_Hg2ClP'
        ]
        extracted_dry_data_molec_cm2_s = [
            ds_drydep[var] for var in dryd_vars_to_extract if var in ds_drydep.data_vars
        ]
        if extracted_dry_data_molec_cm2_s and 'AREA' in ds_drydep.data_vars:
            grid_area_cm2 = ds_drydep['AREA'] / CM2_TO_M2
            dry_deposition_kg_s_per_species = [
                da * (MOLAR_MASS_HG_KG_PER_MOL / AVOGADRO_NUMBER) * grid_area_cm2
                for da in extracted_dry_data_molec_cm2_s
            ]
            total_dry_deposition = sum(dry_deposition_kg_s_per_species)
            total_dry_deposition.attrs['units'] = 'kg/s'
        wetd_conv_vars_to_extract = [
            'WetLossConv_Hg2ORGP', 'WetLossConv_Hg2ClP', 'WetLossConv_HgCl2',
            'WetLossConv_HgOHOH', 'WetLossConv_HgClOH', 'WetLossConv_HgBr2', 'WetLossConv_HgBrOH'
        ]
        wetd_ls_vars_to_extract = [
            'WetLossLS_Hg2ORGP', 'WetLossLS_Hg2ClP', 'WetLossLS_HgCl2',
            'WetLossLS_HgOHOH', 'WetLossLS_HgClOH', 'WetLossLS_HgBr2', 'WetLossLS_HgBrOH'
        ]
        extracted_wet_data_kg_s = []
        for var in wetd_conv_vars_to_extract:
            if var in ds_wetconv.data_vars:
                extracted_wet_data_kg_s.append(ds_wetconv[var])
        for var in wetd_ls_vars_to_extract:
            if var in ds_wetls.data_vars:
                extracted_wet_data_kg_s.append(ds_wetls[var])
        if extracted_wet_data_kg_s:
            total_wet_deposition = sum(extracted_wet_data_kg_s)
            total_wet_deposition.attrs['units'] = 'kg/s'
    finally:
        for ds in (ds_wetconv, ds_wetls, ds_drydep):
            if ds is not None:
                ds.close()
    return total_wet_deposition, total_dry_deposition

def quantify_hg2_to_ssa_transfer(mercurychem_file, statemet_file):
    ds_mercurychem = ds_statemet = None
    hg2_to_ssa_transfer_rate_kg_s = None
    try:
        ds_mercurychem = xr.open_dataset(mercurychem_file)
        ds_statemet = xr.open_dataset(statemet_file)
        if 'Hg2GasToSSA' in ds_mercurychem.data_vars and 'Met_AIRVOL' in ds_statemet.data_vars:
            transfer_rate_molec_cm3_s = ds_mercurychem['Hg2GasToSSA']
            grid_volume_cm3 = ds_statemet['Met_AIRVOL'] / CM3_TO_M3
            hg2_to_ssa_transfer_rate_kg_s = transfer_rate_molec_cm3_s * \
                (MOLAR_MASS_HG_KG_PER_MOL / AVOGADRO_NUMBER) * grid_volume_cm3
            hg2_to_ssa_transfer_rate_kg_s.attrs['units'] = 'kg/s'
    finally:
        for ds in (ds_mercurychem, ds_statemet):
            if ds is not None:
                ds.close()
    return hg2_to_ssa_transfer_rate_kg_s

# --- Minimal File I/O and Variable Analysis ---
def find_geoschem_nc4_files(folder_path):
    files = []
    for fname in os.listdir(folder_path):
        if fname.endswith('.nc4'):
            # Very basic guessing of file type
            if 'WetLossConv' in fname:
                desc = 'WetLossConv'
            elif 'WetLossLS' in fname:
                desc = 'WetLossLS'
            elif 'DryDep' in fname:
                desc = 'DryDep'
            elif 'MercuryChem' in fname:
                desc = 'MercuryChem'
            elif 'StateMet' in fname:
                desc = 'StateMet'
            else:
                desc = 'Other'
            files.append({'filename': fname, 'description': desc})
    return files

def display_file_summary_table(files_info):
    print("\nAvailable .nc4 files:")
    print("{:<3} {:<30} {:<15}".format("No", "Filename", "Description"))
    for idx, file_info in enumerate(files_info, 1):
        print("{:<3} {:<30} {:<15}".format(idx, file_info['filename'], file_info['description']))

def open_geoschem_datasets(folder_path):
    datasets = {}
    for file_info in find_geoschem_nc4_files(folder_path):
        path = os.path.join(folder_path, file_info['filename'])
        try:
            ds = xr.open_dataset(path)
            datasets[file_info['description']] = ds
        except Exception as e:
            print(f"Failed to open {path}: {e}")
    return datasets

def analyze_and_list_variables(opened_datasets):
    all_var_info = []
    for ds_key, ds in opened_datasets.items():
        for var in ds.data_vars:
            units = ds[var].attrs.get('units', 'unknown')
            shape = ds[var].shape
            all_var_info.append({'dataset': ds_key, 'variable': var, 'units': units, 'shape': shape})
    return all_var_info

def display_variable_summary_table(all_variables_info):
    print("\nVariables in loaded datasets:")
    print("{:<15} {:<30} {:<12} {}".format("Dataset", "Variable", "Units", "Shape"))
    for info in all_variables_info:
        print("{:<15} {:<30} {:<12} {}".format(
            info['dataset'], info['variable'], info['units'], info['shape']))

def export_variables_to_csv(all_variables_info):
    import csv
    outname = "variable_summary.csv"
    with open(outname, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['dataset', 'variable', 'units', 'shape'])
        writer.writeheader()
        for row in all_variables_info:
            writer.writerow(row)
    print(f"Variable summary exported to {outname}")

# --- Simple Plotting ---
def plot_variable_global(data, var_name):
    plt.figure(figsize=(8, 4))
    if 'lat' in data.dims and 'lon' in data.dims:
        plt.pcolormesh(data['lon'], data['lat'], np.squeeze(data), shading='auto')
        plt.title(f"Global Distribution of {var_name}")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.colorbar(label=data.attrs.get('units', ''))
        plt.show()
    else:
        print("Cannot plot: variable lacks lat/lon dimensions.")

def plot_variable_antarctic(data, var_name):
    if 'lat' in data.dims and 'lon' in data.dims:
        antarctic_mask = data['lat'] < -60
        plt.figure(figsize=(6, 4))
        plt.pcolormesh(data['lon'], data['lat'][antarctic_mask], np.squeeze(data.sel(lat=data['lat'][antarctic_mask])), shading='auto')
        plt.title(f"Antarctic Distribution of {var_name}")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.colorbar(label=data.attrs.get('units', ''))
        plt.show()
    else:
        print("Cannot plot: variable lacks lat/lon dimensions.")

# --- Main Menu and Workflow ---
def display_main_menu():
    print("\n--- GEOS-Chem Data Analysis Tool ---")
    print("1. Load & Process GEOS-Chem Data Folder")
    print("2. Perform Mercury Analysis")
    print("3. Plotting a selected variable")
    print("4. Exit")
    print("------------------------------------")

def process_data_folder():
    print("\n--- Load & Process GEOS-Chem Data Folder ---")
    folder_path = input("Please provide the path that contains the nc4 files to be processed: ").strip()
    if not folder_path or not os.path.isdir(folder_path):
        print("Invalid path.")
        return None, None, None
    files_info = find_geoschem_nc4_files(folder_path)
    if not files_info:
        print("No .nc4 files found.")
        return None, None, None
    display_file_summary_table(files_info)
    opened_datasets = open_geoschem_datasets(folder_path)
    if not opened_datasets:
        print("No datasets could be opened.")
        return None, None, None
    all_variables_info = analyze_and_list_variables(opened_datasets)
    display_variable_summary_table(all_variables_info)
    if input("Export variable list to CSV? (yes/no): ").strip().lower() == 'yes':
        export_variables_to_csv(all_variables_info)
    return opened_datasets, files_info, folder_path

def perform_mercury_analysis(opened_datasets, files_info, folder_path):
    print("\n--- Perform Mercury Analysis ---")
    if opened_datasets is None or not opened_datasets:
        print("Load data folder first.")
        return
    file_paths_by_description = {f['description']: os.path.join(folder_path, f['filename']) for f in files_info}
    try:
        wetlossconv = file_paths_by_description['WetLossConv']
        wetlossls = file_paths_by_description['WetLossLS']
        drydep = file_paths_by_description['DryDep']
        mercurychem = file_paths_by_description['MercuryChem']
        statemet = file_paths_by_description['StateMet']
    except KeyError as e:
        print(f"Required file missing: {e}")
        return
    twd, tdd = quantify_total_deposition(wetlossconv, wetlossls, drydep)
    print("Wet Deposition, first timestep sum (kg/s):", twd.isel(time=0).sum().item() if twd is not None and 'time' in twd.dims else "N/A")
    print("Dry Deposition, first timestep sum (kg/s):", tdd.isel(time=0).sum().item() if tdd is not None and 'time' in tdd.dims else "N/A")
    hg2ssa = quantify_hg2_to_ssa_transfer(mercurychem, statemet)
    if hg2ssa is not None:
        print("Hg2 to SSA Transfer, first timestep sum (kg/s):", hg2ssa.isel(time=0).sum().item() if 'time' in hg2ssa.dims else hg2ssa.sum().item())
    else:
        print("Hg2 to SSA Transfer calculation failed.")

def plot_selected_variable(opened_datasets):
    print("\n--- Plotting a selected variable ---")
    if opened_datasets is None or not opened_datasets:
        print("Load data folder first.")
        return
    available_vars_map = {}
    display_idx = 1
    for ds_key, ds_obj in opened_datasets.items():
        for var_name in ds_obj.data_vars:
            available_vars_map[display_idx] = (ds_key, var_name, ds_obj[var_name])
            units = ds_obj[var_name].attrs.get('units', 'unitless')
            print(f"{display_idx}. {var_name} (from {ds_key}, units: {units})")
            display_idx += 1
    if not available_vars_map:
        print("No variables to plot.")
        return
    try:
        var_choice = int(input("Select variable number to plot (0 to cancel): ").strip())
        if var_choice == 0:
            return
        _ds_key, selected_var_name, selected_data_array = available_vars_map[var_choice]
        data_to_plot_slice = selected_data_array
        if 'time' in data_to_plot_slice.dims and data_to_plot_slice['time'].size > 1:
            idx = int(input(f"Variable has {data_to_plot_slice['time'].size} time steps. Pick time index (0-based): ").strip())
            data_to_plot_slice = data_to_plot_slice.isel(time=idx)
        if 'lev' in data_to_plot_slice.dims and data_to_plot_slice['lev'].size > 1:
            idx = int(input(f"Variable has {data_to_plot_slice['lev'].size} levels. Pick level index (0-based): ").strip())
            data_to_plot_slice = data_to_plot_slice.isel(lev=idx)
        print("Plot options: 1. Global  2. Antarctic")
        plt_choice = input("Choose plot type (1 or 2): ").strip()
        if plt_choice == '1':
            plot_variable_global(data_to_plot_slice, selected_var_name)
        elif plt_choice == '2':
            plot_variable_antarctic(data_to_plot_slice, selected_var_name)
    except Exception as e:
        print(f"Plotting failed: {e}")

def main():
    print("Welcome to the GEOS-Chem Data Analysis Tool!")
    current_opened_datasets = None
    current_files_info = None
    current_folder_path = None
    while True:
        display_main_menu()
        choice = input("Enter your choice (1-4): ").strip()
        if choice == '1':
            current_opened_datasets, current_files_info, current_folder_path = process_data_folder()
        elif choice == '2':
            perform_mercury_analysis(current_opened_datasets, current_files_info, current_folder_path)
        elif choice == '3':
            plot_selected_variable(current_opened_datasets)
        elif choice == '4':
            print("Exiting. Closing datasets...")
            if current_opened_datasets:
                for ds in current_opened_datasets.values():
                    ds.close()
            break
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()
