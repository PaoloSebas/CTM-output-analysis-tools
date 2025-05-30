# mercury_analysis.py
import xarray as xr
import os
import numpy as np # Often useful for numerical operations

# --- Define Universal Constants ---
# These are fundamental physical constants used in atmospheric chemistry calculations.
# Displaying them explicitly helps the user understand the numerical values involved
# in the unit conversions and scientific calculations.
MOLAR_MASS_HG_G_PER_MOL = 200.59 # grams per mole (Molar Mass of Mercury)
MOLAR_MASS_HG_KG_PER_MOL = MOLAR_MASS_HG_G_PER_MOL / 1000 # kilograms per mole

AVOGADRO_NUMBER = 6.022e23 # molecules per mole (or atoms per mole for elemental species)

# Conversion factors for spatial units
CM2_TO_M2 = 1e-4 # Conversion factor: 1 cm^2 = 1e-4 m^2 (since 1 m = 100 cm, 1 m^2 = 100^2 cm^2 = 10000 cm^2)
CM3_TO_M3 = 1e-6 # Conversion factor: 1 cm^3 = 1e-6 m^3 (since 1 m = 100 cm, 1 m^3 = 100^3 cm^3 = 1000000 cm^3)

# Print out constants for didactic purposes at the module level
print("--- Universal Constants Used in Calculations ---")
print(f"  Molar Mass of Mercury (Hg): {MOLAR_MASS_HG_G_PER_MOL} g/mol ({MOLAR_MASS_HG_KG_PER_MOL} kg/mol)")
print(f"  Avogadro's Number (N_A): {AVOGADRO_NUMBER:.3e} molecules/mol")
print(f"  Unit Conversion: 1 cm^2 = {CM2_TO_M2} m^2")
print(f"  Unit Conversion: 1 cm^3 = {CM3_TO_M3} m^3")
print("----------------------------------------------\n")

def quantify_total_deposition(wetlossconv_file: str, wetlossls_file: str, drydep_file: str):
    """
    Quantifies the total wet and dry deposition of mercury species from GEOS-Chem output.

    This function aims to be didactic, explaining the scientific principles,
    mathematical steps, and unit conversions involved in calculating total
    mercury deposition.

    Args:
        wetlossconv_file (str): Path to the GEOS-Chem WetLossConv_####.nc4 file,
                                containing convective wet deposition fluxes.
        wetlossls_file (str): Path to the GEOS-Chem WetLossLS_####.nc4 file,
                              containing large-scale wet deposition fluxes.
        drydep_file (str): Path to the GEOS-Chem DryDep_####.nc4 file,
                           containing dry deposition fluxes and surface area data.

    Returns:
        tuple[xr.DataArray, xr.DataArray]: A tuple containing:
            - total_wet_deposition (xr.DataArray): Total wet deposition flux in kg/s.
            - total_dry_deposition (xr.DataArray): Total dry deposition flux in kg/s.
            Returns (None, None) if files cannot be opened or critical variables are missing.
    """
    print("\n" + "="*80)
    print("--- QUANTIFYING TOTAL MERCURY DEPOSITION ---".center(80))
    print("="*80)

    ds_wetconv, ds_wetls, ds_drydep = None, None, None # Initialize datasets to None

    try:
        # --- Step 1: Data Acquisition - Open Required Datasets ---
        print("\n[STEP 1/6] Data Acquisition: Opening GEOS-Chem Output Files...")
        print("------------------------------------------------------------------")
        print("Accessing model output files is the first crucial step. We use xarray\n"
              "to open NetCDF files, which are the standard format for atmospheric\n"
              "model data. This allows us to load variables and their associated metadata\n"
              "(such as units, dimensions, and long names).\n")

        # Open WetLossConv file
        print(f"  - Attempting to open Convective Wet Loss file: {os.path.basename(wetlossconv_file)}")
        try:
            ds_wetconv = xr.open_dataset(wetlossconv_file)
            print(f"    SUCCESS: Loaded WetLossConv dataset. Available variables: {list(ds_wetconv.data_vars)}")
        except Exception as e:
            print(f"    ERROR: Could not open WetLossConv file '{wetlossconv_file}'. Please check path and permissions: {e}")
            return None, None

        # Open WetLossLS file
        print(f"  - Attempting to open Large-Scale Wet Loss file: {os.path.basename(wetlossls_file)}")
        try:
            ds_wetls = xr.open_dataset(wetlossls_file)
            print(f"    SUCCESS: Loaded WetLossLS dataset. Available variables: {list(ds_wetls.data_vars)}")
        except Exception as e:
            print(f"    ERROR: Could not open WetLossLS file '{wetlossls_file}'. Please check path and permissions: {e}")
            return None, None

        # Open DryDep file
        print(f"  - Attempting to open Dry Deposition file: {os.path.basename(drydep_file)}")
        try:
            ds_drydep = xr.open_dataset(drydep_file)
            print(f"    SUCCESS: Loaded DryDep dataset. Available variables: {list(ds_drydep.data_vars)}")
        except Exception as e:
            print(f"    ERROR: Could not open DryDep file '{drydep_file}'. Please check path and permissions: {e}")
            return None, None
        
        print("\n[STEP 1 Status]: All required datasets successfully opened.")

        # --- Step 2: Scientific Context & Variable Identification for Dry Deposition ---
        print("\n[STEP 2/6] Dry Deposition: Scientific Context & Variable Identification")
        print("------------------------------------------------------------------------")
        print("Dry deposition is the direct transfer of gases and aerosols from the\n"
              "atmosphere to the Earth's surface in the absence of precipitation.\n"
              "It's a continuous process governed by turbulence and surface properties.\n"
              "In GEOS-Chem, dry deposition fluxes are often provided in units of\n"
              "'molecules per square centimeter per second' (molec cm^-2 s^-1).\n")
        print("To calculate total dry deposition in 'kilograms per second' (kg/s) for the\n"
              "entire model domain (or a specific grid cell), we need to sum the individual\n"
              "species and convert their units from a number flux to a mass flux.\n")
        print("The conversion formula from 'molecules/cm^2/s' to 'kg/s' for a specific grid cell is:\n")
        print("  $$ \\text{Flux}_{\\text{dry, kg/s}} = \\text{Flux}_{\\text{dry, molec/cm}^2\\text{/s}} \\times \\left( \\frac{\\text{M}_{\\text{Hg}}}{\\text{N}_{\\text{A}}} \\right) \\times \\text{Area}_{\\text{grid, cm}^2} $$")
        print("Where:\n"
              f"  - $M_{{Hg}}$: Molar Mass of Mercury ({MOLAR_MASS_HG_KG_PER_MOL} kg/mol)\n"
              f"  - $N_{{A}}$: Avogadro's Number ({AVOGADRO_NUMBER:.3e} molec/mol)\n"
              "  - $\\text{Area}_{\\text{grid, cm}^2}$: Surface area of the grid cell in cm$^2$ (obtained from 'AREA' variable).\n")
        print("Identifying relevant dry deposition species (prefixed with 'DryDep_') from the provided list:")

        # List of dry deposition species to include in total. These are from your provided CSV.
        dryd_vars_to_extract = [
            'DryDep_Hg0',      # Gaseous Elemental Mercury
            'DryDep_HgCl2',    # Gaseous Divalent Mercury (e.g., HgCl2)
            'DryDep_HgOHOH',   # Gaseous Divalent Mercury (e.g., Hg(OH)2)
            'DryDep_HgClOH',   # Gaseous Divalent Mercury (e.g., HgClOH)
            'DryDep_HgBr2',    # Gaseous Divalent Mercury (e.g., HgBr2)
            'DryDep_HgBrOH',   # Gaseous Divalent Mercury (e.g., HgBrOH)
            'DryDep_Hg2ORGP',  # Particulate Organic Divalent Mercury
            'DryDep_Hg2ClP'    # Particulate Inorganic Divalent Mercury (e.g., HgCl2(p))
            # You can add or remove species here based on your specific definition of "total dry deposition".
        ]
        extracted_dry_data_molec_cm2_s = [] # To store DataArrays in their original units (molec/cm^2/s)
        
        for var in dryd_vars_to_extract:
            if var in ds_drydep.data_vars:
                extracted_dry_data_molec_cm2_s.append(ds_drydep[var])
                print(f"  - Found and extracted '{var}' (Shape: {ds_drydep[var].shape}, Units: {ds_drydep[var].attrs.get('units', 'N/A')})")
            else:
                print(f"  - WARNING: Variable '{var}' not found in DryDep dataset. This species will be skipped for the dry deposition sum.")
        
        if not extracted_dry_data_molec_cm2_s:
            print("  - ERROR: No relevant dry deposition variables were found based on the provided list. Cannot calculate total dry deposition.")
            total_dry_deposition = None
        else:
            # Get the grid cell area from the DryDep dataset for unit conversion
            if 'AREA' in ds_drydep.data_vars:
                grid_area_m2 = ds_drydep['AREA'] # The 'AREA' variable is typically in m^2 in GEOS-Chem output.
                print(f"\n  - Extracted 'AREA' from DryDep dataset (Shape: {grid_area_m2.shape}, Units: {grid_area_m2.attrs.get('units', 'N/A')}).")
                # Convert AREA from m^2 to cm^2 to match the units of the flux (molec/cm^2/s).
                # This ensures consistent unit cancellation in the conversion formula.
                grid_area_cm2 = grid_area_m2 / CM2_TO_M2 # Equivalent to grid_area_m2 * 1e4
                print(f"    -> Converted 'AREA' from m^2 to cm^2 for unit consistency in the calculation.")
            else:
                print("  - CRITICAL ERROR: 'AREA' variable (grid cell surface area) not found in DryDep dataset. Cannot convert dry deposition flux from molec/cm^2/s to kg/s.")
                return None, None

            # --- Step 3: Calculation of Total Dry Deposition (in kg/s) ---
            print("\n[STEP 3/6] Dry Deposition: Performing Unit Conversion and Summation")
            print("-------------------------------------------------------------------")
            print("Each dry deposition flux (in molec cm^-2 s^-1) is converted to a mass\n"
                  "flux in 'kg per second' for each grid cell. This involves multiplying\n"
                  "by the molar mass of mercury ($M_{{Hg}}$), dividing by Avogadro's number ($N_{{A}}$),\n"
                  "and multiplying by the surface area of the grid cell ($\\text{Area}_{\\text{grid, cm}^2}$). \n"
                  "Finally, these mass fluxes for all individual mercury species are summed.\n")

            dry_deposition_kg_s_per_species = []
            for da in extracted_dry_data_molec_cm2_s:
                # Apply the conversion formula:
                # (molec/cm^2/s) * (kg/mol) / (molec/mol) * (cm^2) = kg/s
                converted_da = da * (MOLAR_MASS_HG_KG_PER_MOL / AVOGADRO_NUMBER) * grid_area_cm2
                dry_deposition_kg_s_per_species.append(converted_da)
                print(f"  - Converted '{da.name}' to kg/s (Example value for first time step, sum over space: {converted_da.isel(time=0).sum().item():.4e} kg/s)")

            try:
                # Summing all converted dry deposition DataArrays to get the total dry deposition.
                # xarray intelligently aligns dimensions before summing.
                total_dry_deposition = sum(dry_deposition_kg_s_per_species)
                total_dry_deposition.attrs['units'] = 'kg/s'
                total_dry_deposition.attrs['long_name'] = 'Total Dry Deposition of Mercury'
                print(f"\n  - Successfully summed all dry deposition species. Resultant DataArray shape: {total_dry_deposition.shape}")
                
                # Report a summary value for didactic purposes.
                if 'time' in total_dry_deposition.dims and total_dry_deposition['time'].size > 0:
                    summary_value_dry = total_dry_deposition.isel(time=0).sum().item()
                    print(f"  - Example total dry deposition (first time step, sum over all space): {summary_value_dry:.4e} kg/s")
                else: # For data without a 'time' dimension, sum over all existing dimensions.
                    summary_value_dry = total_dry_deposition.sum().item()
                    print(f"  - Example total dry deposition (overall sum over all space): {summary_value_dry:.4e} kg/s")
                print("  - The total dry deposition result is an xarray.DataArray.")
            except Exception as e:
                print(f"  - ERROR: Could not sum dry deposition variables. This might be due to incompatible dimensions or coordinates: {e}")
                total_dry_deposition = None

        # --- Step 4: Scientific Context & Variable Identification for Wet Deposition ---
        print("\n[STEP 4/6] Wet Deposition: Scientific Context & Variable Identification")
        print("-----------------------------------------------------------------------")
        print("Wet deposition is the removal of atmospheric species by precipitation\n"
              "(rain, snow, etc.). This process effectively 'washes out' soluble gases\n"
              "and particulate matter from the atmosphere. It's typically divided into\n"
              "two main categories in atmospheric models:\n"
              "  1. Convective Wet Deposition: Associated with strong, localized updrafts\n"
              "     like thunderstorms and cumulus clouds.\n"
              "  2. Large-Scale Wet Deposition: Associated with broader, more uniform\n"
              "     precipitation events like those from stratiform clouds and frontal systems.\n")
        print("In GEOS-Chem, wet deposition fluxes for mercury species are generally\n"
              "provided directly in 'kilograms per second' (kg/s) per grid cell,\n"
              "simplifying the summation process as no unit conversions are needed.\n")
        print("Identifying relevant wet deposition species from your list:\n"
              "  - Convective Wet Loss (prefixed with 'WetLossConv_Hg*')\n"
              "  - Large-Scale Wet Loss (prefixed with 'WetLossLS_Hg*')\n"
              "Note: Gaseous elemental mercury (Hg0) is generally considered to be\n"
              "poorly soluble in water and thus not significantly removed by wet deposition,\n"
              "which is consistent with its absence from the 'WetLossConv' and 'WetLossLS' lists.")

        # List of wet deposition species to include in total (from your CSV)
        wetd_conv_vars_to_extract = [
            'WetLossConv_Hg2ORGP',
            'WetLossConv_Hg2ClP',
            'WetLossConv_HgCl2',
            'WetLossConv_HgOHOH',
            'WetLossConv_HgClOH',
            'WetLossConv_HgBr2',
            'WetLossConv_HgBrOH'
            # Add or remove species based on your specific definition of "total wet deposition"
        ]
        wetd_ls_vars_to_extract = [
            'WetLossLS_Hg2ORGP',
            'WetLossLS_Hg2ClP',
            'WetLossLS_HgCl2',
            'WetLossLS_HgOHOH',
            'WetLossLS_HgClOH',
            'WetLossLS_HgBr2',
            'WetLossLS_HgBrOH'
            # Add or remove species based on your specific definition of "total wet deposition"
        ]

        extracted_wet_data_kg_s = [] # To store DataArrays (already in kg/s)

        print("\n  - Extracting from Convective Wet Loss dataset:")
        for var in wetd_conv_vars_to_extract:
            if var in ds_wetconv.data_vars:
                extracted_wet_data_kg_s.append(ds_wetconv[var])
                print(f"    - Found and extracted '{var}' (Shape: {ds_wetconv[var].shape}, Units: {ds_wetconv[var].attrs.get('units', 'N/A')})")
            else:
                print(f"    - WARNING: Variable '{var}' not found in WetLossConv dataset. This species will be skipped.")

        print("\n  - Extracting from Large-Scale Wet Loss dataset:")
        for var in wetd_ls_vars_to_extract:
            if var in ds_wetls.data_vars:
                extracted_wet_data_kg_s.append(ds_wetls[var])
                print(f"    - Found and extracted '{var}' (Shape: {ds_wetls[var].shape}, Units: {ds_wetls[var].attrs.get('units', 'N/A')})")
            else:
                print(f"    - WARNING: Variable '{var}' not found in WetLossLS dataset. This species will be skipped.")

        if not extracted_wet_data_kg_s:
            print("  - ERROR: No relevant wet deposition variables were found. Cannot calculate total wet deposition.")
            total_wet_deposition = None
        else:
            # --- Step 5: Calculation of Total Wet Deposition (in kg/s) ---
            print("\n[STEP 5/6] Wet Deposition: Summation of Species")
            print("-----------------------------------------------")
            print("Since wet deposition fluxes are already provided in consistent units (kg/s)\n"
                  "for both convective and large-scale processes, we simply sum the individual\n"
                  "species fluxes to derive the total wet deposition.\n")
            try:
                # Summing all wet deposition DataArrays (they are already in kg/s)
                total_wet_deposition = sum(extracted_wet_data_kg_s)
                total_wet_deposition.attrs['units'] = 'kg/s' # Explicitly set units for clarity
                total_wet_deposition.attrs['long_name'] = 'Total Wet Deposition of Mercury'
                print(f"\n  - Successfully summed all wet deposition species. Resultant DataArray shape: {total_wet_deposition.shape}")

                # Report a summary value for didactic purposes.
                if 'time' in total_wet_deposition.dims and total_wet_deposition['time'].size > 0:
                    summary_value_wet = total_wet_deposition.isel(time=0).sum().item()
                    print(f"  - Example total wet deposition (first time step, sum over all space): {summary_value_wet:.4e} kg/s")
                else: # For data without a 'time' dimension, sum over all existing dimensions.
                    summary_value_wet = total_wet_deposition.sum().item()
                    print(f"  - Example total wet deposition (overall sum over all space): {summary_value_wet:.4e} kg/s")
                print("  - The total wet deposition result is an xarray.DataArray.")
            except Exception as e:
                print(f"  - ERROR: Could not sum wet deposition variables. This might be due to incompatible dimensions or coordinates: {e}")
                total_wet_deposition = None

    finally:
        # --- Step 6: Resource Management - Close Datasets ---
        print("\n[STEP 6/6] Resource Management: Closing Datasets")
        print("-------------------------------------------------")
        print("It's crucial to close datasets after use to free up memory and system\n"
              "resources. This is especially important when processing large files or\n"
              "running multiple analyses to prevent resource exhaustion.\n")
        if ds_wetconv is not None:
            ds_wetconv.close()
            print("  - WetLossConv dataset closed.")
        if ds_wetls is not None:
            ds_wetls.close()
            print("  - WetLossLS dataset closed.")
        if ds_drydep is not None:
            ds_drydep.close()
            print("  - DryDep dataset closed.")
        print("Datasets closure attempt complete.")

    print("\n" + "="*80)
    print("--- TOTAL MERCURY DEPOSITION QUANTIFICATION COMPLETE ---".center(80))
    print("="*80 + "\n")
    return total_wet_deposition, total_dry_deposition


def quantify_hg2_to_ssa_transfer(mercurychem_file: str, statemet_file: str):
    """
    Quantifies the transfer rate of gaseous divalent mercury (Hg(II)) to sea salt aerosol (SSA).
    This function explains the origin and units of the 'Hg2GasToSSA' variable from GEOS-Chem,
    and demonstrates how to convert it to a more common unit like kg/s per grid cell.

    Args:
        mercurychem_file (str): Path to the GEOS-Chem MercuryChem_####.nc4 file,
                                containing the Hg2GasToSSA variable.
        statemet_file (str): Path to the GEOS-Chem StateMet_####.nc4 file,
                             containing meteorological variables like AIRVOL (grid cell volume)
                             which are necessary for unit conversion.

    Returns:
        xr.DataArray: DataArray representing the Hg(II) to SSA transfer rate in kg/s.
                      Returns None if files cannot be opened or critical variables are missing.
    """
    print("\n" + "="*80)
    print("--- QUANTIFYING Hg(II) GAS TO SEA SALT AEROSOL TRANSFER ---".center(80))
    print("="*80)

    ds_mercurychem, ds_statemet = None, None # Initialize datasets to None

    try:
        # --- Step 1: Data Acquisition - Open Required Datasets ---
        print("\n[STEP 1/4] Data Acquisition: Opening GEOS-Chem Output Files...")
        print("---------------------------------------------------------------")
        print("We need to open the MercuryChem file to access the pre-calculated\n"
              "Hg(II) to SSA transfer rate. Additionally, the StateMet file is needed\n"
              "to obtain grid cell volumes, which are essential for converting the\n"
              "volumetric transfer rate into a mass flux (kg/s).\n")

        # Open MercuryChem file
        print(f"  - Attempting to open MercuryChem file: {os.path.basename(mercurychem_file)}")
        try:
            ds_mercurychem = xr.open_dataset(mercurychem_file)
            print(f"    SUCCESS: Loaded MercuryChem dataset. Available variables: {list(ds_mercurychem.data_vars)}")
        except Exception as e:
            print(f"    ERROR: Could not open MercuryChem file '{mercurychem_file}'. Please check path and permissions: {e}")
            return None

        # Open StateMet file (specifically for Met_AIRVOL for unit conversion)
        print(f"  - Attempting to open StateMet file: {os.path.basename(statemet_file)}")
        try:
            ds_statemet = xr.open_dataset(statemet_file)
            print(f"    SUCCESS: Loaded StateMet dataset. Available variables: {list(ds_statemet.data_vars)}")
        except Exception as e:
            print(f"    ERROR: Could not open StateMet file '{statemet_file}'. Please check path and permissions: {e}")
            return None
        
        print("\n[STEP 1 Status]: All required datasets successfully opened.")

        # --- Step 2: Scientific Context & Variable Identification ---
        print("\n[STEP 2/4] Scientific Context: Hg(II) Gas to SSA Transfer")
        print("--------------------------------------------------------")
        print("The transfer of gaseous divalent mercury (Hg(II), or Hg2) to sea salt\n"
              "aerosol (SSA) represents a heterogeneous chemical process in the marine\n"
              "boundary layer, where gas-phase mercury reacts with or adsorbs onto the\n"
              "surface of aerosol particles. This is a key pathway for mercury removal\n"
              "from the gas phase.\n")
        print("In GEOS-Chem, this transfer is calculated internally as part of the\n"
              "mercury chemistry mechanism. The output variable 'Hg2GasToSSA' directly\n"
              "represents this rate, provided in units of 'molecules per cubic centimeter\n"
              "per second' (molec cm^-3 s^-1).\n")
        print("To make this rate physically meaningful for mass budget analyses, we need\n"
              "to convert it from a volumetric molecular rate to a mass flux in 'kilograms\n"
              "per second' (kg/s) for each model grid cell. This requires knowledge of\n"
              "the molar mass of mercury, Avogadro's number, and the volume of the grid cell.\n")
        print("The conversion formula from 'molec/cm^3/s' to 'kg/s' for a grid cell is:\n")
        print("  $$ \\text{Rate}_{\\text{kg/s}} = \\text{Rate}_{\\text{molec/cm}^3\\text{/s}} \\times \\left( \\frac{\\text{M}_{\\text{Hg}}}{\\text{N}_{\\text{A}}} \\right) \\times \\text{Volume}_{\\text{grid, cm}^3} $$")
        print("Where:\n"
              f"  - $M_{{Hg}}$: Molar Mass of Mercury ({MOLAR_MASS_HG_KG_PER_MOL} kg/mol)\n"
              f"  - $N_{{A}}$: Avogadro's Number ({AVOGADRO_NUMBER:.3e} molec/mol)\n"
              "  - $\\text{Volume}_{\\text{grid, cm}^3}$: Volume of the grid cell in cm$^3$ (obtained from 'Met_AIRVOL').\n")

        # Identify the variable for the direct transfer rate from MercuryChem
        hg2_to_ssa_transfer_rate_var_name = 'Hg2GasToSSA'
        transfer_rate_molec_cm3_s = None

        if hg2_to_ssa_transfer_rate_var_name in ds_mercurychem.data_vars:
            transfer_rate_molec_cm3_s = ds_mercurychem[hg2_to_ssa_transfer_rate_var_name]
            print(f"  - Found and extracted '{hg2_to_ssa_transfer_rate_var_name}' (Shape: {transfer_rate_molec_cm3_s.shape}, Units: {transfer_rate_molec_cm3_s.attrs.get('units', 'N/A')})")
        else:
            print(f"  - ERROR: Variable '{hg2_to_ssa_transfer_rate_var_name}' not found in MercuryChem dataset. Cannot quantify transfer.")
            return None

        # Get the grid cell volume from the StateMet dataset for unit conversion
        if 'Met_AIRVOL' in ds_statemet.data_vars:
            grid_volume_m3 = ds_statemet['Met_AIRVOL'] # 'Met_AIRVOL' is typically in m^3 in GEOS-Chem output.
            print(f"\n  - Extracted 'Met_AIRVOL' from StateMet dataset (Shape: {grid_volume_m3.shape}, Units: {grid_volume_m3.attrs.get('units', 'N/A')}).")
            # Convert AIRVOL from m^3 to cm^3 to match the units of the flux (molec/cm^3/s).
            grid_volume_cm3 = grid_volume_m3 / CM3_TO_M3 # Equivalent to grid_volume_m3 * 1e6
            print(f"    -> Converted 'Met_AIRVOL' from m^3 to cm^3 for unit consistency in the calculation.")
        else:
            print("  - CRITICAL ERROR: 'Met_AIRVOL' variable (grid cell volume) not found in StateMet dataset. Cannot convert transfer rate from molec/cm^3/s to kg/s.")
            return None
        
        # --- Step 3: Calculation of Transfer Rate (in kg/s) ---
        print("\n[STEP 3/4] Calculation: Converting Volumetric Rate to Mass Flux (kg/s)")
        print("----------------------------------------------------------------------")
        print("Here, we perform the unit conversion. The direct transfer rate from the\n"
              "model ($Hg2GasToSSA$) in molec cm^-3 s^-1 is converted to a mass flux in\n"
              "kg/s for each grid cell using the formula defined above. This provides\n"
              "a more intuitive measure of the mass of mercury being transferred.\n")

        hg2_to_ssa_transfer_rate_kg_s = None
        if transfer_rate_molec_cm3_s is not None and grid_volume_cm3 is not None:
            try:
                # Apply the conversion formula:
                # (molec/cm^3/s) * (kg/mol) / (molec/mol) * (cm^3) = kg/s
                hg2_to_ssa_transfer_rate_kg_s = transfer_rate_molec_cm3_s * \
                                                 (MOLAR_MASS_HG_KG_PER_MOL / AVOGADRO_NUMBER) * \
                                                 grid_volume_cm3
                
                # Assign meaningful attributes to the resulting DataArray
                hg2_to_ssa_transfer_rate_kg_s.attrs['units'] = 'kg/s'
                hg2_to_ssa_transfer_rate_kg_s.attrs['long_name'] = 'Hg(II) Gas to Sea Salt Aerosol Transfer Rate'

                print(f"\n  - Calculation complete. Resultant DataArray shape: {hg2_to_ssa_transfer_rate_kg_s.shape}")
                
                # Report a summary value for didactic purposes.
                if 'time' in hg2_to_ssa_transfer_rate_kg_s.dims and 'lev' in hg2_to_ssa_transfer_rate_kg_s.dims and \
                   hg2_to_ssa_transfer_rate_kg_s['time'].size > 0 and hg2_to_ssa_transfer_rate_kg_s['lev'].size > 0:
                    summary_value = hg2_to_ssa_transfer_rate_kg_s.isel(time=0, lev=0).sum().item()
                    print(f"  - Example value (first time step, first level, sum over space): {summary_value:.4e} kg/s")
                elif 'time' in hg2_to_ssa_transfer_rate_kg_s.dims and hg2_to_ssa_transfer_rate_kg_s['time'].size > 0:
                    summary_value = hg2_to_ssa_transfer_rate_kg_s.isel(time=0).sum().item()
                    print(f"  - Example value (first time step, sum over space): {summary_value:.4e} kg/s")
                else: # For data without time or level dimensions, sum over all existing dimensions.
                    summary_value = hg2_to_ssa_transfer_rate_kg_s.sum().item()
                    print(f"  - Example value (overall sum over space): {summary_value:.4e} kg/s")
                
                print("  - The calculated transfer rate result is an xarray.DataArray.")

            except Exception as e:
                print(f"  - ERROR: Calculation for Hg(II) to SSA transfer failed. This might be due to incompatible dimensions or coordinates: {e}")
                hg2_to_ssa_transfer_rate_kg_s = None
        else:
            print("  - Skipping Hg(II) to SSA transfer calculation due to missing primary variables.")

    finally:
        # --- Step 4: Resource Management - Close Datasets ---
        print("\n[STEP 4/4] Resource Management: Closing Datasets")
        print("-----------------------------------------------")
        print("Closing datasets is crucial to release memory and system resources.\n"
              "This ensures efficient use of your computer's memory and prevents issues\n"
              "when running subsequent data analyses.\n")
        if ds_mercurychem is not None:
            ds_mercurychem.close()
            print("  - MercuryChem dataset closed.")
        if ds_statemet is not None:
            ds_statemet.close()
            print("  - StateMet dataset closed.")
        print("Datasets closure attempt complete.")

    print("\n" + "="*80)
    print("--- Hg(II) GAS TO SEA SALT AEROSOL TRANSFER QUANTIFICATION COMPLETE ---".center(80))
    print("="*80 + "\n")
    return hg2_to_ssa_transfer_rate_kg_s
