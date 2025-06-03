# mercury_analysis_2.py
import xarray as xr
import numpy as np

# --- Define Universal Constants ---
# These are fundamental physical constants used in atmospheric chemistry calculations.
# Displaying them explicitly helps the user understand the numerical values involved
# in the unit conversions and scientific calculations.
MOLAR_MASS_HG_G_PER_MOL = 200.59  # grams per mole (Molar Mass of Mercury)
MOLAR_MASS_HG_KG_PER_MOL = MOLAR_MASS_HG_G_PER_MOL / 1000  # kilograms per mole

AVOGADRO_NUMBER = 6.022e23  # molecules per mole (or atoms per mole for elemental species)

# Conversion factors for spatial units
CM2_TO_M2 = 1e-4  # Conversion factor: 1 cm^2 = 1e-4 m^2
CM3_TO_M3 = 1e-6  # Conversion factor: 1 cm^3 = 1e-6 m^3

# Print out constants for didactic purposes at the module level
print("\n--- Universal Constants Used in Calculations ---")
print(f"  Molar Mass of Mercury (Hg): {MOLAR_MASS_HG_G_PER_MOL} g/mol ({MOLAR_MASS_HG_KG_PER_MOL} kg/mol)")
print(f"  Avogadro's Number (N_A): {AVOGADRO_NUMBER:.3e} molecules/mol")
print(f"  Unit Conversion: 1 cm^2 = {CM2_TO_M2} m^2")
print(f"  Unit Conversion: 1 cm^3 = {CM3_TO_M3} m^3")
print("----------------------------------------------\n")

def quantify_total_deposition(ds_wetconv: xr.Dataset, ds_wetls: xr.Dataset, ds_drydep: xr.Dataset) -> tuple[xr.DataArray, xr.DataArray]:
    """
    Quantifies the total wet and dry deposition of mercury species from GEOS-Chem output.
    This function aims to be didactic, explaining the scientific principles,
    mathematical steps, and unit conversions involved in calculating total
    mercury deposition.

    Args:
        ds_wetconv (xr.Dataset): The opened GEOS-Chem WetLossConv dataset.
        ds_wetls (xr.Dataset): The opened GEOS-Chem WetLossLS dataset.
        ds_drydep (xr.Dataset): The opened GEOS-Chem DryDep dataset.

    Returns:
        tuple[xr.DataArray, xr.DataArray]: A tuple containing:
            - total_wet_deposition (xr.DataArray): Total wet deposition flux in kg/s.
            - total_dry_deposition (xr.DataArray): Total dry deposition flux in kg/s.
            Returns (None, None) if critical variables are missing.
    """
    print("\n" + "="*80)
    print("--- QUANTIFYING TOTAL MERCURY DEPOSITION ---".center(80))
    print("="*80)

    total_dry_deposition = None
    total_wet_deposition = None

    # --- Step 1: Scientific Context & Variable Identification for Dry Deposition ---
    print("\n[STEP 1/4] Dry Deposition: Scientific Context & Variable Identification")
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

    dryd_vars_to_extract = [
        'DryDep_Hg0', 'DryDep_HgCl2', 'DryDep_HgOHOH', 'DryDep_HgClOH',
        'DryDep_HgBr2', 'DryDep_HgBrOH', 'DryDep_Hg2ORGP', 'DryDep_Hg2ClP'
    ]
    extracted_dry_data_molec_cm2_s = []

    for var in dryd_vars_to_extract:
        if var in ds_drydep.data_vars:
            extracted_dry_data_molec_cm2_s.append(ds_drydep[var])
            print(f"  - Found and extracted '{var}' (Shape: {ds_drydep[var].shape}, Units: {ds_drydep[var].attrs.get('units', 'N/A')})")
        else:
            print(f"  - WARNING: Variable '{var}' not found in DryDep dataset. This species will be skipped for the dry deposition sum.")

    if not extracted_dry_data_molec_cm2_s:
        print("  - ERROR: No relevant dry deposition variables were found based on the provided list. Cannot calculate total dry deposition.")
        return None, None
    
    # Get the grid cell area from the DryDep dataset for unit conversion
    if 'AREA' in ds_drydep.data_vars:
        grid_area_m2 = ds_drydep['AREA'] # The 'AREA' variable is typically in m^2 in GEOS-Chem output.
        print(f"\n  - Extracted 'AREA' from DryDep dataset (Shape: {grid_area_m2.shape}, Units: {grid_area_m2.attrs.get('units', 'N/A')}).")
        grid_area_cm2 = grid_area_m2 / CM2_TO_M2 # Convert AREA from m^2 to cm^2
        print(f"    -> Converted 'AREA' from m^2 to cm^2 for unit consistency in the calculation.")
    else:
        print("  - CRITICAL ERROR: 'AREA' variable (grid cell surface area) not found in DryDep dataset. Cannot convert dry deposition flux from molec/cm^2/s to kg/s.")
        return None, None

    # --- Step 2: Calculation of Total Dry Deposition (in kg/s) ---
    print("\n[STEP 2/4] Dry Deposition: Performing Unit Conversion and Summation")
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
        print(f"  - Converted '{da.name}' to kg/s.")

    try:
        total_dry_deposition = sum(dry_deposition_kg_s_per_species)
        total_dry_deposition.attrs['units'] = 'kg/s'
        total_dry_deposition.attrs['long_name'] = 'Total Dry Deposition of Mercury'
        print(f"\n  - Successfully summed all dry deposition species. Resultant DataArray shape: {total_dry_deposition.shape}")
        if 'time' in total_dry_deposition.dims and total_dry_deposition['time'].size > 0:
            summary_value_dry = total_dry_deposition.isel(time=0).sum().item()
            print(f"  - Example total dry deposition (first time step, sum over all space): {summary_value_dry:.4e} kg/s")
        else:
            summary_value_dry = total_dry_deposition.sum().item()
            print(f"  - Example total dry deposition (overall sum over all space): {summary_value_dry:.4e} kg/s")
        print("  - The total dry deposition result is an xarray.DataArray.")
    except Exception as e:
        print(f"  - ERROR: Could not sum dry deposition variables. This might be due to incompatible dimensions or coordinates: {e}")
        total_dry_deposition = None


    # --- Step 3: Scientific Context & Variable Identification for Wet Deposition ---
    print("\n[STEP 3/4] Wet Deposition: Scientific Context & Variable Identification")
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
          "  - Large-Scale Wet Loss (prefixed with 'WetLossLS_Hg*')\n")

    wetd_conv_vars_to_extract = [
        'WetLossConv_Hg2ORGP', 'WetLossConv_Hg2ClP', 'WetLossConv_HgCl2',
        'WetLossConv_HgOHOH', 'WetLossConv_HgClOH', 'WetLossConv_HgBr2',
        'WetLossConv_HgBrOH'
    ]
    wetd_ls_vars_to_extract = [
        'WetLossLS_Hg2ORGP', 'WetLossLS_Hg2ClP', 'WetLossLS_HgCl2',
        'WetLossLS_HgOHOH', 'WetLossLS_HgClOH', 'WetLossLS_HgBr2',
        'WetLossLS_HgBrOH'
    ]

    extracted_wet_data_kg_s = []

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
        # --- Step 4: Calculation of Total Wet Deposition (in kg/s) ---
        print("\n[STEP 4/4] Wet Deposition: Summation of Species")
        print("-----------------------------------------------")
        print("Since wet deposition fluxes are already provided in consistent units (kg/s)\n"
              "for both convective and large-scale processes, we simply sum the individual\n"
              "species fluxes to derive the total wet deposition.\n")
        try:
            total_wet_deposition = sum(extracted_wet_data_kg_s)
            total_wet_deposition.attrs['units'] = 'kg/s'
            total_wet_deposition.attrs['long_name'] = 'Total Wet Deposition of Mercury'
            print(f"\n  - Successfully summed all wet deposition species. Resultant DataArray shape: {total_wet_deposition.shape}")

            if 'time' in total_wet_deposition.dims and total_wet_deposition['time'].size > 0:
                summary_value_wet = total_wet_deposition.isel(time=0).sum().item()
                print(f"  - Example total wet deposition (first time step, sum over all space): {summary_value_wet:.4e} kg/s")
            else:
                summary_value_wet = total_wet_deposition.sum().item()
                print(f"  - Example total wet deposition (overall sum over all space): {summary_value_wet:.4e} kg/s")
            print("  - The total wet deposition result is an xarray.DataArray.")
        except Exception as e:
            print(f"  - ERROR: Could not sum wet deposition variables. This might be due to incompatible dimensions or coordinates: {e}")
            total_wet_deposition = None

    print("\n" + "="*80)
    print("--- TOTAL MERCURY DEPOSITION QUANTIFICATION COMPLETE ---".center(80))
    print("="*80 + "\n")
    return total_wet_deposition, total_dry_deposition


def quantify_hg2_to_ssa_transfer(ds_mercurychem: xr.Dataset, ds_statemet: xr.Dataset) -> xr.DataArray:
    """
    Quantifies the transfer rate of gaseous divalent mercury (Hg(II)) to sea salt aerosol (SSA).
    This function explains the origin and units of the 'Hg2GasToSSA' variable from GEOS-Chem,
    and demonstrates how to convert it to a more common unit like kg/s per grid cell.

    Args:
        ds_mercurychem (xr.Dataset): The opened GEOS-Chem MercuryChem dataset.
        ds_statemet (xr.Dataset): The opened GEOS-Chem StateMet dataset.

    Returns:
        xr.DataArray: DataArray representing the Hg(II) to SSA transfer rate in kg/s.
                      Returns None if critical variables are missing.
    """
    print("\n" + "="*80)
    print("--- QUANTIFYING Hg(II) GAS TO SEA SALT AEROSOL TRANSFER ---".center(80))
    print("="*80)

    hg2_to_ssa_transfer_rate_kg_s = None

    # --- Step 1: Scientific Context & Variable Identification ---
    print("\n[STEP 1/3] Scientific Context: Hg(II) Gas to SSA Transfer")
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

    hg2_to_ssa_transfer_rate_var_name = 'Hg2GasToSSA'
    transfer_rate_molec_cm3_s = None

    if hg2_to_ssa_transfer_rate_var_name in ds_mercurychem.data_vars:
        transfer_rate_molec_cm3_s = ds_mercurychem[hg2_to_ssa_transfer_rate_var_name]
        print(f"  - Found and extracted '{hg2_to_ssa_transfer_rate_var_name}' (Shape: {transfer_rate_molec_cm3_s.shape}, Units: {transfer_rate_molec_cm3_s.attrs.get('units', 'N/A')})")
    else:
        print(f"  - ERROR: Variable '{hg2_to_ssa_transfer_rate_var_name}' not found in MercuryChem dataset. Cannot quantify transfer.")
        return None

    if 'Met_AIRVOL' in ds_statemet.data_vars:
        grid_volume_m3 = ds_statemet['Met_AIRVOL']
        print(f"\n  - Extracted 'Met_AIRVOL' from StateMet dataset (Shape: {grid_volume_m3.shape}, Units: {grid_volume_m3.attrs.get('units', 'N/A')}).")
        grid_volume_cm3 = grid_volume_m3 / CM3_TO_M3
        print(f"    -> Converted 'Met_AIRVOL' from m^3 to cm^3 for unit consistency in the calculation.")
    else:
        print("  - CRITICAL ERROR: 'Met_AIRVOL' variable (grid cell volume) not found in StateMet dataset. Cannot convert transfer rate from molec/cm^3/s to kg/s.")
        return None
    
    # --- Step 2: Calculation of Transfer Rate (in kg/s) ---
    print("\n[STEP 2/3] Calculation: Converting Volumetric Rate to Mass Flux (kg/s)")
    print("----------------------------------------------------------------------")
    print("Here, we perform the unit conversion. The direct transfer rate from the\n"
          "model ($Hg2GasToSSA$) in molec cm^-3 s^-1 is converted to a mass flux in\n"
          "kg/s for each grid cell using the formula defined above. This provides\n"
          "a more intuitive measure of the mass of mercury being transferred.\n")

    if transfer_rate_molec_cm3_s is not None and grid_volume_cm3 is not None:
        try:
            # Apply the conversion formula:
            # (molec/cm^3/s) * (kg/mol) / (molec/mol) * (cm^3) = kg/s
            hg2_to_ssa_transfer_rate_kg_s = transfer_rate_molec_cm3_s * \
                                            (MOLAR_MASS_HG_KG_PER_MOL / AVOGADRO_NUMBER) * \
                                            grid_volume_cm3
            
            hg2_to_ssa_transfer_rate_kg_s.attrs['units'] = 'kg/s'
            hg2_to_ssa_transfer_rate_kg_s.attrs['long_name'] = 'Hg(II) Gas to Sea Salt Aerosol Transfer Rate'

            print(f"\n  - Calculation complete. Resultant DataArray shape: {hg2_to_ssa_transfer_rate_kg_s.shape}")
            if 'time' in hg2_to_ssa_transfer_rate_kg_s.dims and 'lev' in hg2_to_ssa_transfer_rate_kg_s.dims and \
               hg2_to_ssa_transfer_rate_kg_s['time'].size > 0 and hg2_to_ssa_transfer_rate_kg_s['lev'].size > 0:
                summary_value = hg2_to_ssa_transfer_rate_kg_s.isel(time=0, lev=0).sum().item()
                print(f"  - Example value (first time step, first level, sum over space): {summary_value:.4e} kg/s")
            elif 'time' in hg2_to_ssa_transfer_rate_kg_s.dims and hg2_to_ssa_transfer_rate_kg_s['time'].size > 0:
                summary_value = hg2_to_ssa_transfer_rate_kg_s.isel(time=0).sum().item()
                print(f"  - Example value (first time step, sum over space): {summary_value:.4e} kg/s")
            else:
                summary_value = hg2_to_ssa_transfer_rate_kg_s.sum().item()
                print(f"  - Example value (overall sum over space): {summary_value:.4e} kg/s")
            print("  - The calculated transfer rate result is an xarray.DataArray.")

        except Exception as e:
            print(f"  - ERROR: Calculation for Hg(II) to SSA transfer failed. This might be due to incompatible dimensions or coordinates: {e}")
            hg2_to_ssa_transfer_rate_kg_s = None
    else:
        print("  - Skipping Hg(II) to SSA transfer calculation due to missing primary variables.")

    print("\n" + "="*80)
    print("--- Hg(II) GAS TO SEA SALT AEROSOL TRANSFER QUANTIFICATION COMPLETE ---".center(80))
    print("="*80 + "\n")
    return hg2_to_ssa_transfer_rate_kg_s
