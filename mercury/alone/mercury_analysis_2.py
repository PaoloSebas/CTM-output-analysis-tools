# mercury_analysis_2.py

import xarray as xr
import numpy as np

# --- Constants ---
MOLAR_MASS_HG_KG_PER_MOL = 0.20059  # Molar mass of Mercury in kg/mol
AVOGADRO_NUMBER = 6.022e+23        # Avogadro's number in molec/mol
CM2_TO_M2 = 1e-4                   # Conversion factor from cm^2 to m^2
CM3_TO_M3 = 1e-6                   # Conversion factor from cm^3 to m^3

def quantify_total_deposition(ds_mercurychem: xr.Dataset | None,
                              ds_statemet: xr.Dataset | None,
                              ds_drydep: xr.Dataset | None,
                              ds_wetconv: xr.Dataset | None,
                              ds_wetls: xr.Dataset | None
                             ) -> tuple[xr.DataArray | None, xr.DataArray | None]:
    """
    Quantifies total mercury dry and wet deposition fluxes from GEOS-Chem data.

    Args:
        ds_mercurychem (xr.Dataset | None): The opened GEOS-Chem MercuryChem dataset.
        ds_statemet (xr.Dataset | None): The opened GEOS-Chem StateMet dataset (for AREA).
        ds_drydep (xr.Dataset | None): The opened GEOS-Chem DryDep dataset (for dry deposition species).
        ds_wetconv (xr.Dataset | None): The opened GEOS-Chem WetLossConv dataset (for convective wet deposition species).
        ds_wetls (xr.Dataset | None): The opened GEOS-Chem WetLossLS dataset (for large-scale wet deposition species).

    Returns:
        tuple[xr.DataArray | None, xr.DataArray | None]:
            - total_dry_deposition_kg_s (xr.DataArray): Total dry deposition in kg/s.
            - total_wet_deposition_kg_s (xr.DataArray): Total wet deposition in kg/s.
            Returns None for either if calculation fails or data is missing.
    """
    print("\n" + "="*80)
    print("--- QUANTIFYING TOTAL MERCURY DEPOSITION ---".center(80))
    print("="*80)

    total_dry_deposition = None
    total_wet_deposition = None

    # --- Step 1/4: Dry Deposition: Scientific Context & Variable Identification ---
    print("\n[STEP 1/4] Dry Deposition: Scientific Context & Variable Identification")
    print("--------------------------------------------------------")
    print("Dry deposition is the direct transfer of gases and aerosols from the")
    print("atmosphere to the Earth's surface in the absence of precipitation.")
    print("It's a continuous process governed by turbulence and surface properties.")
    print("In GEOS-Chem, dry deposition fluxes are often provided in units of")
    print("'molecules per square centimeter per second' (molec cm^-2 s^-1).\n")
    print("To calculate total dry deposition in 'kilograms per second' (kg/s) for the")
    print("entire model domain (or a specific grid cell), we need to sum the individual")
    print("species and convert their units from a number flux to a mass flux.\n")
    print("The conversion formula from 'molecules/cm^2/s' to 'kg/s' for a specific grid cell is:\n")
    print("  $$\\text{Flux}_{\\text{dry, kg/s}} = \\text{Flux}_{\\text{dry, molec/cm}^2\\text{/s}} \\times \\left( \\frac{\\text{M}_{\\text{Hg}}}{\\text{N}_{\\text{A}}} \\right) \\times \\text{Area}_{\\text{grid, cm}^2}$$")
    print("Where:")
    print(f"  - $M_{{Hg}}$: Molar Mass of Mercury ({MOLAR_MASS_HG_KG_PER_MOL} kg/mol)")
    print(f"  - $N_{{A}}$: Avogadro's Number ({AVOGADRO_NUMBER:.3e} molec/mol)")
    print("  - $\\text{Area}_{\\text{grid, cm}^2}$: Surface area of the grid cell in cm$^2$ (obtained from 'AREA' variable).\n")

    dry_deposition_species_prefixes = ['DryDep_Hg0', 'DryDep_HgCl2', 'DryDep_HgOHOH', 'DryDep_HgClOH',
                                       'DryDep_HgBr2', 'DryDep_HgBrOH', 'DryDep_Hg2ORGP', 'DryDep_Hg2ClP']
    dry_deposition_kg_s_per_species = []
    grid_area_cm2 = None # Initialize outside conditional block

    if ds_statemet is None:
        print("  - CRITICAL ERROR: **StateMet** dataset not provided. Cannot get 'AREA' variable needed for dry deposition calculation.")
        total_dry_deposition = None
    else:
        if 'AREA' in ds_statemet.data_vars:
            grid_area_m2 = ds_statemet['AREA']
            print(f"\n  - Extracted 'AREA' from StateMet dataset (Shape: {grid_area_m2.shape}, Units: {grid_area_m2.attrs.get('units', 'N/A')}).")
            grid_area_cm2 = grid_area_m2 / CM2_TO_M2 # Convert AREA from m^2 to cm^2
            print(f"    -> Converted 'AREA' from m^2 to cm^2 for unit consistency in the calculation.")
        else:
            print("  - CRITICAL ERROR: 'AREA' variable (grid cell area) not found in StateMet dataset. Cannot calculate dry deposition in kg/s.")
            total_dry_deposition = None # Critical dependency missing

    if ds_drydep is None:
        print("\n  - CRITICAL ERROR: **DryDep** dataset not provided. Cannot find dry deposition species.")
        total_dry_deposition = None # Ensure it's None
    elif grid_area_cm2 is None: # AREA was not found or converted
        print("\n  - CRITICAL ERROR: Grid area ('AREA') is missing or could not be converted. Cannot calculate dry deposition.")
        total_dry_deposition = None
    else:
        print("\nIdentifying relevant dry deposition species (prefixed with 'DryDep_') from the provided list:")
        found_any_dry_species = False
        for var_name in dry_deposition_species_prefixes:
            if var_name in ds_drydep.data_vars: # *** Corrected: Look in ds_drydep ***
                var_data = ds_drydep[var_name]
                print(f"  - Found and extracted '{var_name}' (Shape: {var_data.shape}, Units: {var_data.attrs.get('units', 'N/A')})")
                dry_deposition_kg_s_per_species.append(var_data)
                found_any_dry_species = True
            else:
                print(f"  - Warning: '{var_name}' not found in DryDep dataset. Skipping.")
        
        if not found_any_dry_species:
             print("  - No dry deposition species found in DryDep dataset. Skipping dry deposition calculation.")
             total_dry_deposition = None # No species found, so total dry dep will be None

    # --- Step 2/4: Dry Deposition: Performing Unit Conversion and Summation ---
    print("\n[STEP 2/4] Dry Deposition: Performing Unit Conversion and Summation")
    print("-------------------------------------------------------------------")
    print("Each dry deposition flux (in molec cm^-2 s^-1) is converted to a mass")
    print("flux in 'kg per second' for each grid cell. This involves multiplying")
    print("by the molar mass of mercury ($M_{{Hg}}$), dividing by Avogadro's number ($N_{{A}}$),")
    print("and multiplying by the surface area of the grid cell ($\\text{Area}_{\\text{grid, cm}^2}$). ")
    print("Finally, these mass fluxes for all individual mercury species are summed.\n")

    # Only proceed if we have species and a valid grid area
    if dry_deposition_kg_s_per_species and grid_area_cm2 is not None:
        try:
            converted_dry_species = []
            for species_flux_molec_cm2_s in dry_deposition_kg_s_per_species:
                # Apply the conversion formula:
                # (molec/cm^2/s) * (kg/mol) / (molec/mol) * (cm^2) = kg/s
                converted_flux_kg_s = species_flux_molec_cm2_s * \
                                      (MOLAR_MASS_HG_KG_PER_MOL / AVOGADRO_NUMBER) * \
                                      grid_area_cm2

                converted_flux_kg_s.attrs['units'] = 'kg/s'
                converted_flux_kg_s.attrs['long_name'] = f"Dry Deposition of {species_flux_molec_cm2_s.name} in kg/s"
                converted_dry_species.append(converted_flux_kg_s)
                print(f"  - Converted '{species_flux_molec_cm2_s.name}' to kg/s.")

            total_dry_deposition = sum(converted_dry_species)
            total_dry_deposition.attrs['long_name'] = 'Total Dry Deposition of Mercury'
            total_dry_deposition.attrs['units'] = 'kg/s'

            print(f"\n  - Successfully summed all dry deposition species. Resultant DataArray shape: {total_dry_deposition.shape}")

            if 'time' in total_dry_deposition.dims and total_dry_deposition['time'].size > 0:
                summary_value_dry = total_dry_deposition.isel(time=0).sum().compute().item()
                print(f"  - Example total dry deposition (first time step, sum over space): {summary_value_dry:.4e} kg/s")
            else:
                summary_value_dry = total_dry_deposition.sum().compute().item()
                print(f"  - Example total dry deposition (overall sum over space): {summary_value_dry:.4e} kg/s")
            print("  - The total dry deposition result is an xarray.DataArray.")

        except Exception as e:
            print(f"  - ERROR: Could not sum dry deposition variables. This might be due to incompatible dimensions or coordinates: {e}")
            total_dry_deposition = None
    else:
        # This else block is reached if dry_deposition_kg_s_per_species is empty or grid_area_cm2 is None
        print("  - Skipping dry deposition calculation due to missing species data or area.")
        total_dry_deposition = None # Ensure it's None if skipped

    # --- Step 3/4: Wet Deposition: Scientific Context & Variable Identification ---
    print("\n[STEP 3/4] Wet Deposition: Scientific Context & Variable Identification")
    print("-------------------------------------------------------")
    print("Wet deposition is the removal of atmospheric species by precipitation")
    print("(rain, snow, etc.). This process effectively 'washes out' soluble gases")
    print("and particulate matter from the atmosphere. It's typically divided into")
    print("two main categories in atmospheric models:")
    print("  1. Convective Wet Deposition: Associated with strong, localized updrafts")
    print("     like thunderstorms and cumulus clouds.")
    print("  2. Large-Scale Wet Deposition: Associated with broader, more uniform")
    print("     precipitation events like those from stratiform clouds and frontal systems.\n")
    print("In GEOS-Chem, wet deposition fluxes for mercury species are generally")
    print("provided directly in 'kilograms per second' (kg/s) per grid cell,")
    print("simplifying the summation process as no unit conversions are needed.\n")
    print("\nIdentifying relevant wet deposition species from your list:")
    print("  - Convective Wet Loss (prefixed with 'WetLossConv_Hg*')")
    print("  - Large-Scale Wet Loss (prefixed with 'WetLossLS_Hg*')\n")

    convective_wet_loss_species = ['WetLossConv_Hg2ORGP', 'WetLossConv_Hg2ClP', 'WetLossConv_HgCl2',
                                   'WetLossConv_HgOHOH', 'WetLossConv_HgClOH', 'WetLossConv_HgBr2',
                                   'WetLossConv_HgBrOH']
    large_scale_wet_loss_species = ['WetLossLS_Hg2ORGP', 'WetLossLS_Hg2ClP', 'WetLossLS_HgCl2',
                                    'WetLossLS_HgOHOH', 'WetLossLS_HgClOH', 'WetLossLS_HgBr2',
                                    'WetLossLS_HgBrOH']

    all_wet_deposition_species = []

    if ds_wetconv is None:
        print("  - CRITICAL ERROR: **WetLossConv** dataset not provided. Cannot find convective wet deposition species.")
    else:
        print("  - Extracting from Convective Wet Loss dataset:")
        found_any_conv_species = False
        for var_name in convective_wet_loss_species:
            if var_name in ds_wetconv.data_vars: # *** Corrected: Look in ds_wetconv ***
                var_data = ds_wetconv[var_name]
                print(f"    - Found and extracted '{var_name}' (Shape: {var_data.shape}, Units: {var_data.attrs.get('units', 'N/A')})")
                all_wet_deposition_species.append(var_data)
                found_any_conv_species = True
            else:
                print(f"    - Warning: '{var_name}' not found in WetLossConv dataset. Skipping.")
        if not found_any_conv_species:
             print("  - No convective wet deposition species found in WetLossConv dataset.")

    if ds_wetls is None:
        print("\n  - CRITICAL ERROR: **WetLossLS** dataset not provided. Cannot find large-scale wet deposition species.")
    else:
        print("\n  - Extracting from Large-Scale Wet Loss dataset:")
        found_any_ls_species = False
        for var_name in large_scale_wet_loss_species:
            if var_name in ds_wetls.data_vars: # *** Corrected: Look in ds_wetls ***
                var_data = ds_wetls[var_name]
                print(f"    - Found and extracted '{var_name}' (Shape: {var_data.shape}, Units: {var_data.attrs.get('units', 'N/A')})")
                all_wet_deposition_species.append(var_data)
                found_any_ls_species = True
            else:
                print(f"    - Warning: '{var_name}' not found in WetLossLS dataset. Skipping.")
        if not found_any_ls_species:
             print("  - No large-scale wet deposition species found in WetLossLS dataset.")

    # --- Step 4/4: Wet Deposition: Summation of Species ---
    print("\n[STEP 4/4] Wet Deposition: Summation of Species")
    print("-----------------------------------------------")
    print("Since wet deposition fluxes are already provided in consistent units (kg/s)")
    print("for both convective and large-scale processes, we simply sum the individual")
    print("species fluxes to derive the total wet deposition.\n")

    if all_wet_deposition_species:
        try:
            total_wet_deposition = sum(all_wet_deposition_species)
            total_wet_deposition.attrs['long_name'] = 'Total Wet Deposition of Mercury'
            total_wet_deposition.attrs['units'] = 'kg/s'

            print(f"\n  - Successfully summed all wet deposition species. Resultant DataArray shape: {total_wet_deposition.shape}")

            if 'time' in total_wet_deposition.dims and total_wet_deposition['time'].size > 0:
                summary_value_wet = total_wet_deposition.isel(time=0).sum().compute().item()
                print(f"  - Example total wet deposition (first time step, sum over space): {summary_value_wet:.4e} kg/s")
            else:
                summary_value_wet = total_wet_deposition.sum().compute().item()
                print(f"  - Example total wet deposition (overall sum over space): {summary_value_wet:.4e} kg/s")
            print("  - The total wet deposition result is an xarray.DataArray.")

        except Exception as e:
            print(f"  - ERROR: Could not sum wet deposition variables. This might be due to incompatible dimensions or coordinates: {e}")
            total_wet_deposition = None
    else:
        print("  - Skipping wet deposition calculation due to missing species data.")
        total_wet_deposition = None # Ensure it's None if skipped

    print("\n" + "="*80)
    print("--- TOTAL MERCURY DEPOSITION QUANTIFICATION COMPLETE ---".center(80))
    print("="*80 + "\n")

    return total_dry_deposition, total_wet_deposition


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

    # Initialize the variable correctly
    hg2_to_ssa_transfer_rate_kg_s = None # Correct variable name for the result

    # --- Step 1/3: Scientific Context & Variable Identification ---
    print("\n[STEP 1/3] Scientific Context: Hg(II) Gas to SSA Transfer")
    print("--------------------------------------------------------")
    print("The transfer of gaseous divalent mercury (Hg(II), or Hg2) to sea salt")
    print("aerosol (SSA) represents a heterogeneous chemical process in the marine")
    print("boundary layer, where gas-phase mercury reacts with or adsorbs onto the")
    print("surface of aerosol particles. This is a key pathway for mercury removal")
    print("from the gas phase.\n")
    print("In GEOS-Chem, this transfer is calculated internally as part of the")
    print("mercury chemistry mechanism. The output variable 'Hg2GasToSSA' directly")
    print("represents this rate, provided in units of 'molecules per cubic centimeter")
    print("per second' (molec cm^-3 s^-1).\n")
    print("To make this rate physically meaningful for mass budget analyses, we need")
    print("to convert it from a volumetric molecular rate to a mass flux in 'kilograms")
    print("per second' (kg/s) for each model grid cell. This requires knowledge of")
    print("the molar mass of mercury, Avogadro's number, and the volume of the grid cell.\n")
    print("The conversion formula from 'molec/cm^3/s' to 'kg/s' for a grid cell is:\n")
    print("  $$\\text{Rate}_{\\text{kg/s}} = \\text{Rate}_{\\text{molec/cm}^3\\text{/s}} \\times \\left( \\frac{\\text{M}_{\\text{Hg}}}{\\text{N}_{\\text{A}}} \\right) \\times \\text{Volume}_{\\text{grid, cm}^3}$$")
    print("Where:")
    print(f"  - $M_{{Hg}}$: Molar Mass of Mercury ({MOLAR_MASS_HG_KG_PER_MOL} kg/mol)")
    print(f"  - $N_{{A}}$: Avogadro's Number ({AVOGADRO_NUMBER:.3e} molec/mol)")
    print("  - $\\text{Volume}_{\\text{grid, cm}^3}$: Volume of the grid cell in cm$^3$ (obtained from 'Met_AIRVOL').\n")

    hg2_to_ssa_transfer_rate_var_name = 'Hg2GasToSSA'
    transfer_rate_molec_cm3_s = None

    if hg2_to_ssa_transfer_rate_var_name in ds_mercurychem.data_vars:
        transfer_rate_molec_cm3_s = ds_mercurychem[hg2_to_ssa_transfer_rate_var_name]
        print(f"  - Found and extracted '{hg2_to_ssa_transfer_rate_var_name}' (Shape: {transfer_rate_molec_cm3_s.shape}, Units: {transfer_rate_molec_cm3_s.attrs.get('units', 'N/A')})")
    else:
        print(f"  - ERROR: Variable '{hg2_to_ssa_transfer_rate_var_name}' not found in MercuryChem dataset. Cannot quantify transfer.")
        return None

    grid_volume_m3 = None
    if 'Met_AIRVOL' in ds_statemet.data_vars:
        grid_volume_m3 = ds_statemet['Met_AIRVOL']
        print(f"\n  - Extracted 'Met_AIRVOL' from StateMet dataset (Shape: {grid_volume_m3.shape}, Units: {grid_volume_m3.attrs.get('units', 'N/A')}).")
        # Convert AIRVOL from m^3 to cm^3 for unit consistency
        grid_volume_cm3 = grid_volume_m3 / CM3_TO_M3
        print(f"    -> Converted 'Met_AIRVOL' from m^3 to cm^3 for unit consistency in the calculation.")
    else:
        print("  - CRITICAL ERROR: 'Met_AIRVOL' variable (grid cell volume) not found in StateMet dataset. Cannot convert transfer rate from molec/cm^3/s to kg/s.")
        return None
    
    # --- Step 2/3: Calculation of Transfer Rate (in kg/s) ---
    print("\n[STEP 2/3] Calculation: Converting Volumetric Rate to Mass Flux (kg/s)")
    print("----------------------------------------------------------------------")
    print("Here, we perform the unit conversion. The direct transfer rate from the")
    print("model ($Hg2GasToSSA$) in molec cm^-3 s^-1 is converted to a mass flux in")
    print("kg/s for each grid cell using the formula defined above. This provides")
    print("a more intuitive measure of the mass of mercury being transferred.\n")

    if transfer_rate_molec_cm3_s is not None and grid_volume_m3 is not None:
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
                summary_value = hg2_to_ssa_transfer_rate_kg_s.isel(time=0, lev=0).sum().compute().item()
                print(f"  - Example value (first time step, first level, sum over space): {summary_value:.4e} kg/s")
            elif 'time' in hg2_to_ssa_transfer_rate_kg_s.dims and hg2_to_ssa_transfer_rate_kg_s['time'].size > 0:
                summary_value = hg2_to_ssa_transfer_rate_kg_s.isel(time=0).sum().compute().item()
                print(f"  - Example value (first time step, sum over space): {summary_value:.4e} kg/s")
            else:
                summary_value = hg2_to_ssa_transfer_rate_kg_s.sum().compute().item()
                print(f"  - Example value (overall sum over space): {summary_value:.4e} kg/s")
            print("  - The calculated transfer rate result is an xarray.DataArray.")

        except Exception as e:
            print(f"  - ERROR: Calculation for Hg(II) to SSA transfer failed. This might be due to incompatible dimensions or coordinates: {e}")
            hg2_to_ssa_transfer_rate_kg_s = None # Ensure it's set to None on error
    else:
        print("  - Skipping Hg(II) to SSA transfer calculation due to missing primary variables.")

    print("\n" + "="*80)
    print("--- Hg(II) GAS TO SEA SALT AEROSOL TRANSFER QUANTIFICATION COMPLETE ---".center(80))
    print("="*80 + "\n")
    return hg2_to_ssa_transfer_rate_kg_s
