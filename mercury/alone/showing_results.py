# showing_results.py

import xarray as xr

def display_mercury_results(total_dry_deposition_kg_s: xr.DataArray | None,
                            total_wet_deposition_kg_s: xr.DataArray | None,
                            hg2_to_ssa_transfer_rate_kg_s: xr.DataArray | None):
    """
    Displays the results of the mercury deposition and transfer analysis.

    Args:
        total_dry_deposition_kg_s (xr.DataArray | None): Total dry deposition in kg/s.
        total_wet_deposition_kg_s (xr.DataArray | None): Total wet deposition in kg/s.
        hg2_to_ssa_transfer_rate_kg_s (xr.DataArray | None): Hg(II) to SSA transfer rate in kg/s.
    """
    if (total_dry_deposition_kg_s is None and
        total_wet_deposition_kg_s is None and
        hg2_to_ssa_transfer_rate_kg_s is None):
        print("No results to display. Analysis may have failed or produced no valid data.")
        return

    print("\n" + "="*80)
    print("                      --- MERCURY ANALYSIS SUMMARY ---                      ".center(80))
    print("="*80)

    # Helper function to get spatial sum and then time average
    def calculate_summary_rate(data_array: xr.DataArray, name: str):
        if data_array is None:
            print(f"\n--- {name}: Calculation Failed or Data Not Available ---")
            return

        print(f"\n--- {name} ---")
        
        # Determine spatial dimensions to sum over (lat, lon, and potentially lev)
        # 'lev' is treated as a spatial dimension for summation if present, to get total column rate
        spatial_dims = [dim for dim in data_array.dims if dim not in ['time']]
        
        try:
            # Step 1: Sum over all spatial dimensions (lat, lon, and 'lev' if present)
            # This results in a DataArray with only the 'time' dimension (if it exists)
            # and represents the total rate for the entire domain (and vertical column) at each time step.
            spatial_sum_da = data_array.sum(dim=spatial_dims)

            # Step 2: Average over the time dimension if it exists and has data
            avg_rate = None
            if 'time' in spatial_sum_da.dims and spatial_sum_da['time'].size > 0:
                avg_rate = spatial_sum_da.mean(dim='time').compute().item()
                print(f"  - Average Total {name} Rate (over simulation period): {avg_rate:.4e} kg/s")
            else:
                # If no time dimension, it's a single snapshot or climatology.
                # The spatial_sum_da is already a scalar or has no 'time' dim left.
                avg_rate = spatial_sum_da.compute().item()
                print(f"  - Total {name} Rate (single snapshot, summed over space): {avg_rate:.4e} kg/s")

            print(f"  - DataArray shape: {data_array.shape}")
            print(f"  - DataArray units: {data_array.attrs.get('units', 'N/A')}")

        except Exception as e:
            print(f"  - ERROR: Could not compute summary for {name}: {e}")
            print(f"  - DataArray shape: {data_array.shape}")
            print(f"  - DataArray units: {data_array.attrs.get('units', 'N/A')}")

    # Display results using the helper function
    calculate_summary_rate(total_dry_deposition_kg_s, "Dry Deposition")
    calculate_summary_rate(total_wet_deposition_kg_s, "Wet Deposition")
    calculate_summary_rate(hg2_to_ssa_transfer_rate_kg_s, "Hg(II) Gas to Sea Salt Aerosol Transfer")

    print("\n" + "="*80)
    print("                          --- SUMMARY COMPLETE ---                          ".center(80))
    print("="*80)
