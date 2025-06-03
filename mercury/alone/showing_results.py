# showing_results.py
import xarray as xr
from tabulate import tabulate
import numpy as np

# Import plotting functions
from plotting import plot_variable_global, plot_variable_antarctic

def display_mercury_results(wet_dep_result: xr.DataArray = None,
                            dry_dep_result: xr.DataArray = None,
                            hg2_ssa_transfer_result: xr.DataArray = None):
    """
    Displays the textual summary of the mercury analysis results.
    It prints aggregated values for total wet/dry deposition or Hg(II) to SSA transfer.

    Args:
        wet_dep_result (xr.DataArray, optional): Total wet deposition DataArray.
        dry_dep_result (xr.DataArray, optional): Total dry deposition DataArray.
        hg2_ssa_transfer_result (xr.DataArray, optional): Hg(II) to SSA transfer DataArray.
    """
    print("\n--- Summary of Mercury Analysis Results ---")

    results_table_data = []
    headers = ["Analysis Type", "Time Period", "Total (kg/s)"]

    if dry_dep_result is not None and wet_dep_result is not None:
        print("\nDidactic Note: For multi-time step data, we will display the spatially-integrated sum")
        print("for each time step. For single time steps, we will display the overall sum.")
        print("These values represent the total mass of mercury deposited/transferred across the entire domain.")

        # Handle Total Dry Deposition
        if 'time' in dry_dep_result.dims:
            for time_idx, time_val in enumerate(dry_dep_result['time'].values):
                # Sum over all non-time dimensions
                spatial_sum = dry_dep_result.isel(time=time_idx).sum().item()
                results_table_data.append([
                    "Total Dry Deposition",
                    f"Time Step {time_idx+1} ({time_val})",
                    f"{spatial_sum:.4e}"
                ])
        else: # No time dimension
            spatial_sum = dry_dep_result.sum().item()
            results_table_data.append([
                "Total Dry Deposition",
                "Overall",
                f"{spatial_sum:.4e}"
            ])

        # Handle Total Wet Deposition
        if 'time' in wet_dep_result.dims:
            for time_idx, time_val in enumerate(wet_dep_result['time'].values):
                # Sum over all non-time dimensions
                spatial_sum = wet_dep_result.isel(time=time_idx).sum().item()
                results_table_data.append([
                    "Total Wet Deposition",
                    f"Time Step {time_idx+1} ({time_val})",
                    f"{spatial_sum:.4e}"
                ])
        else: # No time dimension
            spatial_sum = wet_dep_result.sum().item()
            results_table_data.append([
                "Total Wet Deposition",
                "Overall",
                f"{spatial_sum:.4e}"
            ])

    elif hg2_ssa_transfer_result is not None:
        print("\nDidactic Note: For multi-time step/level data, we will display the spatially-integrated sum")
        print("for each time step and level. For single time steps, we will display the overall sum.")
        print("These values represent the total mass of mercury transferred across the entire domain for each instance.")

        # Handle Hg(II) to SSA Transfer
        if 'time' in hg2_ssa_transfer_result.dims and 'lev' in hg2_ssa_transfer_result.dims:
            for time_idx, time_val in enumerate(hg2_ssa_transfer_result['time'].values):
                for lev_idx, lev_val in enumerate(hg2_ssa_transfer_result['lev'].values):
                    spatial_sum = hg2_ssa_transfer_result.isel(time=time_idx, lev=lev_idx).sum().item()
                    results_table_data.append([
                        "Hg(II) to SSA Transfer",
                        f"Time {time_idx+1} (Lev {lev_val})",
                        f"{spatial_sum:.4e}"
                    ])
        elif 'time' in hg2_ssa_transfer_result.dims: # Only time dimension
            for time_idx, time_val in enumerate(hg2_ssa_transfer_result['time'].values):
                spatial_sum = hg2_ssa_transfer_result.isel(time=time_idx).sum().item()
                results_table_data.append([
                    "Hg(II) to SSA Transfer",
                    f"Time Step {time_idx+1} ({time_val})",
                    f"{spatial_sum:.4e}"
                ])
        else: # No time or level dimension
            spatial_sum = hg2_ssa_transfer_result.sum().item()
            results_table_data.append([
                "Hg(II) to SSA Transfer",
                "Overall",
                f"{spatial_sum:.4e}"
            ])
    else:
        print("No valid results found to display.")
        return

    print(tabulate(results_table_data, headers=headers, tablefmt="grid"))
    print("-------------------------------------------\n")


def offer_and_execute_plotting(wet_dep_result: xr.DataArray = None,
                               dry_dep_result: xr.DataArray = None,
                               hg2_ssa_transfer_result: xr.DataArray = None):
    """
    Offers plotting options based on the available analysis results and
    executes the chosen plots.

    Args:
        wet_dep_result (xr.DataArray, optional): Total wet deposition DataArray.
        dry_dep_result (xr.DataArray, optional): Total dry deposition DataArray.
        hg2_ssa_transfer_result (xr.DataArray, optional): Hg(II) to SSA transfer DataArray.
    """
    plot_options = []
    available_data_for_plot = {}
    plot_counter = 1

    print("\n--- Plotting Options ---")
    print("Didactic Note: For variables with multiple time steps or vertical levels,")
    print("we will plot the *first time step* and *first level* (if applicable)")
    print("to provide a representative spatial distribution.\n")

    if dry_dep_result is not None and wet_dep_result is not None:
        # For multi-dimensional results, select the first time step and first level
        # if 'time' dimension exists, slice to a single time point for 2D plotting.
        # If 'lev' dimension exists, sum or select a level for 2D plotting.
        if 'time' in dry_dep_result.dims:
            dry_dep_2d = dry_dep_result.isel(time=0)
            print(f"  - Preparing 'Total Dry Deposition' for plotting (using first time step: {dry_dep_result['time'].values[0]}).")
        else:
            dry_dep_2d = dry_dep_result
            print("  - Preparing 'Total Dry Deposition' for plotting (no time dimension).")

        if 'lev' in dry_dep_2d.dims: # If it's still 3D (time,lev,lat,lon), e.g. after time slicing
            dry_dep_2d = dry_dep_2d.isel(lev=0) # Select first level
            print(f"  - Further preparing 'Total Dry Deposition' for plotting (using first level: {dry_dep_result['lev'].values[0] if 'lev' in dry_dep_result.dims else 'N/A'}).")

        plot_options.append(f"[{plot_counter}] Total Dry Deposition (Global)")
        available_data_for_plot[str(plot_counter)] = (dry_dep_2d, "Total Dry Deposition", plot_variable_global)
        plot_counter += 1

        if 'time' in wet_dep_result.dims:
            wet_dep_2d = wet_dep_result.isel(time=0)
            print(f"  - Preparing 'Total Wet Deposition' for plotting (using first time step: {wet_dep_result['time'].values[0]}).")
        else:
            wet_dep_2d = wet_dep_result
            print("  - Preparing 'Total Wet Deposition' for plotting (no time dimension).")

        if 'lev' in wet_dep_2d.dims:
            wet_dep_2d = wet_dep_2d.isel(lev=0) # Select first level
            print(f"  - Further preparing 'Total Wet Deposition' for plotting (using first level: {wet_dep_result['lev'].values[0] if 'lev' in wet_dep_result.dims else 'N/A'}).")

        plot_options.append(f"[{plot_counter}] Total Wet Deposition (Global)")
        available_data_for_plot[str(plot_counter)] = (wet_dep_2d, "Total Wet Deposition", plot_variable_global)
        plot_counter += 1

        plot_options.append(f"[{plot_counter}] Total Dry Deposition (Antarctic Focus)")
        available_data_for_plot[str(plot_counter)] = (dry_dep_2d, "Total Dry Deposition", plot_variable_antarctic)
        plot_counter += 1

        plot_options.append(f"[{plot_counter}] Total Wet Deposition (Antarctic Focus)")
        available_data_for_plot[str(plot_counter)] = (wet_dep_2d, "Total Wet Deposition", plot_variable_antarctic)
        plot_counter += 1


    elif hg2_ssa_transfer_result is not None:
        if 'time' in hg2_ssa_transfer_result.dims:
            hg2_ssa_2d = hg2_ssa_transfer_result.isel(time=0)
            print(f"  - Preparing 'Hg(II) to SSA Transfer' for plotting (using first time step: {hg2_ssa_transfer_result['time'].values[0]}).")
        else:
            hg2_ssa_2d = hg2_ssa_transfer_result
            print("  - Preparing 'Hg(II) to SSA Transfer' for plotting (no time dimension).")

        if 'lev' in hg2_ssa_2d.dims:
            hg2_ssa_2d = hg2_ssa_2d.isel(lev=0) # Select first level
            print(f"  - Further preparing 'Hg(II) to SSA Transfer' for plotting (using first level: {hg2_ssa_transfer_result['lev'].values[0] if 'lev' in hg2_ssa_transfer_result.dims else 'N/A'}).")

        plot_options.append(f"[{plot_counter}] Hg(II) Gas to Sea Salt Aerosol Transfer (Global)")
        available_data_for_plot[str(plot_counter)] = (hg2_ssa_2d, "Hg(II) to SSA Transfer", plot_variable_global)
        plot_counter += 1

        plot_options.append(f"[{plot_counter}] Hg(II) Gas to Sea Salt Aerosol Transfer (Antarctic Focus)")
        available_data_for_plot[str(plot_counter)] = (hg2_ssa_2d, "Hg(II) to SSA Transfer", plot_variable_antarctic)
        plot_counter += 1

    if not plot_options:
        print("No valid results available to plot.")
        return

    while True:
        print("\nAvailable Plots:")
        for option_str in plot_options:
            print(f"  {option_str}")
        print(f"  [{plot_counter}] Plot All Available")
        print(f"  [0] Go back / No more plots")

        plot_choice = input("Enter plot choice (e.g., 1, 2, or '0' to quit plotting): ").strip()

        if plot_choice == '0':
            print("Exiting plotting menu.")
            break
        elif plot_choice == str(plot_counter): # Plot All
            print("Plotting all available results...")
            for key in available_data_for_plot:
                data, name, plot_func = available_data_for_plot[key]
                print(f"Generating plot for: {name} with {plot_func.__name__}...")
                plot_func(data, name)
        elif plot_choice in available_data_for_plot:
            data, name, plot_func = available_data_for_plot[plot_choice]
            print(f"Generating plot for: {name} with {plot_func.__name__}...")
            plot_func(data, name)
        else:
            print("Invalid plot choice. Please try again.")
