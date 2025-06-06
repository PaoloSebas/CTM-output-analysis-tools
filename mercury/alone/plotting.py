# plotting.py
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker 
import numpy as np 
from matplotlib.path import Path 

def plot_variable_global(data_array: xr.DataArray, var_name: str = None):
    """
    Plots a 2D global distribution of an xarray.DataArray.
    Assumes time/level slicing has already been performed if needed.

    Args:
        data_array (xr.DataArray): The 2D DataArray (lat, lon) to plot.
        var_name (str, optional): The name of the variable for the plot title and colorbar.
                                   Defaults to data_array.name if not provided.
    """
    if 'lon' not in data_array.dims or 'lat' not in data_array.dims:
        print(f"Error: DataArray '{data_array.name}' must have 'lon' and 'lat' dimensions for global plotting.")
        return

    # Use the name from the DataArray or the provided var_name for the title
    display_name = var_name if var_name else data_array.name
    units = data_array.attrs.get('units', 'unitless')
    long_name = data_array.attrs.get('long_name', display_name) # Use long_name if available

    print(f"\nGenerating Global Plot for '{display_name}'...")

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree()) # Standard global projection
    
    # Plotting the data using pcolormesh for continuous color mapping
    # Infer vmin/vmax if not explicitly set, or let xarray handle it
    plot = data_array.plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(), # Data coordinates
        x='lon', y='lat',
        add_colorbar=True,
        cbar_kwargs={'label': f'{long_name} ({units})'},
        cmap='viridis', # Choose an appropriate colormap (e.g., 'viridis', 'jet', 'cmo.balance')
        extend='both' # Extends colorbar for values outside vmin/vmax
    )

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.7)
    ax.add_feature(cfeature.LAKES, alpha=0.5, facecolor='lightcyan')
    ax.add_feature(cfeature.RIVERS, alpha=0.5, edgecolor='blue')
    ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    
    ax.set_title(f"Global Distribution of {display_name}")
    
    # Add gridlines with labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = cticker.LongitudeFormatter()
    gl.yformatter = cticker.LatitudeFormatter() # <--- CORRECTED LINE (removed .ticker)

    plt.tight_layout() # Adjust layout to prevent labels overlapping
    plt.show()

def plot_variable_antarctic(data_array: xr.DataArray, var_name: str = None):
    """
    Plots a 2D distribution focused on the Antarctic region using South Polar Stereographic projection.
    Assumes time/level slicing has already been performed if needed.

    Args:
        data_array (xr.DataArray): The 2D DataArray (lat, lon) to plot.
        var_name (str, optional): The name of the variable for the plot title and colorbar.
                                   Defaults to data_array.name if not provided.
    """
    if 'lon' not in data_array.dims or 'lat' not in data_array.dims:
        print(f"Error: DataArray '{data_array.name}' must have 'lon' and 'lat' dimensions for Antarctic plotting.")
        return

    display_name = var_name if var_name else data_array.name
    units = data_array.attrs.get('units', 'unitless')
    long_name = data_array.attrs.get('long_name', display_name)

    print(f"\nGenerating Antarctic Plot for '{display_name}'...")

    fig = plt.figure(figsize=(10, 10))
    # Use a stereographic projection centered on the South Pole
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo()) 
    
    # Set the extent to focus on the Antarctic region
    # lon_min, lon_max, lat_min, lat_max
    ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree()) # Adjusted slightly to show more context

    # Plotting the data
    plot = data_array.plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(), # Data coordinates are still PlateCarree
        x='lon', y='lat',
        add_colorbar=True,
        cbar_kwargs={'label': f'{long_name} ({units})'},
        cmap='viridis',
        extend='both'
    )

    ax.coastlines()
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')
    ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.7)

    ax.set_title(f"Antarctic Distribution of {display_name}")
    
    # Add gridlines and labels specific to polar projection
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False 
    gl.right_labels = False
    gl.xlabels_top = False
    gl.xlabels_bottom = True
    gl.ylabels_left = True
    gl.ylabels_right = False
    
    gl.xformatter = cticker.LongitudeFormatter()
    gl.yformatter = cticker.LatitudeFormatter() # <--- CORRECTED LINE (removed .ticker)
    
    # Create a circular boundary path in normalized axes coordinates (0-1)
    circle_center = (0.5, 0.5) # Center of the subplot in axes coordinates
    circle_radius = 0.5       # Radius to fill the subplot (from 0 to 1)

    # Generate points for a circle
    theta = np.linspace(0, 2*np.pi, 100) # 100 points to draw the circle
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    # Scale and translate these points to the desired center and radius
    circle_path = Path(verts * circle_radius + circle_center)

    # Set this path as the boundary, transforming it with ax.transAxes
    ax.set_boundary(circle_path, transform=ax.transAxes)
    
    plt.tight_layout()
    plt.show()
