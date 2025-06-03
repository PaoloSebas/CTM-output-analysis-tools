#IMPORT necessary Libraries

import xarray as xr
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import gcpy.plot as gcplot
from gcpy.plot.single_panel import single_panel
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean
from matplotlib.backends.backend_pdf import PdfPages

# Read data from source nc4 files

#COMPLETE DATASET

ds_i = xr.open_dataset(
    '/mnt/d/PPO/NEW_SET_2025/SCENARIO_1/SCENARIO1_GEOSChem.SpeciesConc.20190701_0000z.nc4'
)

print('This is the dataset I want investigate:', ds_i)
print('####################################################################################')


#Defining the species I want to plot

Species_X='HO2'
Species_X_units = 'SpeciesConcVV_HO2'
X = ds_i.variables[Species_X_units]
print(X)

print('########################## SUBSETTING ON Species:', Species_X, ': ##################################')
Xds_i = ds_i['SpeciesConcVV_HO2']

print('```````````````````````````````````````````````````````````````````')
print('This is the dataset after SUBSETTING on HO2', Xds_i)
print('```````````````````````````````````````````````````````````````````')

# SUM only over time
Xdsi_time = (Xds_i.sum(dim=['time'], keep_attrs=True))
Xdsi_time_ppb = Xdsi_time*1e+9
Xdsi_burden = np.sum(Xdsi_time_ppb.values)

print('Burden ppb in one month:', Xdsi_burden)

#diff_time = ((Xdsf_time - Xdsi_time)/(Xdsi_time))*100
max_plot = np.max(Xdsi_time_ppb.values)
latitudes = ds_i.variables['lat']
longitudes = ds_i.variables['lon']
levels = ds_i.variables['lev']

print(Xdsi_time_ppb.shape)

max_index = np.unravel_index(np.argmax(Xdsi_time_ppb.values), Xdsi_time_ppb.shape)
max_lev = levels[max_index[0]]
max_latitude = latitudes[max_index[1]]
max_longitude = longitudes[max_index[2]]
print('Max lat, max long', max_lev, max_latitude, max_longitude )


min_plot = np.min(Xdsi_time_ppb.values)

print('Maximum value in plot(ppb):',max_plot)
print('Minimum value in plot(ppb):',min_plot)

def max_modulus(a, b):
    if abs(a) > abs(b):
        return a
    else:
        return b

max_plot_x = max_modulus(max_plot,min_plot)
print('max_plot_x:', max_plot_x)

#print(f"Sum  of diff_time on all values: {np.sum(diff_time.values)}")

abs_max=abs(max_plot_x)
print('The scale max value in the plot is:',abs_max)


single_panel(
  Xdsi_time_ppb.isel(lev=0),
  title=Species_X + ' (ppb) - GEOS-Chem v14.3.0 - no modifications 07/2019 - LEV=0',
  #comap=plt.get_cmap("Purples"),
  #log_color_scale=False,
    #vmin= -abs_max,
    #vmax= abs_max,
     vmin= 0,
     vmax= 0.016,
   extent=[-180, 180, -90, 90]
)

plt.show()
                                                      
