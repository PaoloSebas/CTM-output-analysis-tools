import xarray as xr
import matplotlib.pyplot as plt
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean

# Read data
gcc_ds = xr.open_dataset(
    '/mnt/d/PPO/NEW_SET_2025/SCENARIO_1/GEOSChem.Restart.fullchem.20190701_0000z.nc4'
)
modified_ds = xr.open_dataset(
    '/mnt/d/PPO/NEW_SET_2025/SCENARIO_1/SCENARIO1_GEOSChem.StateMet.20190701_0000z.nc4'
)

# Plot comparison of surface ozone over the North Pacific
compare_single_level(
    gcc_ds,
    'GEOS-Chem - Restart file',
    modified_ds,
    'Spin-up',
    varlist=['SpeciesConcVV_HO2'],
    extra_title_txt='Surface'
)
plt.show()
