from typing import Tuple
import xarray as xr
import numpy as np
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def nearest_point(da: xr.DataArray, lat: float, lon: float) -> xr.DataArray:
    """
    Seleziona il punto più vicino (lat, lon) da una DataArray con coordinate 'lat' e 'lon'.
    Restituisce la DataArray selezionata (ancora con eventuale dimensione 'month' o 'time').
    """
    if 'lat' not in da.coords or 'lon' not in da.coords:
        raise ValueError("La DataArray non ha coordinate 'lat' e 'lon'.")
    return da.sel(lat=lat, lon=lon, method='nearest')

def safe_list_files(path: Path, pattern: str):
    """
    Safe wrapper per path.glob(pattern), ritorna lista ordinata di Path.
    """
    files = sorted(path.glob(pattern))
    if not files:
        logger.debug(f"safe_list_files: nessun file trovato per {path}/{pattern}")
    return files

def ensure_time_month(ds: xr.Dataset) -> xr.Dataset:
    """
    Se il dataset ha 'time' ma non 'month', aggiunge coordinate month per comodità.
    (utile quando si vuole groupby('time.month') in seguito)
    """
    if 'time' in ds.coords and 'month' not in ds.coords:
        try:
            ds = ds.assign_coords(month=ds['time'].dt.month)
        except Exception:
            logger.debug("ensure_time_month: impossibile creare coords 'month' (time non datelike)")
    return ds
