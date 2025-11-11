from pathlib import Path
from typing import Optional, List
import xarray as xr
import logging
from .config import cfg
from .utils import safe_list_files, ensure_time_month

logger = logging.getLogger(__name__)

def build_file_pattern(year: str, type_name: str, data_type: str) -> str:
    """
    Costruisce il pattern per glob (stringa) come nel file originale.
    """
    return f"GEOSChem.{data_type}.{year}*.nc4"

def get_dataset(year: str, type_name: str, data_type: str, data_dir: Path = None,
                use_dask: bool = True, preprocess=None) -> Optional[xr.Dataset]:
    """
    Carica uno o più file xarray per anno/tipo/data_type.
    Ritorna xr.Dataset o None se non trova file o errori.
    """
    data_dir = data_dir or cfg.data_dir
    short_key = f"{cfg.perturbation_type_map.get(type_name, type_name)}_{year}_{cfg.data_type_map.get(data_type, data_type)}"
    logger.info(f"Ricerca file per {short_key} in {data_dir}/{year}/{type_name}")

    folder = Path(data_dir) / year / type_name
    pattern = build_file_pattern(year, type_name, data_type)
    files = safe_list_files(folder, pattern)

    if not files:
        logger.warning(f"get_dataset: nessun file trovato per {short_key}")
        return None

    try:
        open_kwargs = dict(combine='by_coords')
        if use_dask:
            # non specifico chunk qui; l'utente può modificare se necessario
            open_kwargs['parallel'] = True
        ds = xr.open_mfdataset([str(p) for p in files], **open_kwargs)
        if preprocess is not None:
            ds = preprocess(ds)
        # aggiungo coord month utile dopo
        ds = ensure_time_month(ds)
        logger.info(f"Caricati {len(files)} file per {short_key}")
        return ds
    except Exception as e:
        logger.exception(f"Errore caricamento dataset {short_key}: {e}")
        return None
