

import os
from netCDF4 import Dataset

def list_nc4_files(folder_path):
    """Restituisce una lista di path dei file .nc4 nella cartella."""
    return [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.nc4')]

def load_nc4_file(file_path):
    """Apre un file .nc4 e restituisce l’oggetto Dataset."""
    return Dataset(file_path, mode='r')

# Aggiungeremo altre funzioni secondo necessità!
