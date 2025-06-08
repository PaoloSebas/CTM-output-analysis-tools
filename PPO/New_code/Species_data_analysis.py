### main.py 
### 08062025 Paolo & Copilot in GitHub

import nc4_utils

def main():
    folder = input("Inserisci il percorso della cartella che contiene i file .nc4: ").strip()
    files = nc4_utils.list_nc4_files(folder)
    if not files:
        print("Nessun file .nc4 trovato nella cartella indicata.")
        return
    for file in files:
        nc_data = nc4_utils.load_nc4_file(file)
        print(f"File {file} caricato con successo!")
        # Qui puoi aggiungere altre operazioni sui dati

if __name__ == "__main__":
    main()
