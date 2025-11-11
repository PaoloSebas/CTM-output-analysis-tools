from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

@dataclass
class Config:
    # Percorso base dati (modifica al tuo path locale)
    data_dir: Path = Path("D:/MERCURY_ANU/DATA_ANALYSIS/v14_6_0/Preind_outputs/New_diagnostic_island/")

    # anni e tipi di perturbazione (modificabili)
    years: List[str] = field(default_factory=lambda: ["2018"])
    perturbation_types: List[str] = field(default_factory=lambda: [
        "NO_PERTUBS", "OCEAN_LESS", "OCEAN_MORE",
        "SSA_CONSTANT", "SSA_ZERO", "WIND_LESS", "WIND_MORE", "OX_ALL"
    ])

    # mappe corte per chiavi
    perturbation_type_map: Dict[str, str] = field(default_factory=lambda: {
        "NO_PERTUBS": "NP",
        "OCEAN_LESS": "OL",
        "OCEAN_MORE": "OM",
        "SSA_CONSTANT": "SSAC",
        "SSA_ZERO": "SSA0",
        "WIND_LESS": "WL",
        "WIND_MORE": "WM",
        "OX_ALL": "OX"
    })

    data_type_map: Dict[str, str] = field(default_factory=lambda: {
        "DryDep": "DD",
        "WetLossLS": "WLLS",
        "WetLossConv": "WLConv",
        "MercuryChem": "HgChem",
        "MercuryEmis": "HgEmis",
        "MercuryOcean": "HgOcean",
        "SpeciesConc": "SpCon",
        "StateMet": "Met",
        "Budget": "Budget"
    })

    # costanti e fattori conversione (usati in seguito)
    MW_Hg: float = 200.59
    NA: float = 6.023e23
    s_in_month: float = 2.628e6
    s_in_yr: float = 3.154e7
    ug_g: float = 1e6
    cm2_m2: float = 1e4

    @property
    def cf_units_dd_year(self) -> float:
        # from molecule/cm2*s to micrograms/m2*year
        return (self.MW_Hg / self.NA) * self.ug_g * self.cm2_m2 * self.s_in_yr

# istanza di default importabile
cfg = Config()
