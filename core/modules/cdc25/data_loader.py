import pandas as pd

def load_cdc25_anchor_set() -> pd.DataFrame:
    data = [
        {"compound_name": "NSC95397 (Lead)", "vina_score": -8.4, "gnina_score": 0.82, "wetlab_ic50_uM": 0.22},
        {"compound_name": "BN82002",          "vina_score": -7.1, "gnina_score": 0.65, "wetlab_ic50_uM": 2.40},
        {"compound_name": "Compound 5",        "vina_score": -6.8, "gnina_score": 0.58, "wetlab_ic50_uM": 5.10},
        {"compound_name": "IRC-083864",       "vina_score": -8.1, "gnina_score": 0.76, "wetlab_ic50_uM": 0.85},
        {"compound_name": "DA-30038",         "vina_score": -6.3, "gnina_score": 0.49, "wetlab_ic50_uM": 12.50}
    ]
    return pd.DataFrame(data)

def load_cdc25_screening_candidates() -> pd.DataFrame:
    data = [
        {"compound_name": "Novel_CDC25_Inh_01", "vina_score": -8.8, "gnina_score": 0.85},
        {"compound_name": "Novel_CDC25_Inh_02", "vina_score": -7.9, "gnina_score": 0.71},
        {"compound_name": "Novel_CDC25_Inh_03", "vina_score": -6.5, "gnina_score": 0.52},
        {"compound_name": "Novel_CDC25_Inh_04", "vina_score": -8.2, "gnina_score": 0.79}
    ]
    return pd.DataFrame(data)
