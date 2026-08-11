import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compute_feature_contributions(model_dict: dict, candidate_df: pd.DataFrame):
    if not model_dict.get("calibrated"):
        return None
    
    model = model_dict["model"]
    X = candidate_df[['vina_score', 'gnina_score']].values
    preds = model.predict(X)
    
    vina_contrib = X[:, 0] * model.coef_[0]
    gnina_contrib = X[:, 1] * model.coef_[1]
    
    results = candidate_df.copy()
    results['predicted_pIC50'] = preds
    results['predicted_IC50_uM'] = (10 ** (-preds)) * 1e6
    results['vina_impact'] = vina_contrib
    results['gnina_impact'] = gnina_contrib
    
    fig, ax = plt.subplots(figsize=(6, 3))
    top_row = results.iloc[0]
    features = ['AutoDock Vina', 'GNINA ML Score']
    impacts = [top_row['vina_impact'], top_row['gnina_impact']]
    
    colors = ['#1f77b4' if v >= 0 else '#d62728' for v in impacts]
    ax.barh(features, impacts, color=colors)
    ax.set_xlabel("Contribution to Predicted pIC50")
    ax.set_title(f"Feature Breakdown: {top_row['compound_name']}")
    ax.grid(True, linestyle="--", alpha=0.3)
    plt.tight_layout()
    
    return {
        "ranked_df": results.sort_values("predicted_pIC50", ascending=False),
        "impact_plot": fig
    }
