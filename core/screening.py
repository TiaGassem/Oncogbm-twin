import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_squared_error, r2_score

def calibrate_anchor_model(df_anchor: pd.DataFrame):
    df = df_anchor.copy()
    df['pIC50'] = -np.log10(df['wetlab_ic50_uM'] * 1e-6)
    
    X = df[['vina_score', 'gnina_score']].values
    y = df['pIC50'].values
    
    n_samples = len(df)
    if n_samples < 4:
        return {
            "calibrated": False,
            "warning": f"Anchor set size (N={n_samples}) is too small. Showing raw docking rankings."
        }
    
    loo = LeaveOneOut()
    y_true, y_pred = [], []
    
    for train_idx, test_idx in loo.split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        
        model = Ridge(alpha=1.0)
        model.fit(X_tr, y_tr)
        pred = model.predict(X_te)
        
        y_true.append(y_te[0])
        y_pred.append(pred[0])
        
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    
    final_model = Ridge(alpha=1.0)
    final_model.fit(X, y)
    
    return {
        "calibrated": True,
        "n_anchors": n_samples,
        "loocv_rmse": rmse,
        "loocv_r2": r2,
        "model": final_model,
        "coefficients": {
            "vina_weight": final_model.coef_[0],
            "gnina_weight": final_model.coef_[1],
            "intercept": final_model.intercept_
        }
    }
