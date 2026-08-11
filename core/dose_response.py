import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def four_parameter_logistic(x, a, b, c, d):
    return d + (a - d) / (1.0 + (np.maximum(x, 1e-12) / c) ** b)

def fit_4pl_curve(concentrations_uM: list, viability_pct: list):
    x = np.array(concentrations_uM, dtype=float)
    y = np.array(viability_pct, dtype=float)
    
    if len(x) < 4:
        return {"success": False, "error": "Need at least 4 concentration points for 4PL fit."}
    
    p0 = [min(y), 1.0, np.median(x), max(y)]
    bounds = ([0.0, 0.1, 1e-6, 0.0], [100.0, 10.0, max(x) * 10, 150.0])
    
    try:
        popt, _ = curve_fit(four_parameter_logistic, x, y, p0=p0, bounds=bounds, maxfev=10000)
        a, b, c, d = popt
        
        fig, ax = plt.subplots(figsize=(6, 4))
        x_dense = np.logspace(np.log10(min(x) * 0.5), np.log10(max(x) * 2), 200)
        y_dense = four_parameter_logistic(x_dense, a, b, c, d)
        
        ax.scatter(x, y, color="#1f77b4", label="Experimental Data", zorder=3)
        ax.plot(x_dense, y_dense, color="#d62728", linestyle="--", label=f"4PL Fit (IC50 = {c:.3f} µM)")
        ax.axhline(50, color="gray", linestyle=":", alpha=0.7)
        ax.set_xscale("log")
        ax.set_xlabel("Concentration (µM)")
        ax.set_ylabel("Cell Viability (%)")
        ax.set_title("4PL Dose-Response Curve")
        ax.legend()
        ax.grid(True, which="both", alpha=0.2)
        plt.tight_layout()
        
        return {
            "success": True,
            "ic50_uM": c,
            "hill_slope": b,
            "min_viability": a,
            "max_viability": d,
            "figure": fig
        }
    except Exception as e:
        return {"success": False, "error": str(e)}
