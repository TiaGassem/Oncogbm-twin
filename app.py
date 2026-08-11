import json
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_squared_error, r2_score
import streamlit as st

# ==========================================
# 1. PAGE CONFIGURATION & DISCLAIMER
# ==========================================
st.set_page_config(page_title="GBM-Twin Platform", layout="wide", page_icon="🧬")

st.warning(
    "**RESEARCH & EDUCATIONAL USE ONLY.** GBM-Twin is a computational platform designed for target prioritization and drug discovery research. It is **not** a clinical decision-support tool."
)

st.title("🧬 GBM-Twin: Computational Oncology Platform for Glioblastoma")
st.caption("A modular, free-to-use platform combining public gene exploration with anchor-calibrated small-molecule screening.")

# ==========================================
# 2. DATA DEEPLINKS & PUBLIC TARGET DICTIONARY
# ==========================================
GBM_TARGETS = {
    "CDC25A": {"uniprot": "P30304", "gene": "CDC25A", "type": "Cell Cycle Control"},
    "CDC25B": {"uniprot": "P30305", "gene": "CDC25B", "type": "Cell Cycle Control"},
    "CDC25C": {"uniprot": "P30307", "gene": "CDC25C", "type": "Cell Cycle Control"},
    "EGFR":   {"uniprot": "P00533", "gene": "EGFR",   "type": "Receptor Tyrosine Kinase"},
    "PDGFRA": {"uniprot": "P16234", "gene": "PDGFRA", "type": "Receptor Tyrosine Kinase"},
    "PTEN":   {"uniprot": "P60484", "gene": "PTEN",   "type": "Tumor Suppressor"},
    "TP53":   {"uniprot": "P04637", "gene": "TP53",   "type": "Tumor Suppressor"},
    "IDH1":   {"uniprot": "O75874", "gene": "IDH1",   "type": "Metabolic Enzyme"},
    "IDH2":   {"uniprot": "P48735", "gene": "IDH2",   "type": "Metabolic Enzyme"},
    "MGMT":   {"uniprot": "P16455", "gene": "MGMT",   "type": "DNA Repair"},
    "ATRX":   {"uniprot": "P46100", "gene": "ATRX",   "type": "Chromatin Remodeling"},
    "CDKN2A": {"uniprot": "Q8N726", "gene": "CDKN2A", "type": "Cell Cycle Inhibitor"},
    "CDKN2B": {"uniprot": "P42773", "gene": "CDKN2B", "type": "Cell Cycle Inhibitor"},
    "RB1":    {"uniprot": "P06400", "gene": "RB1",    "type": "Tumor Suppressor"},
    "NF1":    {"uniprot": "P21359", "gene": "NF1",    "type": "Ras GTPase Activator"},
    "TERT":   {"uniprot": "O14746", "gene": "TERT",   "type": "Telomere Maintenance"},
    "PIK3CA": {"uniprot": "P42336", "gene": "PIK3CA", "type": "Kinase Signalling"},
    "CDK4":   {"uniprot": "P11802", "gene": "CDK4",   "type": "Cell Cycle Kinase"},
    "CDK6":   {"uniprot": "Q00534", "gene": "CDK6",   "type": "Cell Cycle Kinase"},
    "MDM2":   {"uniprot": "Q00987", "gene": "MDM2",   "type": "p53 Regulator"},
    "MDM4":   {"uniprot": "O15151", "gene": "MDM4",   "type": "p53 Regulator"}
}

# ==========================================
# 3. LIVE REST API FETCHERS
# ==========================================
@st.cache_data(ttl=86400)
def fetch_uniprot_summary(uniprot_id: str) -> dict:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return {"status": "error", "message": f"UniProt API HTTP {response.status_code}"}
        data = response.json()
        rec_name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
        seq = data.get("sequence", {}).get("value", "")
        return {
            "status": "success",
            "accession": uniprot_id,
            "full_name": rec_name,
            "length": len(seq),
            "sequence_preview": f"{seq[:30]}...{seq[-10:]}" if seq else "N/A"
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}

@st.cache_data(ttl=86400)
def fetch_cbioportal_mutation_stats(gene_symbol: str) -> dict:
    url = f"https://www.cbioportal.org/api/studies/gbm_tcga_pan_can_atlas_2018/genes/{gene_symbol}/mutations"
    headers = {"Accept": "application/json"}
    try:
        response = requests.get(url, headers=headers, timeout=10)
        if response.status_code != 200:
            return {"status": "error", "message": f"cBioPortal API HTTP {response.status_code}"}
        mutations = response.json()
        sample_muts = [
            f"{m.get('proteinChange', 'Variant')} ({m.get('mutationType', 'Unknown')})"
            for m in mutations[:5]
        ]
        return {
            "status": "success",
            "total_mutations": len(mutations),
            "sample_variants": sample_muts if sample_muts else ["No recurrent variants recorded"]
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==========================================
# 4. STATISTICAL & SCREENING ENGINES
# ==========================================
def four_parameter_logistic(x, a, b, c, d):
    return d + (a - d) / (1.0 + (np.maximum(x, 1e-12) / c) ** b)

def fit_4pl_curve(concentrations_uM: list, viability_pct: list):
    x = np.array(concentrations_uM, dtype=float)
    y = np.array(viability_pct, dtype=float)
    if len(x) < 4:
        return {"success": False, "error": "Need at least 4 points for 4PL fit."}
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
        return {"success": True, "ic50_uM": c, "hill_slope": b, "figure": fig}
    except Exception as e:
        return {"success": False, "error": str(e)}

def load_cdc25_anchor_set() -> pd.DataFrame:
    return pd.DataFrame([
        {"compound_name": "NSC95397 (Lead)", "vina_score": -8.4, "gnina_score": 0.82, "wetlab_ic50_uM": 0.22},
        {"compound_name": "BN82002",          "vina_score": -7.1, "gnina_score": 0.65, "wetlab_ic50_uM": 2.40},
        {"compound_name": "Compound 5",        "vina_score": -6.8, "gnina_score": 0.58, "wetlab_ic50_uM": 5.10},
        {"compound_name": "IRC-083864",       "vina_score": -8.1, "gnina_score": 0.76, "wetlab_ic50_uM": 0.85},
        {"compound_name": "DA-30038",         "vina_score": -6.3, "gnina_score": 0.49, "wetlab_ic50_uM": 12.50}
    ])

def load_cdc25_screening_candidates() -> pd.DataFrame:
    return pd.DataFrame([
        {"compound_name": "Novel_CDC25_Inh_01", "vina_score": -8.8, "gnina_score": 0.85},
        {"compound_name": "Novel_CDC25_Inh_02", "vina_score": -7.9, "gnina_score": 0.71},
        {"compound_name": "Novel_CDC25_Inh_03", "vina_score": -6.5, "gnina_score": 0.52},
        {"compound_name": "Novel_CDC25_Inh_04", "vina_score": -8.2, "gnina_score": 0.79}
    ])

def calibrate_anchor_model(df_anchor: pd.DataFrame):
    df = df_anchor.copy()
    df['pIC50'] = -np.log10(df['wetlab_ic50_uM'] * 1e-6)
    X = df[['vina_score', 'gnina_score']].values
    y = df['pIC50'].values
    n_samples = len(df)
    if n_samples < 4:
        return {"calibrated": False, "warning": f"Anchor set size (N={n_samples}) is too small."}
    loo = LeaveOneOut()
    y_true, y_pred = [], []
    for train_idx, test_idx in loo.split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        model = Ridge(alpha=1.0)
        model.fit(X_tr, y_tr)
        y_true.append(y_te[0])
        y_pred.append(model.predict(X_te)[0])
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    final_model = Ridge(alpha=1.0)
    final_model.fit(X, y)
    return {"calibrated": True, "n_anchors": n_samples, "loocv_rmse": rmse, "loocv_r2": r2, "model": final_model}

def compute_feature_contributions(model_dict: dict, candidate_df: pd.DataFrame):
    if not model_dict.get("calibrated"):
        return None
    model = model_dict["model"]
    X = candidate_df[['vina_score', 'gnina_score']].values
    preds = model.predict(X)
    results = candidate_df.copy()
    results['predicted_pIC50'] = preds
    results['predicted_IC50_uM'] = (10 ** (-preds)) * 1e6
    results['vina_impact'] = X[:, 0] * model.coef_[0]
    results['gnina_impact'] = X[:, 1] * model.coef_[1]
    
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
    return {"ranked_df": results.sort_values("predicted_pIC50", ascending=False), "impact_plot": fig}

# ==========================================
# 5. USER INTERFACE DASHBOARD
# ==========================================
st.sidebar.header(" Execution Mode")
mode = st.sidebar.radio("Select Workflow Mode:", ["Mode A: Guided / Manual (Free / Web UI)", "Mode B: Automated / API (Neurosnap)"])

tab1, tab2, tab3 = st.tabs(["📖 Platform Overview", "🌐 Layer 1: GBM Protein Explorer", "🎯 Layer 2: CDC25 Module"])

with tab1:
    st.markdown("### Two-Layer Architecture")
    col1, col2 = st.columns(2)
    with col1:
        st.info("**Layer 1 — Public Explorer:** Live, query-on-demand connections (UniProt, cBioPortal) across 21 GBM targets.")
    with col2:
        st.success("**Layer 2 — Calibrated Screening:** Target-specific anchor-bridge calibration (CDC25 Phosphatase Inhibitors).")

with tab2:
    st.markdown("### Layer 1: GBM Target Explorer")
    target_gene = st.selectbox("Select Target Gene:", list(GBM_TARGETS.keys()))
    meta = GBM_TARGETS[target_gene]
    st.markdown(f"**Category:** `{meta['type']}` | **UniProt ID:** `{meta['uniprot']}`")
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("UniProt Sequence Info")
        u_data = fetch_uniprot_summary(meta['uniprot'])
        st.json(u_data)
    with c2:
        st.subheader("TCGA Mutations (cBioPortal)")
        c_data = fetch_cbioportal_mutation_stats(meta['gene'])
        if c_data['status'] == 'success':
            st.metric("Total Mutation Records", c_data['total_mutations'])
            for v in c_data['sample_variants']:
                st.write(f"- {v}")
        else:
            st.error(c_data['message'])

with tab3:
    st.markdown("### Layer 2: CDC25 Phosphatase Discovery")
    subtab_a, subtab_b, subtab_c = st.tabs(["1. Wet-Lab 4PL Fitting", "2. Anchor Calibration", "3. Mode A Manual Input"])
    
    with subtab_a:
        st.subheader("Wet-Lab Viability Fitting (4PL)")
        c_str = st.text_input("Concentrations (µM):", "0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0")
        v_str = st.text_input("Cell Viability (%):", "98.2, 91.5, 78.4, 32.1, 12.8, 4.2, 1.1")
        if st.button("Fit 4PL Curve"):
            try:
                concs = [float(x.strip()) for x in c_str.split(",")]
                viabs = [float(x.strip()) for x in v_str.split(",")]
                res = fit_4pl_curve(concs, viabs)
                if res['success']:
                    st.success(f"Calculated IC50: {res['ic50_uM']:.4f} µM")
                    st.pyplot(res['figure'])
                else:
                    st.error(res['error'])
            except Exception as e:
                st.error(f"Error parsing data: {e}")
                
    with subtab_b:
        st.subheader("Anchor-Set Calibration")
        df_anchor = load_cdc25_anchor_set()
        st.dataframe(df_anchor, use_container_width=True)
        calib_res = calibrate_anchor_model(df_anchor)
        if calib_res['calibrated']:
            st.metric("LOOCV RMSE (pIC50)", f"{calib_res['loocv_rmse']:.3f}")
            df_cand = load_cdc25_screening_candidates()
            explain_res = compute_feature_contributions(calib_res, df_cand)
            st.dataframe(explain_res['ranked_df'][['compound_name', 'vina_score', 'gnina_score', 'predicted_IC50_uM']], use_container_width=True)
            st.pyplot(explain_res['impact_plot'])
            
    with subtab_c:
        st.subheader("Mode A: External Tool Parser")
        sample_json = json.dumps({"compound": "NSC95397_Derivative", "target": "CDC25A", "affinity_score": 0.89}, indent=2)
        user_json = st.text_area("Paste Tool Output JSON:", sample_json, height=180)
        if st.button("Parse Result"):
            try:
                st.json(json.loads(user_json))
            except Exception as e:
                st.error(f"Invalid JSON: {e}")
