import streamlit as st
import pandas as pd
import json

from core.explorer import GBM_TARGETS, fetch_uniprot_summary, fetch_cbioportal_mutation_stats
from core.dose_response import fit_4pl_curve
from core.screening import calibrate_anchor_model
from core.explainability import compute_feature_contributions
from modules.cdc25.data_loader import load_cdc25_anchor_set, load_cdc25_screening_candidates

st.set_page_config(page_title="GBM-Twin Platform", layout="wide", page_icon="🧬")

st.warning(
    "**RESEARCH & EDUCATIONAL USE ONLY.** GBM-Twin is a computational platform designed for target prioritization and drug discovery research. It is **not** a validated clinical decision-support tool."
)

st.title("🧬 GBM-Twin: Computational Oncology Platform for Glioblastoma")
st.caption("A modular, free-to-use platform combining public gene exploration with anchor-calibrated small-molecule screening.")

st.sidebar.header(" Execution Mode")
mode = st.sidebar.radio(
    "Select Workflow Mode:",
    ["Mode A: Guided / Manual (Free / Web UI)", "Mode B: Automated / API (Neurosnap)"]
)

if "Mode A" in mode:
    st.sidebar.info("💡 **Mode A Active:** Run docking/complexes on web UIs ( Neurosnap Chai-1/Boltz-2), then upload or paste JSON/CSV results manually below.")
else:
    st.sidebar.warning("⚡ **Mode B Active:** API integration mode. Requires configured Neurosnap/cBioPortal API keys.")

tab1, tab2, tab3 = st.tabs([
    "📖 Platform Overview", 
    "🌐 Layer 1: GBM Protein Explorer", 
    "🎯 Layer 2: CDC25 Inhibitor Module (Thesis)"
])

with tab1:
    st.markdown("### Two-Layer Architecture")
    col1, col2 = st.columns(2)
    with col1:
        st.info("""
        **Layer 1 — Public GBM Protein Explorer**
        * Live, query-on-demand connections (UniProt, cBioPortal, PDB).
        * Covers all 21 key glioblastoma target proteins.
        * Requires zero local database downloads or wet-lab data.
        """)
    with col2:
        st.success("""
        **Layer 2 — Calibrated Screening Modules**
        * Target-specific anchor-bridge calibration.
        * Calibrates molecular docking descriptors using real wet-lab IC50 values.
        * Currently active: **CDC25 Phosphatase Inhibitor Discovery** (Thesis Module).
        """)

with tab2:
    st.markdown("### Layer 1: GBM Target Explorer (Query-on-Demand)")
    target_gene = st.selectbox("Select GBM Target Gene:", list(GBM_TARGETS.keys()))
    meta = GBM_TARGETS[target_gene]
    
    st.markdown(f"**Category:** `{meta['type']}` | **UniProt ID:** `{meta['uniprot']}`")
    
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("UniProt Sequence & Structural Metadata")
        u_data = fetch_uniprot_summary(meta['uniprot'])
        if u_data['status'] == 'success':
            st.json(u_data)
        else:
            st.error(u_data['message'])
            
    with c2:
        st.subheader("TCGA Glioblastoma Mutations (cBioPortal)")
        c_data = fetch_cbioportal_mutation_stats(meta['gene'])
        if c_data['status'] == 'success':
            st.metric("Total Mutation Records", c_data['total_mutations'])
            st.write("**Top Recurrent Protein Changes:**")
            for v in c_data['sample_variants']:
                st.write(f"- {v}")
        else:
            st.error(c_data['message'])

with tab3:
    st.markdown("### Layer 2: CDC25 Phosphatase Discovery (Anchor-Bridge Model)")
    
    subtab_a, subtab_b, subtab_c = st.tabs([
        "1. Wet-Lab 4PL Fitting", 
        "2. Anchor Calibration & LOOCV", 
        "3. Mode A Manual Result Input"
    ])
    
    with subtab_a:
        st.subheader("Wet-Lab Viability Fitting (4PL)")
        st.markdown("Enter experimental crystal violet assay concentrations (µM) and cell viability (%) for NSC95397:")
        
        default_conc = "0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0"
        default_viab = "98.2, 91.5, 78.4, 32.1, 12.8, 4.2, 1.1"
        
        c_str = st.text_input("Concentrations (µM, comma-separated):", default_conc)
        v_str = st.text_input("Cell Viability (%, comma-separated):", default_viab)
        
        if st.button("Fit 4PL Dose-Response Curve"):
            try:
                concs = [float(x.strip()) for x in c_str.split(",")]
                viabs = [float(x.strip()) for x in v_str.split(",")]
                res = fit_4pl_curve(concs, viabs)
                
                if res['success']:
                    st.success(f"**Calculated IC50:** {res['ic50_uM']:.4f} µM | **Hill Slope:** {res['hill_slope']:.2f}")
                    st.pyplot(res['figure'])
                else:
                    st.error(res['error'])
            except Exception as e:
                st.error(f"Invalid input format: {e}")
                
    with subtab_b:
        st.subheader("Anchor-Set Calibration (Wet-Lab IC50 + Docking Descriptors)")
        df_anchor = load_cdc25_anchor_set()
        st.write("**Active CDC25 Anchor Set:**")
        st.dataframe(df_anchor, use_container_width=True)
        
        calib_res = calibrate_anchor_model(df_anchor)
        if calib_res['calibrated']:
            st.metric("Leave-One-Out CV RMSE (pIC50)", f"{calib_res['loocv_rmse']:.3f}")
            st.metric("Leave-One-Out CV R²", f"{calib_res['loocv_r2']:.3f}")
            
            df_cand = load_cdc25_screening_candidates()
            explain_res = compute_feature_contributions(calib_res, df_cand)
            
            st.subheader("Calibrated Prioritization Ranking")
            st.dataframe(explain_res['ranked_df'][['compound_name', 'vina_score', 'gnina_score', 'predicted_IC50_uM']], use_container_width=True)
            st.pyplot(explain_res['impact_plot'])
        else:
            st.warning(calib_res['warning'])

    with subtab_c:
        st.subheader("Mode A: Manual Result Parser")
        st.markdown("Paste JSON results from external tools (Chai-1 / Boltz-2 / Neurosnap) for single-compound evaluation:")
        
        sample_json = json.dumps({
            "compound": "NSC95397_Derivative",
            "target": "CDC25A",
            "chai1_affinity_score": 0.89,
            "boltz2_plddt": 88.4,
            "estimated_kd_nM": 142.5
        }, indent=2)
        
        user_json = st.text_area("Paste Tool Output JSON:", sample_json, height=180)
        
        if st.button("Parse & Display Result"):
            try:
                parsed = json.loads(user_json)
                st.success("Successfully Parsed External Tool Output!")
                st.json(parsed)
            except Exception as e:
                st.error(f"Invalid JSON: {e}")
