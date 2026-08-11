# Oncogbm-twin: Glioblastoma Precision Oncology & In Silico Discovery Workbench

Developed and Maintained by **Tasnim Gassem**  
*Platform Version:* 9.5 | *License:* MIT Academic Research License  
*Target Domain:* Computational Oncology, Translational Neuro-Oncology, & In Silico Drug Design  

---

## Executive Summary

**Oncogbm-twin** is an open-access computational oncology visualization portal, workflow engine, and multi-API decision-support workbench tailored for Glioblastoma Multiforme (GBM) research. The platform unifies population-scale genomic expression profiling, structural molecular docking and molecular dynamics (MD) protocols, QSAR-based organ toxicity predictors, blood-brain barrier (BBB) permeability models, live cellular invasion networks, and in vitro 4-Parameter Logistic (4PL) regression analytics into a unified web-based interface.

Designed to eliminate fragmentation across isolated computational tools, **Oncogbm-twin** provides a reproducible end-to-end pipeline that bridges transcriptomic target validation with small-molecule candidate evaluation.

 ┌─────────────────────────────────────────────────────────────────────────────────────────┐
│                        ONCOGBM-TWIN INTEGRATED WORKFLOW PIPELINE                        │
└─────────────────────────────────────────────────────────────────────────────────────────┘
                                            │
                                            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────┐
│ WORKSTATION I: Genomic & Survival Analytics                                             │
├─────────────────────────────────────────────────────────────────────────────────────────┤
│  • Cohort Analysis     : TCGA Primary GBM (N=163) vs. GTEx Normal Brain (N=207)          │
│  • Clinical Outcome    : Cox Hazard Ratios & Kaplan-Meier Survival Curves               │
│  • Network Co-Expr     : Pearson Correlation Matrix (r) Across Biomarkers               │
│  • Mutation Profiling  : Real-Time Somatic Variant Retrieval via cBioPortal REST API    │
└─────────────────────────────────────────────────────────────────────────────────────────┘
                                            │
                                            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────┐
│ WORKSTATION II: Structural Docking & 100 ns Molecular Dynamics Protocols                │
├─────────────────────────────────────────────────────────────────────────────────────────┤
│  • Pocket Docking      : Active Site Screening Thresholds (ΔG ≤ -6.0 kcal/mol)           │
│  • MD Simulations      : Solvated Explicit TIP3P 100 ns Protocols (CHARMM-GUI / GROMACS)  │
│  • Trajectory Profiling: Cα Backbone RMSD (< 2.0 Å) & Residue RMSF Fluctuation Engines  │
└─────────────────────────────────────────────────────────────────────────────────────────┘
                                            │
                                            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────┐
│ WORKSTATION III: ProTox-3 Toxicity & ADMET BBB Permeability Predictor                   │
├─────────────────────────────────────────────────────────────────────────────────────────┤
│  • Acute Oral Toxicity : OECD Guideline 423 Hazard Classification (GHS Classes 1–6)     │
│  • Organ Safety        : Deterministic Toxicity Profiling (CNS, Liver, Cardiotoxicity)  │
│  • Chemical Structure  : SMILES Graph Parsing & Property Retrieval via PubChem API      │
│  • BBB Permeability    : SwissADME BOILED-Egg Mapping (TPSA vs. WLOGP Boundaries)        │
└─────────────────────────────────────────────────────────────────────────────────────────┘
                                            │
                                            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────┐
│ WORKSTATION IV: Invasion Pathways, 4PL Kinetics & Master Literature Hub                 │
├─────────────────────────────────────────────────────────────────────────────────────────┤
│  • Cellular Invasion   : Real-Time KEGG Infiltration Pathway Search (EMT & ECM Cleavage) │
│  • In Vitro Kinetics   : Sigmoidal 4-Parameter Logistic (4PL) IC50 Curve Fitting        │
│  • Citation Export     : Master BibTeX Repository & Reference Code Exporter             │
└─────────────────────────────────────────────────────────────────────────────────────────┘
## Key Features & Methodological Implementation

### Workstation I — Genomic & Survival Analytics
* **Cohort Validation:** Integrates normalized RNA-seq transcript expression profiles from primary TCGA Glioblastoma tumors ($N=163$) and GTEx non-diseased cortical controls ($N=207$) using $\log_2(\text{TPM} + 1)$ metrics.
* **Prognostic Modeling:** Evaluates overall patient survival probability using the Cox Proportional Hazards Model and Log-rank test statistics.
* **Mutational Landscape:** Programmatically queries cBioPortal OpenAPI endpoints to retrieve recurrent missense variants (including `EGFRvIII`, `TP53` R273H, and `IDH1` R132H).

### Workstation II — Docking & Molecular Dynamics
* **In Silico Screening Guidelines:** Establishes thermodynamic hit selection thresholds ($\Delta G \le -6.0\text{ kcal/mol}$, corresponding to $K_d \le 40\ \mu\text{M}$) and hydrogen-bonding donor-acceptor distances ($\le 3.2\text{ \AA}$).
* **Trajectory Profiling:** Provides standardized protocol specifications for explicit solvation ($0.15\text{ M}$ $\text{NaCl}$ in TIP3P water box) and equilibrium trajectory analysis ($C_\alpha$ backbone RMSD $< 2.0\text{ \AA}$).

### Workstation III — Toxicity & ADMET Profiling
* **Acute Toxicity Engine:** Categorizes estimated oral $\text{LD}_{50}$ parameters ($\text{mg/kg}$) into OECD Globally Harmonized System (GHS) hazard classes 1 through 6.
* **BBB Penetration Modeling:** Dynamically parses topological polar surface area ($\text{TPSA}$) and lipophilicity ($\text{WLOGP}$) from canonical SMILES strings to mathematically project candidate molecules onto the SwissADME BOILED-Egg permeation space ($\text{TPSA} < 75\text{ \AA}^2, 0.5 < \text{WLOGP} < 3.5$).
* **Database Interoperability:** Includes one-click external routing to NCBI PubChem entries.

### Workstation IV — Cellular Kinetics & Pathway Search
* **KEGG Invasion Querying:** Queries the KEGG REST API to identify signaling cascades involved in epithelial-mesenchymal transition (EMT), focal adhesion, and matrix metalloproteinase ($\text{MMP2}/\text{MMP9}$) activation.
* **4PL Non-Linear Regression:** Implements `scipy.optimize.curve_fit` to calculate half-maximal inhibitory concentrations ($\text{IC}_{50}$), Hill slopes ($b$), and coefficients of determination ($R^2$) from user-provided in vitro dose-response datasets:
  $$y = d + \frac{a - d}{1 + \left(\frac{x}{c}\right)^b}$$

---

## Installation & Setup Guide

### Prerequisites
* Python 3.9 or higher
* `pip` package installer

### Environment Configuration

1. **Clone the Repository:**
   ```bash
   git clone [https://github.com/TiaGassem/Oncogbm-twin.git](https://github.com/TiaGassem/Oncogbm-twin.git)
   cd Oncogbm-twin
Create and Activate a Virtual Environment:Linux/macOS:Bashpython3 -m venv venv
source venv/bin/activate
Windows:DOSpython -m venv venv
venv\Scripts\activate
Install Dependencies:Create a requirements.txt file (or use the one provided) containing:Plaintextstreamlit>=1.30.0
requests>=2.31.0
numpy>=1.24.0
pandas>=2.0.0
matplotlib>=3.7.0
scipy>=1.10.0
Install via pip:Bashpip install -r requirements.txt
Launch the Web Workbench:Bashstreamlit run app.py
External Data Sources & Web API IntegrationsThe platform programmatically connects with the following public biomedical repositories:Resource NameData DomainAccess ProtocolNCBI PubChemSMILES, Molecular Weight, TPSA, WLOGP, H-BondsPUG-REST APIcBioPortalTCGA Glioblastoma Somatic Mutation FrequenciesOpen REST APIKEGG DBGlioblastoma Invasion & Cell Migration PathwaysKEGG REST APIRCSB PDBMacromolecular Receptor 3D StructuresDirect HTTPSProTox 3.0Acute Toxicity & Organ Safety ReferencesPublic Web ServerSwissADMEBOILED-Egg Gastrointestinal/BBB BoundariesMethodological StandardCitation & Academic ReferencesIf you utilize the Oncogbm-twin platform or its analytical workflows in your research or thesis, please cite the underlying methodologies as follows:Extrait de code@article{banerjee2024protox,
  title={ProTox 3.0: a webserver for the prediction of toxicities of small molecules},
  author={Banerjee, Preeti and Kemmler, Eva and Dunkel, Mathias and Preissner, Robert},
  journal={Nucleic Acids Research},
  volume={52},
  number={W1},
  pages={W513--W520},
  year={2024},
  doi={10.1093/nar/gkae303}
}

@article{daina2016boiled,
  title={A BOILED-Egg To Predict Gastrointestinal Absorption and Brain Penetration of Small Molecules},
  author={Daina, Antoine and Zoete, Vincent},
  journal={ChemMedChem},
  volume={11},
  number={11},
  pages={1117--1121},
  year={2016},
  doi={10.1002/cmdc.201600182}
}

@article{meng2011molecular,
  title={Molecular Docking: a powerful approach for structure-based drug discovery},
  author={Meng, Xiao-Yin and Zhang, Hong-Xing and Mezei, Mihaly and Cui, Meng},
  journal={Current Computer-Aided Drug Design},
  volume={7},
  number={2},
  pages={146--157},
  year={2011},
  pmid={21532826}
}

@article{tcga2008comprehensive,
  title={Comprehensive genomic characterization defines human glioblastoma genes and core pathways},
  author={{TCGA Research Network}},
  journal={Nature},
  volume={455},
  number={7216},
  pages={1061--1068},
  year={2008},
  pmid={18772890}
}


License & Intellectual PropertyDistributed under the MIT Academic Research License.Copyright © 2026 Tasnim Gassem. All rights reserved.
