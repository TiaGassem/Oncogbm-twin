import requests
import streamlit as st

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
        seq_len = data.get("sequence", {}).get("length", 0)
        
        return {
            "status": "success",
            "accession": uniprot_id,
            "full_name": rec_name,
            "length": seq_len,
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
