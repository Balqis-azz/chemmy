import streamlit as st
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import stmol

# ===================== FUNGSI BARU ===================== #
def show_3d_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.warning("SMILES tidak valid!")
            return
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        mb = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(width=400, height=400)
        view.addModel(mb, 'mol')
        view.setStyle({'stick': {}})
        view.zoomTo()
        view.spin()
        stmol.showmol(view, height=400, width=400)
        
    except Exception as e:
        st.error(f"Error: {e}")

# ===================== KODE YANG SUDAH ADA ===================== #
# (Fungsi load_elements(), create_electronegativity_chart(), analyze_smiles(), 
#  dan bagian UI tetap sama seperti sebelumnya)

with tab2:
    st.header("🧬 Analisis Struktur Molekul")
    col1, col2 = st.columns([1, 1])

    with col1:
        # [Input SMILES dan analisis yang sudah ada...]

    with col2:
        # ✅ Tambahkan visualisasi di sini
        if 'smiles' in locals() and smiles:
            st.subheader("🌐 Model 3D")
            show_3d_molecule(smiles)
        
        # [Panduan SMILES yang sudah ada...]
