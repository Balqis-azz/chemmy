import streamlit as st
import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import plotly.express as px  # Missing import

st.set_page_config(page_title="Chemmy", layout="wide")

# Load data unsur
with open("PeriodicTableJSON.json") as file:
    elements = json.load(file)["elements"]

df = pd.DataFrame(elements)

st.title("🧪 ChemExplorer 2.0")

tab1, tab2 = st.tabs(["🔬 Tabel Periodik", "🧬 Visualisasi Molekul"])

# ========= TAB 1 ============
with tab1:
    st.header("📘 Info Unsur")
    
    selected_symbol = st.sidebar.selectbox("Pilih simbol unsur", df["symbol"].sort_values())
    elemen = df[df["symbol"] == selected_symbol].iloc[0]
    
    st.subheader(f"{elemen['name']} ({elemen['symbol']})")
    st.markdown(f"**Nomor Atom:** {elemen['number']}")
    st.markdown(f"**Golongan:** {elemen.get('group', '-')}")
    st.markdown(f"**Elektronegativitas:** {elemen.get('electronegativity', '-')}")
    st.markdown(f"**Konfigurasi Elektron:** {elemen.get('electron_configuration', '-')}")
    
    st.subheader("📈 Tren Elektronegativitas")
    
    # Filter out None values for electronegativity
    df_filtered = df.dropna(subset=['electronegativity'])
    
    fig = px.line(df_filtered.sort_values("number"), 
                  x="number", 
                  y="electronegativity",
                  labels={"number": "Nomor Atom", "electronegativity": "Elektronegativitas"},
                  title="Elektronegativitas vs Nomor Atom")
    st.plotly_chart(fig, use_container_width=True)

# ========= TAB 2 ==============
with tab2:
    st.header("🧬 Visualisasi Struktur Molekul")
    
    col1, col2 = st.columns([2, 3])
    
    with col1:
        input_type = st.radio("Input molekul via:", ["SMILES", "Nama (pakai SMILES manual)"])
        
        if input_type == "SMILES":
            smiles = st.text_input("Masukkan SMILES (contoh: C1=CC=CC=C1 untuk benzena):", value="CCO")
        else:
            smiles = st.text_input("Masukkan nama atau SMILES-nya (sementara gunakan SMILES):", value="CCO")
        
        if st.button("Tampilkan Molekul"):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    st.success("Struktur dikenali!")
                    # Generate molecule image
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img)
                    
                    # Display molecular properties
                    st.subheader("Properti Molekul")
                    st.markdown(f"**Formula Molekul:** {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
                    st.markdown(f"**Berat Molekul:** {Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}")
                    st.markdown(f"**Jumlah Atom:** {mol.GetNumAtoms()}")
                    st.markdown(f"**Jumlah Ikatan:** {mol.GetNumBonds()}")
                else:
                    st.error("Struktur tidak valid atau tidak dikenali.")
            except Exception as e:
                st.error(f"Error: {str(e)}")
    
    with col2:
        st.markdown("""
        #### Contoh SMILES:
        - **Benzena**: C1=CC=CC=C1
        - **Ethanol**: CCO
        - **Glukosa**: C(C1C(C(C(C(O1)O)O)O)O)O
        - **Asam asetat**: CC(=O)O
        - **Aspirin**: CC(=O)OC1=CC=CC=C1C(=O)O
        - **Kafein**: CN1C=NC2=C1C(=O)N(C(=O)N2C)C
        - **Air**: O
        - **Metana**: C
        """)
        
        st.markdown("""
        #### Tips:
        - Gunakan SMILES yang valid untuk hasil terbaik
        - Aplikasi akan menampilkan struktur 2D molekul
        - Informasi tambahan seperti berat molekul akan ditampilkan
        """)

# Footer
st.markdown("---")
st.markdown("🧪 **ChemExplorer 2.0** - Eksplorasi Kimia Interaktif")
