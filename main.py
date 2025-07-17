import streamlit as st
import pandas as pd
import json
import plotly.express as px
import plotly.graph_objects as go
import requests
from io import BytesIO

st.set_page_config(page_title="Chemmy", layout="wide")

# Load data unsur
@st.cache_data
def load_elements():
    # Sample periodic table data - replace with your JSON file
    sample_elements = [
        {"symbol": "H", "name": "Hydrogen", "number": 1, "group": 1, "electronegativity": 2.20, "electron_configuration": "1s1"},
        {"symbol": "He", "name": "Helium", "number": 2, "group": 18, "electronegativity": None, "electron_configuration": "1s2"},
        {"symbol": "Li", "name": "Lithium", "number": 3, "group": 1, "electronegativity": 0.98, "electron_configuration": "[He] 2s1"},
        {"symbol": "Be", "name": "Beryllium", "number": 4, "group": 2, "electronegativity": 1.57, "electron_configuration": "[He] 2s2"},
        {"symbol": "B", "name": "Boron", "number": 5, "group": 13, "electronegativity": 2.04, "electron_configuration": "[He] 2s2 2p1"},
        {"symbol": "C", "name": "Carbon", "number": 6, "group": 14, "electronegativity": 2.55, "electron_configuration": "[He] 2s2 2p2"},
        {"symbol": "N", "name": "Nitrogen", "number": 7, "group": 15, "electronegativity": 3.04, "electron_configuration": "[He] 2s2 2p3"},
        {"symbol": "O", "name": "Oxygen", "number": 8, "group": 16, "electronegativity": 3.44, "electron_configuration": "[He] 2s2 2p4"},
        {"symbol": "F", "name": "Fluorine", "number": 9, "group": 17, "electronegativity": 3.98, "electron_configuration": "[He] 2s2 2p5"},
        {"symbol": "Ne", "name": "Neon", "number": 10, "group": 18, "electronegativity": None, "electron_configuration": "[He] 2s2 2p6"},
    ]
    
    try:
        with open("PeriodicTableJSON.json") as file:
            elements = json.load(file)["elements"]
    except FileNotFoundError:
        st.warning("PeriodicTableJSON.json not found. Using sample data.")
        elements = sample_elements
    
    return pd.DataFrame(elements)

# Molecular data dictionary
MOLECULES = {
    "CCO": {"name": "Ethanol", "formula": "C2H6O", "mw": 46.07},
    "C1=CC=CC=C1": {"name": "Benzene", "formula": "C6H6", "mw": 78.11},
    "CC(=O)O": {"name": "Acetic Acid", "formula": "C2H4O2", "mw": 60.05},
    "O": {"name": "Water", "formula": "H2O", "mw": 18.02},
    "C": {"name": "Methane", "formula": "CH4", "mw": 16.04},
    "CC(=O)OC1=CC=CC=C1C(=O)O": {"name": "Aspirin", "formula": "C9H8O4", "mw": 180.16},
    "C(C1C(C(C(C(O1)O)O)O)O)O": {"name": "Glucose", "formula": "C6H12O6", "mw": 180.16},
}

def get_molecule_image_url(smiles):
    """Generate molecule image URL using external service"""
    # Using a free chemical structure service
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles"
    try:
        # This is a simplified approach - in practice you'd want more robust handling
        encoded_smiles = smiles.replace("=", "%3D").replace("(", "%28").replace(")", "%29")
        url = f"{base_url}/{encoded_smiles}/PNG"
        return url
    except:
        return None

df = load_elements()

st.title("🧪 ChemExplorer 2.0")

tab1, tab2 = st.tabs(["🔬 Tabel Periodik", "🧬 Visualisasi Molekul"])

# ========= TAB 1 ============
with tab1:
    st.header("📘 Info Unsur")
    
    selected_symbol = st.sidebar.selectbox("Pilih simbol unsur", df["symbol"].sort_values())
    elemen = df[df["symbol"] == selected_symbol].iloc[0]
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader(f"{elemen['name']} ({elemen['symbol']})")
        st.markdown(f"**Nomor Atom:** {elemen['number']}")
        st.markdown(f"**Golongan:** {elemen.get('group', '-')}")
        st.markdown(f"**Elektronegativitas:** {elemen.get('electronegativity', '-')}")
        st.markdown(f"**Konfigurasi Elektron:** {elemen.get('electron_configuration', '-')}")
    
    with col2:
        st.subheader("📈 Tren Elektronegativitas")
        
        # Filter out None values for electronegativity
        df_filtered = df.dropna(subset=['electronegativity'])
        
        fig = px.line(df_filtered.sort_values("number"), 
                      x="number", 
                      y="electronegativity",
                      labels={"number": "Nomor Atom", "electronegativity": "Elektronegativitas"},
                      title="Elektronegativitas vs Nomor Atom")
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

# ========= TAB 2 ==============
with tab2:
    st.header("🧬 Visualisasi Struktur Molekul")
    
    col1, col2 = st.columns([2, 3])
    
    with col1:
        input_type = st.radio("Input molekul via:", ["SMILES", "Pilih dari contoh"])
        
        if input_type == "SMILES":
            smiles = st.text_input("Masukkan SMILES:", value="CCO")
        else:
            molecule_options = list(MOLECULES.keys())
            molecule_names = [f"{MOLECULES[smiles]['name']} ({smiles})" for smiles in molecule_options]
            selected_idx = st.selectbox("Pilih molekul:", range(len(molecule_names)), 
                                      format_func=lambda x: molecule_names[x])
            smiles = molecule_options[selected_idx]
        
        if st.button("Tampilkan Molekul"):
            if smiles:
                st.success("Memproses struktur molekul...")
                
                # Try to get molecule info
                if smiles in MOLECULES:
                    mol_info = MOLECULES[smiles]
                    st.subheader(f"📋 {mol_info['name']}")
                    st.markdown(f"**SMILES:** {smiles}")
                    st.markdown(f"**Formula Molekul:** {mol_info['formula']}")
                    st.markdown(f"**Berat Molekul:** {mol_info['mw']} g/mol")
                else:
                    st.subheader("📋 Molekul")
                    st.markdown(f"**SMILES:** {smiles}")
                
                # Try to display image from external service
                st.subheader("🖼️ Struktur 2D")
                image_url = get_molecule_image_url(smiles)
                if image_url:
                    try:
                        st.image(image_url, caption=f"Struktur molekul: {smiles}")
                    except:
                        st.info("Tidak dapat menampilkan gambar struktur. Silakan coba molekul lain.")
                else:
                    st.info("Gambar struktur tidak tersedia untuk SMILES ini.")
                
                # Simple molecular analysis
                st.subheader("📊 Analisis Sederhana")
                carbon_count = smiles.count('C')
                oxygen_count = smiles.count('O')
                nitrogen_count = smiles.count('N')
                
                if carbon_count > 0:
                    st.markdown(f"**Perkiraan atom Karbon:** {carbon_count}")
                if oxygen_count > 0:
                    st.markdown(f"**Perkiraan atom Oksigen:** {oxygen_count}")
                if nitrogen_count > 0:
                    st.markdown(f"**Perkiraan atom Nitrogen:** {nitrogen_count}")
                
                # Ring detection
                if '1' in smiles or '2' in smiles or '3' in smiles:
                    st.markdown("**Struktur:** Mengandung cincin/ring")
                else:
                    st.markdown("**Struktur:** Struktur rantai")
            else:
                st.error("Silakan masukkan SMILES yang valid.")
    
    with col2:
        st.markdown("""
        #### Contoh SMILES:
        - **Benzena**: C1=CC=CC=C1
        - **Ethanol**: CCO
        - **Glukosa**: C(C1C(C(C(C(O1)O)O)O)O)O
        - **Asam asetat**: CC(=O)O
        - **Aspirin**: CC(=O)OC1=CC=CC=C1C(=O)O
        - **Air**: O
        - **Metana**: C
        """)
        
        st.markdown("""
        #### Tips:
        - Gunakan SMILES yang valid untuk hasil terbaik
        - Aplikasi akan mencoba menampilkan struktur dari layanan eksternal
        - Pilih dari contoh molekul untuk hasil yang lebih baik
        """)
        
        st.markdown("""
        #### Catatan:
        - Visualisasi struktur menggunakan layanan eksternal
        - Beberapa struktur mungkin tidak tersedia
        - Analisis molekul bersifat sederhana
        """)

# Footer
st.markdown("---")
st.markdown("🧪 **ChemExplorer 2.0** - Eksplorasi Kimia Interaktif (Tanpa RDKit)")
