import streamlit as st
import pandas as pd
import json

st.set_page_config(page_title="Chemmy", layout="wide")

@st.cache_data
def load_elements():
    try:
        with open("PeriodicTableJSON.json") as file:
            data = json.load(file)
            
        # Handle different possible structures
        if "elements" in data:
            elements = data["elements"]
        else:
            elements = data
            
        # Debug: Show first element structure
        if elements:
            st.sidebar.write("📋 Struktur data yang ditemukan:")
            st.sidebar.write(f"Kolom tersedia: {list(elements[0].keys())}")
            
        return pd.DataFrame(elements)
        
    except FileNotFoundError:
        st.warning("PeriodicTableJSON.json tidak ditemukan. Menggunakan data sampel.")
        # Sample periodic table data
        sample_elements = [
            {"symbol": "H", "name": "Hydrogen", "number": 1, "group": 1, "electronegativity": 2.20, "electron_configuration": "1s1", "atomic_mass": 1.008},
            {"symbol": "He", "name": "Helium", "number": 2, "group": 18, "electronegativity": None, "electron_configuration": "1s2", "atomic_mass": 4.003},
            {"symbol": "Li", "name": "Lithium", "number": 3, "group": 1, "electronegativity": 0.98, "electron_configuration": "[He] 2s1", "atomic_mass": 6.94},
            {"symbol": "Be", "name": "Beryllium", "number": 4, "group": 2, "electronegativity": 1.57, "electron_configuration": "[He] 2s2", "atomic_mass": 9.01},
            {"symbol": "B", "name": "Boron", "number": 5, "group": 13, "electronegativity": 2.04, "electron_configuration": "[He] 2s2 2p1", "atomic_mass": 10.81},
            {"symbol": "C", "name": "Carbon", "number": 6, "group": 14, "electronegativity": 2.55, "electron_configuration": "[He] 2s2 2p2", "atomic_mass": 12.01},
            {"symbol": "N", "name": "Nitrogen", "number": 7, "group": 15, "electronegativity": 3.04, "electron_configuration": "[He] 2s2 2p3", "atomic_mass": 14.01},
            {"symbol": "O", "name": "Oxygen", "number": 8, "group": 16, "electronegativity": 3.44, "electron_configuration": "[He] 2s2 2p4", "atomic_mass": 16.00},
            {"symbol": "F", "name": "Fluorine", "number": 9, "group": 17, "electronegativity": 3.98, "electron_configuration": "[He] 2s2 2p5", "atomic_mass": 19.00},
            {"symbol": "Ne", "name": "Neon", "number": 10, "group": 18, "electronegativity": None, "electron_configuration": "[He] 2s2 2p6", "atomic_mass": 20.18},
            {"symbol": "Na", "name": "Sodium", "number": 11, "group": 1, "electronegativity": 0.93, "electron_configuration": "[Ne] 3s1", "atomic_mass": 22.99},
            {"symbol": "Mg", "name": "Magnesium", "number": 12, "group": 2, "electronegativity": 1.31, "electron_configuration": "[Ne] 3s2", "atomic_mass": 24.31},
            {"symbol": "Al", "name": "Aluminum", "number": 13, "group": 13, "electronegativity": 1.61, "electron_configuration": "[Ne] 3s2 3p1", "atomic_mass": 26.98},
            {"symbol": "Si", "name": "Silicon", "number": 14, "group": 14, "electronegativity": 1.90, "electron_configuration": "[Ne] 3s2 3p2", "atomic_mass": 28.09},
            {"symbol": "P", "name": "Phosphorus", "number": 15, "group": 15, "electronegativity": 2.19, "electron_configuration": "[Ne] 3s2 3p3", "atomic_mass": 30.97},
            {"symbol": "S", "name": "Sulfur", "number": 16, "group": 16, "electronegativity": 2.58, "electron_configuration": "[Ne] 3s2 3p4", "atomic_mass": 32.06},
            {"symbol": "Cl", "name": "Chlorine", "number": 17, "group": 17, "electronegativity": 3.16, "electron_configuration": "[Ne] 3s2 3p5", "atomic_mass": 35.45},
            {"symbol": "Ar", "name": "Argon", "number": 18, "group": 18, "electronegativity": None, "electron_configuration": "[Ne] 3s2 3p6", "atomic_mass": 39.95},
            {"symbol": "K", "name": "Potassium", "number": 19, "group": 1, "electronegativity": 0.82, "electron_configuration": "[Ar] 4s1", "atomic_mass": 39.10},
            {"symbol": "Ca", "name": "Calcium", "number": 20, "group": 2, "electronegativity": 1.00, "electron_configuration": "[Ar] 4s2", "atomic_mass": 40.08},
        ]
        return pd.DataFrame(sample_elements)
    
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return pd.DataFrame()  # Return empty DataFrame on error

# Molecular data dictionary
MOLECULES = {
    "CCO": {"name": "Ethanol", "formula": "C₂H₆O", "mw": 46.07, "description": "Alkohol yang umum digunakan"},
    "C1=CC=CC=C1": {"name": "Benzene", "formula": "C₆H₆", "mw": 78.11, "description": "Senyawa aromatik dasar"},
    "CC(=O)O": {"name": "Acetic Acid", "formula": "C₂H₄O₂", "mw": 60.05, "description": "Asam asetat (cuka)"},
    "O": {"name": "Water", "formula": "H₂O", "mw": 18.02, "description": "Air"},
    "C": {"name": "Methane", "formula": "CH₄", "mw": 16.04, "description": "Gas alam utama"},
    "CC(=O)OC1=CC=CC=C1C(=O)O": {"name": "Aspirin", "formula": "C₉H₈O₄", "mw": 180.16, "description": "Obat penghilang rasa sakit"},
    "C(C1C(C(C(C(O1)O)O)O)O)O": {"name": "Glucose", "formula": "C₆H₁₂O₆", "mw": 180.16, "description": "Gula sederhana"},
    "CC(C)O": {"name": "Isopropanol", "formula": "C₃H₈O", "mw": 60.10, "description": "Alkohol isopropil"},
    "CCO": {"name": "Ethanol", "formula": "C₂H₆O", "mw": 46.07, "description": "Alkohol etil"},
}

def analyze_smiles(smiles):
    """Simple SMILES analysis"""
    analysis = {}
    
    # Count elements
    analysis['carbon'] = smiles.count('C')
    analysis['oxygen'] = smiles.count('O')
    analysis['nitrogen'] = smiles.count('N')
    analysis['sulfur'] = smiles.count('S')
    analysis['phosphorus'] = smiles.count('P')
    
    # Detect rings
    analysis['has_ring'] = any(char.isdigit() for char in smiles)
    
    # Detect double bonds
    analysis['has_double_bond'] = '=' in smiles
    
    # Detect aromatic
    analysis['aromatic'] = 'c' in smiles or (analysis['has_ring'] and analysis['has_double_bond'])
    
    return analysis

def create_electronegativity_chart(df):
    """Create electronegativity chart using matplotlib"""
    # Check if electronegativity column exists
    if 'electronegativity' not in df.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'Data elektronegativitas tidak tersedia', 
                horizontalalignment='center', verticalalignment='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title('Elektronegativitas vs Nomor Atom')
        return fig
    
    df_filtered = df.dropna(subset=['electronegativity'])
    
    if df_filtered.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'Tidak ada data elektronegativitas yang valid', 
                horizontalalignment='center', verticalalignment='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title('Elektronegativitas vs Nomor Atom')
        return fig
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(df_filtered['number'], df_filtered['electronegativity'], 'bo-', linewidth=2, markersize=6)
    ax.set_xlabel('Nomor Atom')
    ax.set_ylabel('Elektronegativitas')
    ax.set_title('Elektronegativitas vs Nomor Atom')
    ax.grid(True, alpha=0.3)
    
    return fig

df = load_elements()

st.title("🧪 ChemExplorer 2.0")
st.sidebar.markdown("### 🔧 Navigasi")

tab1, tab2 = st.tabs(["🔬 Tabel Periodik", "🧬 Analisis Molekul"])

# ========= TAB 1 ============
with tab1:
    st.header("📘 Informasi Unsur")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        selected_symbol = st.selectbox("Pilih simbol unsur:", df["symbol"].sort_values())
        elemen = df[df["symbol"] == selected_symbol].iloc[0]
        
        st.subheader(f"🔬 {elemen['name']} ({elemen['symbol']})")
        
        info_container = st.container()
        with info_container:
            st.markdown(f"**Nomor Atom:** {elemen['number']}")
            st.markdown(f"**Golongan:** {elemen.get('group', 'Tidak diketahui')}")
            
            # Handle different possible column names for electronegativity
            electronegativity = elemen.get('electronegativity') or elemen.get('electronegativity_pauling') or elemen.get('electronegativity_allred_rochow')
            st.markdown(f"**Elektronegativitas:** {electronegativity if electronegativity else 'Tidak diketahui'}")
            
            # Handle different possible column names for atomic mass
            atomic_mass = elemen.get('atomic_mass') or elemen.get('atomic_weight') or elemen.get('standard_atomic_weight')
            st.markdown(f"**Massa Atom:** {atomic_mass if atomic_mass else 'Tidak diketahui'}")
            
            # Handle different possible column names for electron configuration
            electron_config = elemen.get('electron_configuration') or elemen.get('electron_configuration_semantic') or elemen.get('electronic_configuration')
            st.markdown(f"**Konfigurasi Elektron:** `{electron_config if electron_config else 'Tidak diketahui'}`")
        
        # Additional info
        if electronegativity and electronegativity != 'Tidak diketahui':
            try:
                en_value = float(electronegativity)
                if en_value > 3.0:
                    st.info("🔸 Unsur sangat elektronegatif")
                elif en_value > 2.0:
                    st.info("🔹 Unsur elektronegatif")
                else:
                    st.info("🔹 Unsur elektropositif")
            except (ValueError, TypeError):
                pass
    
    with col2:
        st.subheader("📈 Tren Elektronegativitas")
        
        # Create and display the chart
        fig = create_electronegativity_chart(df)
        st.pyplot(fig)
        
        # Show periodic trends
        st.subheader("📊 Tren Periodik")
        
        # Check if electronegativity column exists and has data
        if 'electronegativity' in df.columns:
            df_filtered = df.dropna(subset=['electronegativity'])
            
            if not df_filtered.empty:
                max_en = df_filtered.loc[df_filtered['electronegativity'].idxmax()]
                min_en = df_filtered.loc[df_filtered['electronegativity'].idxmin()]
                
                st.markdown(f"**Elektronegativitas Tertinggi:** {max_en['name']} ({max_en['electronegativity']})")
                st.markdown(f"**Elektronegativitas Terendah:** {min_en['name']} ({min_en['electronegativity']})")
            else:
                st.info("Data elektronegativitas tidak tersedia untuk analisis tren")
        else:
            st.info("Data elektronegativitas tidak tersedia dalam dataset")

# ========= TAB 2 ==============
with tab2:
    st.header("🧬 Analisis Struktur Molekul")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📝 Input Molekul")
        
        input_method = st.radio("Metode input:", ["Pilih dari contoh", "Input SMILES manual"])
        
        if input_method == "Pilih dari contoh":
            molecule_options = list(MOLECULES.keys())
            molecule_labels = [f"{MOLECULES[smiles]['name']} ({MOLECULES[smiles]['formula']})" for smiles in molecule_options]
            
            selected_idx = st.selectbox("Pilih molekul:", range(len(molecule_labels)), 
                                      format_func=lambda x: molecule_labels[x])
            smiles = molecule_options[selected_idx]
            
            # Show molecule info
            mol_info = MOLECULES[smiles]
            st.info(f"**Deskripsi:** {mol_info['description']}")
            
        else:
            smiles = st.text_input("Masukkan SMILES:", value="CCO", help="Contoh: CCO untuk ethanol")
        
        if st.button("🔍 Analisis Molekul", type="primary"):
            if smiles:
                analysis = analyze_smiles(smiles)
                
                st.success("✅ Analisis berhasil!")
                
                # Basic info
                st.subheader("📋 Informasi Dasar")
                if smiles in MOLECULES:
                    mol_info = MOLECULES[smiles]
                    st.markdown(f"**Nama:** {mol_info['name']}")
                    st.markdown(f"**Formula:** {mol_info['formula']}")
                    st.markdown(f"**Berat Molekul:** {mol_info['mw']} g/mol")
                
                st.markdown(f"**SMILES:** `{smiles}`")
                
                # Composition analysis
                st.subheader("🔬 Komposisi Unsur")
                composition = []
                if analysis['carbon'] > 0:
                    composition.append(f"Karbon: {analysis['carbon']} atom")
                if analysis['oxygen'] > 0:
                    composition.append(f"Oksigen: {analysis['oxygen']} atom")
                if analysis['nitrogen'] > 0:
                    composition.append(f"Nitrogen: {analysis['nitrogen']} atom")
                if analysis['sulfur'] > 0:
                    composition.append(f"Sulfur: {analysis['sulfur']} atom")
                
                if composition:
                    for comp in composition:
                        st.markdown(f"• {comp}")
                else:
                    st.markdown("• Tidak ada unsur yang terdeteksi dalam notasi SMILES")
                
                # Structural features
                st.subheader("🏗️ Fitur Struktural")
                if analysis['has_ring']:
                    st.markdown("• ✅ Mengandung struktur cincin")
                else:
                    st.markdown("• ➖ Struktur rantai terbuka")
                
                if analysis['has_double_bond']:
                    st.markdown("• ✅ Mengandung ikatan rangkap")
                else:
                    st.markdown("• ➖ Hanya ikatan tunggal")
                
                if analysis['aromatic']:
                    st.markdown("• ✅ Kemungkinan aromatik")
                else:
                    st.markdown("• ➖ Non-aromatik")
            else:
                st.error("❌ Silakan masukkan SMILES yang valid")
    
    with col2:
        st.subheader("📚 Panduan SMILES")
        
        with st.expander("🔤 Notasi Dasar"):
            st.markdown("""
            **Atom:**
            - C = Karbon
            - O = Oksigen  
            - N = Nitrogen
            - S = Sulfur
            - P = Fosfor
            
            **Ikatan:**
            - - = Ikatan tunggal (opsional)
            - = = Ikatan rangkap dua
            - # = Ikatan rangkap tiga
            """)
        
        with st.expander("🔄 Cincin"):
            st.markdown("""
            **Cincin menggunakan angka:**
            - C1CCCCC1 = Sikloheksana
            - C1=CC=CC=C1 = Benzena
            - Angka menunjukkan koneksi cincin
            """)
        
        with st.expander("🧪 Contoh Molekul"):
            st.markdown("""
            **Molekul Sederhana:**
            - CCO = Ethanol
            - CC(=O)O = Asam asetat
            - O = Air
            - C = Metana
            
            **Molekul Kompleks:**
            - C1=CC=CC=C1 = Benzena
            - CC(=O)OC1=CC=CC=C1C(=O)O = Aspirin
            """)
        
        st.subheader("ℹ️ Informasi")
        st.info("""
        Aplikasi ini menggunakan analisis SMILES sederhana untuk memberikan informasi dasar tentang molekul. 
        Untuk visualisasi struktur yang lebih akurat, gunakan software kimia khusus seperti ChemDraw atau Avogadro.
        """)

# Footer
st.markdown("---")
st.markdown("🧪 **ChemExplorer 2.0** - Eksplorasi Kimia Interaktif")
st.markdown("*Dibuat dengan ❤️ menggunakan Streamlit*")
