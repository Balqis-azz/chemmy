import streamlit as st
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np

# Tambahkan import untuk 3D molekul
try:
    import stmol
    import py3Dmol
    STMOL_AVAILABLE = True
except ImportError:
    STMOL_AVAILABLE = False
    st.warning("⚠️ stmol tidak terinstal. Jalankan: pip install stmol")

st.set_page_config(page_title="Chemmy", layout="wide")

# ===================== LOAD DATA ===================== #
@st.cache_data
def load_elements():
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
    ]
    try:
        with open("PeriodicTableJSON.json") as file:
            elements = json.load(file)["elements"]
    except Exception:
        st.warning("⚠️ File JSON tidak ditemukan atau tidak valid. Menggunakan data contoh.")
        elements = sample_elements
    return pd.DataFrame(elements)

# ===================== 3D MOLECULE VIEWER ===================== #
def show_molecule_3d(smiles, molecule_name="Molekul"):
    """Tampilkan molekul dalam 3D menggunakan py3Dmol"""
    if not STMOL_AVAILABLE:
        st.error("stmol tidak tersedia. Install dengan: pip install stmol")
        return
    
    try:
        # Buat viewer 3D
        viewer = py3Dmol.view(width=500, height=400)
        
        # Tambahkan molekul dari SMILES
        viewer.addModel(smiles, 'smi')
        
        # Set style visualisasi
        viewer.setStyle({'stick': {'radius': 0.2}, 'sphere': {'radius': 0.4}})
        viewer.setBackgroundColor('white')
        viewer.zoomTo()
        viewer.spin(True)  # Rotasi otomatis
        
        # Tampilkan di Streamlit
        st.subheader(f"🔬 Struktur 3D: {molecule_name}")
        stmol.showmol(viewer, height=400, width=500)
        
        # Kontrol tambahan
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("🔄 Reset View"):
                viewer.zoomTo()
        with col2:
            if st.button("⏸️ Stop Rotation"):
                viewer.spin(False)
        with col3:
            if st.button("▶️ Start Rotation"):
                viewer.spin(True)
                
    except Exception as e:
        st.error(f"Gagal menampilkan molekul 3D: {e}")

# ===================== CHART ===================== #
def create_electronegativity_chart(df):
    try:
        df['number'] = pd.to_numeric(df['number'], errors='coerce')
        df['electronegativity'] = pd.to_numeric(df['electronegativity'], errors='coerce')
        df_filtered = df.dropna(subset=['number', 'electronegativity'])
        if df_filtered.empty:
            st.warning("⚠️ Tidak ada data elektronegativitas yang valid.")
            return None
        df_filtered = df_filtered.sort_values("number")
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(df_filtered['number'], df_filtered['electronegativity'], 'bo-', linewidth=2, markersize=6)
        ax.set_xlabel('Nomor Atom')
        ax.set_ylabel('Elektronegativitas')
        ax.set_title('Elektronegativitas vs Nomor Atom')
        ax.grid(True, alpha=0.3)
        return fig
    except Exception as e:
        st.error(f"Gagal membuat grafik: {e}")
        return None

# ===================== ANALISIS SMILES ===================== #
def analyze_smiles(smiles):
    analysis = {}
    analysis['carbon'] = smiles.count('C')
    analysis['oxygen'] = smiles.count('O')
    analysis['nitrogen'] = smiles.count('N')
    analysis['sulfur'] = smiles.count('S')
    analysis['phosphorus'] = smiles.count('P')
    analysis['has_ring'] = any(char.isdigit() for char in smiles)
    analysis['has_double_bond'] = '=' in smiles
    analysis['aromatic'] = 'c' in smiles or (analysis['has_ring'] and analysis['has_double_bond'])
    return analysis

# ===================== MOLEKUL ===================== #
MOLECULES = {
    "CCO": {"name": "Ethanol", "formula": "C₂H₆O", "mw": 46.07, "description": "Alkohol yang umum digunakan"},
    "C1=CC=CC=C1": {"name": "Benzene", "formula": "C₆H₆", "mw": 78.11, "description": "Senyawa aromatik dasar"},
    "CC(=O)O": {"name": "Acetic Acid", "formula": "C₂H₄O₂", "mw": 60.05, "description": "Asam asetat (cuka)"},
    "O": {"name": "Water", "formula": "H₂O", "mw": 18.02, "description": "Air"},
    "C": {"name": "Methane", "formula": "CH₄", "mw": 16.04, "description": "Gas alam utama"},
    "CC(=O)OC1=CC=CC=C1C(=O)O": {"name": "Aspirin", "formula": "C₉H₈O₄", "mw": 180.16, "description": "Obat penghilang rasa sakit"},
    "C(C1C(C(C(C(O1)O)O)O)O)O": {"name": "Glucose", "formula": "C₆H₁₂O₆", "mw": 180.16, "description": "Gula sederhana"},
    "CC(C)O": {"name": "Isopropanol", "formula": "C₃H₈O", "mw": 60.10, "description": "Alkohol isopropil"},
}

# ===================== LOAD ===================== #
df = load_elements()

# ===================== UI START ===================== #
st.title("🧪 ChemExplorer 3D")
st.markdown("*Eksplorasi Kimia dengan Visualisasi 3D Interaktif*")

tab1, tab2 = st.tabs(["🔬 Tabel Periodik", "🧬 Analisis Molekul 3D"])

# ---------------- TAB 1 ----------------
with tab1:
    st.header("📘 Informasi Unsur")
    col1, col2 = st.columns([1, 2])

    with col1:
        selected_symbol = st.selectbox("Pilih simbol unsur:", df["symbol"].sort_values())
        elemen = df[df["symbol"] == selected_symbol].iloc[0]
        st.subheader(f"🔬 {elemen['name']} ({elemen['symbol']})")
        st.markdown(f"**Nomor Atom:** {elemen['number']}")
        st.markdown(f"**Golongan:** {elemen.get('group', '-')}")
        st.markdown(f"**Elektronegativitas:** {elemen.get('electronegativity', '-')}")
        st.markdown(f"**Massa Atom:** {elemen.get('atomic_mass', '-')}")
        st.markdown(f"**Konfigurasi Elektron:** `{elemen.get('electron_configuration', '-')}`")

        if elemen.get('electronegativity') is not None:
            if elemen['electronegativity'] > 3.0:
                st.info("🔸 Unsur sangat elektronegatif")
            elif elemen['electronegativity'] > 2.0:
                st.info("🔹 Unsur elektronegatif")
            else:
                st.info("🔹 Unsur elektropositif")

    with col2:
        st.subheader("📈 Tren Elektronegativitas")
        fig = create_electronegativity_chart(df)
        if fig:
            st.pyplot(fig)

        df_filtered = df.dropna(subset=['electronegativity'])
        if not df_filtered.empty:
            max_en = df_filtered.loc[df_filtered['electronegativity'].idxmax()]
            min_en = df_filtered.loc[df_filtered['electronegativity'].idxmin()]
            st.markdown(f"**Elektronegativitas Tertinggi:** {max_en['name']} ({max_en['electronegativity']})")
            st.markdown(f"**Elektronegativitas Terendah:** {min_en['name']} ({min_en['electronegativity']})")

# ---------------- TAB 2 ----------------
with tab2:
    st.header("🧬 Analisis Struktur Molekul 3D")
    
    col1, col2 = st.columns([1, 1])

    with col1:
        input_method = st.radio("Metode input:", ["Pilih dari contoh", "Input SMILES manual"])
        
        if input_method == "Pilih dari contoh":
            molecule_options = list(MOLECULES.keys())
            molecule_labels = [f"{MOLECULES[smiles]['name']} ({MOLECULES[smiles]['formula']})" for smiles in molecule_options]
            selected_idx = st.selectbox("Pilih molekul:", range(len(molecule_labels)), format_func=lambda x: molecule_labels[x])
            smiles = molecule_options[selected_idx]
            mol_info = MOLECULES[smiles]
            st.info(f"**Deskripsi:** {mol_info['description']}")
        else:
            smiles = st.text_input("Masukkan SMILES:", value="CCO")

        if st.button("🔍 Analisis & Visualisasi 3D", type="primary"):
            if smiles:
                analysis = analyze_smiles(smiles)
                st.success("✅ Analisis berhasil!")
                
                # Info molekul
                if smiles in MOLECULES:
                    mol_info = MOLECULES[smiles]
                    st.markdown(f"**Nama:** {mol_info['name']}")
                    st.markdown(f"**Formula:** {mol_info['formula']}")
                    st.markdown(f"**Berat Molekul:** {mol_info['mw']} g/mol")
                st.markdown(f"**SMILES:** `{smiles}`")

                # Analisis komposisi
                st.subheader("🔬 Komposisi Unsur")
                atom_counts = []
                for elem, count in analysis.items():
                    if elem in ['carbon', 'oxygen', 'nitrogen', 'sulfur', 'phosphorus'] and count > 0:
                        atom_counts.append(f"{elem.capitalize()}: {count}")
                        st.markdown(f"- {elem.capitalize()}: {count} atom")

                # Fitur struktural
                st.subheader("🏗️ Fitur Struktural")
                st.markdown("• " + ("✅ Mengandung struktur cincin" if analysis['has_ring'] else "➖ Struktur rantai terbuka"))
                st.markdown("• " + ("✅ Mengandung ikatan rangkap" if analysis['has_double_bond'] else "➖ Hanya ikatan tunggal"))
                st.markdown("• " + ("✅ Kemungkinan aromatik" if analysis['aromatic'] else "➖ Non-aromatik"))
                
            else:
                st.error("❌ Silakan masukkan SMILES yang valid")

    with col2:
        st.subheader("📚 Panduan SMILES")
        with st.expander("🔤 Notasi Dasar"):
            st.markdown("""
            - C = Karbon  
            - O = Oksigen  
            - N = Nitrogen  
            - S = Sulfur  
            - P = Fosfor  
            - = = Ikatan rangkap dua  
            - # = Ikatan rangkap tiga
            """)
        with st.expander("🔄 Cincin"):
            st.markdown("""
            Gunakan angka untuk menutup cincin:  
            - C1CCCCC1 = Sikloheksana  
            - C1=CC=CC=C1 = Benzena
            """)
        with st.expander("🧪 Contoh Molekul"):
            st.markdown("""
            - CCO = Etanol  
            - O = Air  
            - C = Metana  
            - C1=CC=CC=C1 = Benzena
            """)

    # Visualisasi 3D (full width)
    if 'smiles' in locals() and smiles and st.session_state.get('show_3d', False):
        st.markdown("---")
        molecule_name = MOLECULES.get(smiles, {}).get('name', 'Molekul')
        show_molecule_3d(smiles, molecule_name)

# Tombol untuk toggle 3D viewer
if st.button("🔬 Tampilkan/Sembunyikan Visualisasi 3D"):
    st.session_state['show_3d'] = not st.session_state.get('show_3d', False)
    st.rerun()

# Footer
st.markdown("---")
st.markdown("🧪 **ChemExplorer 3D** - Eksplorasi Kimia dengan Visualisasi 3D")
st.markdown("*Dibuat dengan ❤️ menggunakan Streamlit & py3Dmol*")

# Instalasi info
if not STMOL_AVAILABLE:
    st.sidebar.markdown("## 📦 Instalasi")
    st.sidebar.code("pip install stmol", language="bash")
    st.sidebar.markdown("Restart aplikasi setelah instalasi.")
