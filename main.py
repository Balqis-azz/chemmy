import streamlit as st
import pandas as pd
import plotly.express as px
import json
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

st.set_page_config(page_title="ChemExplorer 2.0", layout="wide")

# Load data unsur
with open("periodic_table.json") as file:
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
    st.markdown(f"*Nomor Atom:* {elemen['number']}")
    st.markdown(f"*Golongan:* {elemen.get('group', '-')}")
    st.markdown(f"*Elektronegativitas:* {elemen.get('electronegativity', '-')}")
    st.markdown(f"*Konfigurasi Elektron:* {elemen.get('electron_configuration', '-')}")

    st.subheader("📈 Tren Elektronegativitas")
    fig = px.line(df.sort_values("number"), x="number", y="electronegativity",
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
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                st.success("Struktur dikenali!")
                st.image(Draw.MolToImage(mol, size=(300, 300)))
            else:
                st.error("Struktur tidak valid atau tidak dikenali.")

    with col2:
        st.markdown("""
        #### Contoh SMILES:
        - *Benzena*: C1=CC=CC=C1
        - *Ethanol*: CCO
        - *Glukosa*: C(C1C(C(C(C(O1)O)O)O)O)O
        - *Asam asetat*: CC(=O)O
        """)

---

## 🚀 LANGKAH 4 – Jalankan Aplikasinya

```bash
streamlit run main.py
