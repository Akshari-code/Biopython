import streamlit as st
from bio_engine import (
    analyze_sequence_metrics,
    extract_coding_sequence,
    fetch_genbank_record,
    parse_fasta_text,
)

# 1. Configure webpage title and layout
st.set_page_config(page_title="BioScout | Gene & Mutation Analyzer", page_icon="🧬", layout="wide")

st.title("🧬 BioScout: Automated Gene & Mutation Impact Analyzer")
st.markdown("Fetch live sequences from **NCBI** or upload a local **FASTA** file to analyze coding regions and simulate mutations.")

# Cache NCBI fetches so repeated clicks load instantaneously
@st.cache_data(show_spinner=False)
def cached_ncbi_fetch(acc_id):
    record = fetch_genbank_record(acc_id)
    cds = extract_coding_sequence(record)
    return record.id, record.description, record.annotations.get("organism", "Unknown"), str(cds)

# 2. Sidebar Input Controls
st.sidebar.header("1. Sequence Input Source")
input_mode = st.sidebar.radio("Choose Input Method:", ["NCBI Accession ID", "Upload FASTA File"])

seq_id = None
seq_desc = ""
organism = "Custom / Uploaded"
active_dna = None

try:
    if input_mode == "NCBI Accession ID":
        acc_input = st.sidebar.text_input("Enter NCBI Accession ID:", value="NM_000518", help="Try NM_000518 (Human HBB) or NM_007294 (BRCA1)")
        if acc_input:
            with st.spinner(f"Fetching {acc_input} from NCBI..."):
                seq_id, seq_desc, organism, active_dna = cached_ncbi_fetch(acc_input)

    elif input_mode == "Upload FASTA File":
        uploaded_file = st.sidebar.file_uploader("Upload a .fasta or .fa file", type=["fasta", "fa", "txt"])
        if uploaded_file is not None:
            raw_text = uploaded_file.getvalue().decode("utf-8")
            record = parse_fasta_text(raw_text)
            seq_id = record.id
            seq_desc = record.description
            active_dna = str(record.seq)

except ValueError as e:
    st.error(f"⚠️ {e}")

# 3. Main Dashboard Display
if active_dna:
    st.success(f"**Loaded Sequence:** `{seq_id}` | **Organism:** *{organism}*")
    if seq_desc:
        st.caption(seq_desc)

    # Compute metrics from bio_engine
    metrics = analyze_sequence_metrics(active_dna)

    # Display clean KPI cards in 3 columns
    col1, col2, col3 = st.columns(3)
    col1.metric("Coding Sequence Length", f"{metrics['length']} bp")
    col2.metric("GC Content", f"{metrics['gc_content']}%")
    col3.metric("Translated Protein Length", f"{metrics['length'] // 3} aa")

    # Expandable raw sequence viewer
    with st.expander("View Raw Coding DNA Sequence"):
        st.code(active_dna, language="text")
else:
    st.info("👈 Select an input source in the left sidebar to load a genetic sequence.")
