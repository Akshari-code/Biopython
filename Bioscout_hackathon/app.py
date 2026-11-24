import streamlit as st
from bio_engine import (
    analyze_sequence_metrics,
    extract_coding_sequence,
    fetch_genbank_record,
    parse_fasta_text,)
import base64
import streamlit as st
from bio_engine import (
    analyze_sequence_metrics,
    batch_scan_mutations,
    extract_coding_sequence,
    fetch_genbank_record,
    parse_fasta_text,
    simulate_point_mutation,)


# 1. Configure webpage title and layout
st.set_page_config(page_title="BioScout | Gene & Mutation Analyzer", layout="wide")

st.title("BioScout: Automated Gene & Mutation Impact Analyzer")
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


if active_dna:
    st.success(f"**Loaded Sequence:** `{seq_id}` | **Organism:** *{organism}*")
    if seq_desc:
        st.caption(seq_desc)

    metrics = analyze_sequence_metrics(active_dna)

    # Top KPI Cards
    col1, col2, col3 = st.columns(3)
    col1.metric("Coding Sequence Length", f"{metrics['length']} bp")
    col2.metric("GC Content", f"{metrics['gc_content']}%")
    col3.metric("Translated Protein Length", f"{metrics['length'] // 3} aa")

    st.divider()

    # Create Interactive Workspace Tabs
    tab1, tab2, tab3 = st.tabs([
        "📊 Nucleotide Composition",
        "🧬 Point Mutation Simulator",
        "📋 Batch Variant Scanner"
    ])

    # --- TAB 1: SEQUENCE OVERVIEW & CHART ---
    with tab1:
        st.subheader("Base Pair Distribution")
        st.bar_chart(metrics["base_counts"])
        with st.expander("View Raw Coding DNA Sequence"):
            st.code(active_dna, language="text")

    # --- TAB 2: SINGLE POINT MUTATION SIMULATOR ---
    with tab2:
        st.subheader("Simulate a Single Nucleotide Polymorphism (SNP)")
        c1, c2 = st.columns(2)

        with c1:
            pos_input = st.number_input(
                "Select 1-Based DNA Position:",
                min_value=1,
                max_value=len(active_dna),
                value=min(20, len(active_dna)),
                step=1
            )
        with c2:
            current_base = active_dna[pos_input - 1]
            st.info(f"Reference Base at Position **{pos_input}**: **`{current_base}`**")
            new_base_input = st.selectbox("Select Mutated Base:", ["A", "T", "G", "C"], index=1)

        if st.button("Run Mutation Simulation", type="primary"):
            try:
                report = simulate_point_mutation(active_dna, pos_input, new_base_input)

                # Color-coded alert based on severity
                if "Nonsense" in report["impact"] or "Missense" in report["impact"]:
                    st.warning(f"**{report['impact']}** — {report['detail']}")
                elif "Silent" in report["impact"]:
                    st.info(f"**{report['impact']}** — {report['detail']}")
                else:
                    st.success(f"**{report['impact']}** — {report['detail']}")

                with st.expander("Compare Wild-Type vs. Mutated Protein Sequence"):
                    st.text("Wild-Type Protein:")
                    st.code(report["wt_protein"], language="text")
                    st.text("Mutated Protein:")
                    st.code(report["mut_protein"], language="text")
            except ValueError as e:
                st.error(str(e))

    # --- TAB 3: BATCH VARIANT SCANNER ---
    with tab3:
        st.subheader("Scan Multiple Variants Simultaneously")
        st.markdown("Enter one mutation per line in the format **`Position,NewBase`** (e.g., `20,T`):")

        default_batch = "9,T\n20,T\n118,T"
        batch_text = st.text_area("Variant List:", value=default_batch, height=120)

        if st.button("Scan All Variants"):
            parsed_variants = []
            for line in batch_text.strip().splitlines():
                if "," in line:
                    parts = line.split(",")
                    try:
                        parsed_variants.append((int(parts.strip()), parts.strip().upper()))
                    except ValueError:
                        st.warning(f"Skipped invalid line: `{line}`")

            if parsed_variants:
                results_df = batch_scan_mutations(active_dna, parsed_variants)
                st.dataframe(results_df, use_container_width=True)

                # One-click CSV export
                csv_data = results_df.to_csv(index=False).encode("utf-8")
                st.download_button(
                    label="📥 Download Clinical Variant Report (.csv)",
                    data=csv_data,
                    file_name=f"{seq_id}_variant_report.csv",
                    mime="text/csv")
