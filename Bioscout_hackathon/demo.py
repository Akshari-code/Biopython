import streamlit as st

st.title("Mini DNA Analyzer")
dna = st.text_input("Enter DNA Sequence:", "ATGCGTAC")
st.write("Sequence Length:", len(dna))
st.bar_chart({"A": dna.count("A"), "T": dna.count("T"), "G": dna.count("G"), "C": dna.count("C")})
