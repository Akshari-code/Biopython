import pandas as pd
from Bio import Entrez
from Bio import SeqIO

# 1. NCBI Identity Setup (Mandatory)
Entrez.email = "dm.akshari20@gmail.com" 

def fetch_sequence_from_ncbi(accession_id):
    print(f"Connecting to NCBI to fetch: {accession_id}...")

    # 2. Fetch data via the Entrez API
    with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
        # 3. Parse the live stream directly with SeqIO
        record = SeqIO.read(handle, "fasta")

    return record

if __name__ == "__main__":
    # Accession number for Human Hemoglobin Subunit Beta (HBB) mRNA
    hbb_accession = "NM_000518"

    try:
        hbb_record = fetch_sequence_from_ncbi(hbb_accession)

        print(f"\nSuccessfully Retrieved Sequence!")
        print(f"ID: {hbb_record.id}")
        print(f"Description: {hbb_record.description}")
        print(f"Sequence Length: {len(hbb_record.seq)} base pairs")

        # Calculate composition for visual summary
        seq_str = str(hbb_record.seq).upper()
        composition = {base: seq_str.count(base) for base in ["A", "T", "G", "C"]}
        print("\nNucleotide Base Counts:", composition)

    except Exception as e:
        print(f"An error occurred during database retrieval: {e}")
