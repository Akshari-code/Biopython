from Bio import Entrez
from Bio import SeqIO

Entrez.email = "dm.akshari20@gmail.com" 

def inspect_genbank_record(accession_id):
    print(f"Downloading full GenBank record for {accession_id}...")

    # 1. Fetch in 'gb' (GenBank) format 
    with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    print(f"\nOrganism: {record.annotations.get('organism')}")
    print(f"Total Sequence Length: {len(record.seq)} bp")
    print(f"Total Annotated Features: {len(record.features)}")
    print("-" * 50)

    # 2. Loop through biological features to find the Coding Sequence (CDS)
    for feature in record.features:
        if feature.type == "CDS":
            print(f"Found Feature: {feature.type}")
            print(f"Exact Location: {feature.location}")

            # Extract the coding DNA slice using the feature location
            cds_dna = feature.extract(record.seq)
            print(f"Coding DNA Length: {len(cds_dna)} bp")

            # Access official NCBI qualifiers (metadata dictionary)
            product_name = feature.qualifiers.get("product", ["Unknown"])[0]
            protein_seq = feature.qualifiers.get("translation", ["Not provided"])[0]

            print(f"Protein Product: {product_name}")
            print(f"Official Protein Sequence:\n{protein_seq}")

if __name__ == "__main__":
    try:
        inspect_genbank_record("NM_000518")
    except Exception as e:
        print(f"Error fetching GenBank record: {e}")
