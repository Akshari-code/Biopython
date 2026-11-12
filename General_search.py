from Bio import Entrez
from Bio import SeqIO

Entrez.email = "dm.akshari20@gmail.com"

def search_and_fetch_metadata(search_term, max_results=3):
    print(f"Searching NCBI Nucleotide database for: '{search_term}'...")

    # 1. Run the text search query
    with Entrez.esearch(db="nucleotide", term=search_term, retmax=max_results) as search_handle:
        search_results = Entrez.read(search_handle)

    # 2. Extract the list of matching IDs found
    id_list = search_results["IdList"]
    print(f"Found {search_results['Count']} total matches. Fetching details for top {len(id_list)} IDs...")

    # 3. Use the IDs to fetch minimal structural record details
    for ncbi_id in id_list:
        with Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text") as fetch_handle:
            record = SeqIO.read(fetch_handle, "fasta")
            print("-" * 50)
            print(f"Accession ID: {record.id}")
            print(f"Description:  {record.description}")
            print(f"Sequence Len: {len(record.seq)} bp")

if __name__ == "__main__":
    query = "insulin"
    try:
        search_and_fetch_metadata(query)
    except Exception as e:
        print(f"An error occurred: {e}")
