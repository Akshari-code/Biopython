from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
Entrez.email = "dm.akshari20@gmail.com"

def fetch_genbank_record(accession_id):
    """Fetches a GenBank record from NCBI and returns the SeqRecord."""
    with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
    return record
def extract_coding_sequence(record):
    """Extracts the first Coding Sequence (CDS) DNA slice from a GenBank SeqRecord."""
    for feature in record.features:
        if feature.type == "CDS":
            return feature.extract(record.seq)
    # Fallback to full sequence if no explicit CDS annotation is present
    return record.seq
def analyze_sequence_metrics(seq_obj):
    """Returns length, GC content %, and A/T/G/C counts as a dictionary."""
    seq_str = str(seq_obj).upper()
    length = len(seq_str)
    gc_percent = round(gc_fraction(seq_obj) * 100, 2)

    base_counts = {
        "A": seq_str.count("A"),
        "T": seq_str.count("T"),
        "G": seq_str.count("G"),
        "C": seq_str.count("C")
    }

    return {
        "length": length,
        "gc_content": gc_percent,
        "base_counts": base_counts
    }

def simulate_point_mutation(dna_seq_str, position_1based, new_base):
    """Mutates a single base (1-based index) and classifies the protein impact."""
    dna_clean = str(dna_seq_str).upper()
    new_base = new_base.upper()
    idx = position_1based - 1
    
    if idx < 0 or idx >= len(dna_clean):
        raise ValueError(f"Position {position_1based} is out of bounds for sequence length {len(dna_clean)}.")
    if new_base not in ["A", "T", "G", "C"]:
        raise ValueError("New base must be A, T, G, or C.")
        
    original_base = dna_clean[idx]
    mutated_dna = dna_clean[:idx] + new_base + dna_clean[idx + 1:]
    
    wt_protein = str(Seq(dna_clean).translate(to_stop=True))
    mut_protein = str(Seq(mutated_dna).translate(to_stop=True))
    
    codon_index = idx // 3
    codon_num = codon_index + 1
    if original_base == new_base:
        impact = "No Change (Identical Base)"
        detail = f"Position {position_1based} remained {original_base}."
    elif len(mut_protein) < len(wt_protein):
        impact = "Nonsense Mutation (Premature Stop)"
        detail = f"Protein truncated at amino acid position {len(mut_protein) + 1}."
    elif wt_protein == mut_protein:
        impact = "Silent Mutation (Synonymous)"
        aa = wt_protein[codon_index] if codon_index < len(wt_protein) else "Stop"
        detail = f"Codon {codon_num} changed ({original_base}->{new_base}), but amino acid remained '{aa}'."
    else:
        impact = "Missense Mutation (Non-synonymous)"
        wt_aa = wt_protein[codon_index]
        mut_aa = mut_protein[codon_index]
        detail = f"Amino acid at position {codon_num} changed from '{wt_aa}' to '{mut_aa}'."
        
    return {
        "position": position_1based,
        "original_base": original_base,
        "mutated_base": new_base,
        "impact": impact,
        "detail": detail,
        "wt_protein": wt_protein,
        "mut_protein": mut_protein
    }

if __name__ == "__main__":
    test_record = fetch_genbank_record("NM_000518")
    cds_seq = extract_coding_sequence(test_record)
    mutation_report = simulate_point_mutation(cds_seq, 20, "T")
    
    print("\n--- Mutation Analysis Report ---")
    print(f"Classification: {mutation_report['impact']}")
    print(f"Summary:        {mutation_report['detail']}")