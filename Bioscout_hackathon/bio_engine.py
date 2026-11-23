from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import re
from Bio.SeqUtils import gc_fraction, nt_search
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
import io
from Bio import SeqIO

def parse_fasta_text(fasta_content):
    """
    Parses raw FASTA formatted text into a Biopython SeqRecord.
    Works seamlessly with both local file reads and Streamlit file uploads.
    """
    fasta_io = io.StringIO(fasta_content.strip())
    records = list(SeqIO.parse(fasta_io, "fasta"))

    if not records:
        raise ValueError("No valid FASTA sequence found. Make sure the text starts with a '>' header line.")

    # Return the first sequence record in the file
    return records[0]
if __name__ == "__main__":
    sample_fasta = """
    >Custom_Test_Gene
    ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG
    """
    custom_record = parse_fasta_text(sample_fasta)
    metrics = analyze_sequence_metrics(custom_record.seq)
    mut_report = simulate_point_mutation(custom_record.seq, 20, "T")

    print(f"Loaded Custom ID: {custom_record.id}")
    print(f"Length: {metrics['length']} bp | GC Content: {metrics['gc_content']}%")
    print(f"Mutation Check: {mut_report['impact']} -> {mut_report['detail']}")
def batch_scan_mutations(dna_seq_str, mutations_list):
    """
    Scans a list of (position_1based, new_base) tuples against a DNA sequence
    and returns a clean Pandas DataFrame summarizing the impact of each variant.
    """
    results = []

    for pos, new_base in mutations_list:
        try:
            report = simulate_point_mutation(dna_seq_str, pos, new_base)
            results.append({
                "Position": report["position"],
                "Original_Base": report["original_base"],
                "Mutated_Base": report["mutated_base"],
                "Impact": report["impact"],
                "Summary": report["detail"]
            })
        except ValueError as e:
            # Record invalid positions cleanly without crashing the whole batch
            results.append({
                "Position": pos,
                "Original_Base": "-",
                "Mutated_Base": new_base,
                "Impact": "Error",
                "Summary": str(e)
            })
    return pd.DataFrame(results)
if __name__ == "__main__":
    test_record = fetch_genbank_record("NM_000518")
    cds_seq = extract_coding_sequence(test_record)

    # Test a batch of mutations: (1-based position, new base)
    sample_variants = [
        (9, "T"),     # Silent mutation in Codon 3 (CAT -> CAC/CAT)
        (20, "T"),    # Missense mutation (Sickle Cell: GAG -> GTG)
        (118, "T"),   # Nonsense mutation (Creates a premature stop codon)
        (999, "G")    # Out-of-bounds error check
    ]

    batch_df = batch_scan_mutations(cds_seq, sample_variants)
    print("\n--- Batch Variant Scan Report ---")
    print(batch_df.to_string(index=False))
def search_motif(dna_seq_str, motif_str):
    """
    Finds all 1-based starting positions of a target motif (including overlaps
    and IUPAC ambiguous codes) inside a DNA sequence.
    """
    dna_clean = str(dna_seq_str).upper().strip()
    motif_clean = str(motif_str).upper().strip()

    if not motif_clean:
        return []

    results = nt_search(dna_clean, motif_clean)
    # results[0] is the motif pattern; results[1:] are 0-based indices
    return [idx + 1 for idx in results[1:]]

def find_all_orfs(dna_seq_str, min_protein_length=10):
    """
    Scans all 6 reading frames (3 forward, 3 reverse complement) for Open Reading Frames
    starting with 'M' and ending with a stop codon, filtered by minimum amino acid length.
    """
    dna_seq = Seq(str(dna_seq_str).upper().strip())
    pattern = re.compile(r'(?=(M[^*]*\*))')
    orf_results = []
    strands = [
        ("Forward (+)", dna_seq),
        ("Reverse (-)", dna_seq.reverse_complement())
    ]

    for strand_label, seq in strands:
        for frame in range(3):
            translated = str(seq[frame:].translate())
            for match in pattern.findall(translated):
                protein = match[:-1]  # Drop trailing '*'
                if len(protein) >= min_protein_length:
                    orf_results.append({
                        "Strand": strand_label,
                        "Frame": frame + 1,
                        "Length_AA": len(protein),
                        "Protein_Sequence": protein
                    })

    return pd.DataFrame(orf_results)
if __name__ == "__main__":
    test_record = fetch_genbank_record("NM_000518")
    cds_seq = extract_coding_sequence(test_record)

    # 1. Test Motif Finder
    motif_hits = search_motif(cds_seq, "CTGAGG")
    print(f"\nMotif 'CTGAGG' found at 1-based positions: {motif_hits}")

    # 2. Test 6-Frame ORF Scanner
    orfs_df = find_all_orfs(cds_seq, min_protein_length=20)
    print("\n--- Detected Open Reading Frames (>= 20 AA) ---")
    print(orfs_df.to_string(index=False))

