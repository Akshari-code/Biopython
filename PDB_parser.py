from Bio.PDB import PDBList, PDBParser
import warnings
from Bio import BiopythonWarning

# Suppress minor PDB header warnings for cleaner console output
warnings.simplefilter('ignore', BiopythonWarning)

def inspect_3d_structure(pdb_id):
    print(f"Downloading 3D structure for PDB ID: {pdb_id}...")

    # 1. Initialize PDBList and download the structure file in standard 'pdb' format
    pdbl = PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")

    # 2. Initialize the parser and load the structure
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, pdb_file_path)

    # 3. Traverse the SMCRA hierarchy (Structure -> Model -> Chain -> Residue)
    model = structure[0]  # Select the primary 3D model
    print(f"\nStructure '{pdb_id}' successfully loaded!")
    print(f"Number of Polypeptide Chains: {len(model)}")
    print("-" * 50)

    for chain in model:
        # Filter out water molecules ('W') and heteroatoms to count pure amino acids
        amino_acids = [res for res in chain if res.id[0] == " "]
        print(f"Chain {chain.id}: Contains {len(amino_acids)} amino acid residues")

        # Inspect the first residue and its Alpha-Carbon ('CA') 3D coordinates
        if amino_acids:
            first_res = amino_acids[0]
            if 'CA' in first_res:
                ca_atom = first_res['CA']
                coords = ca_atom.get_coord()
                print(f"   -> First Residue: {first_res.get_resname()} (Pos {first_res.id[1]}) | Alpha-Carbon XYZ: {coords}")

if __name__ == "__main__":
    try:
        # '6WZG' is Human secretin hormone
        inspect_3d_structure("6WZG")
    except Exception as e:
        print(f"An error occurred: {e}")
