from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

def find_binding_pocket_residues(pdb_file, chain_id="A", ligand_name="HEM", cutoff=4.0):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    # FIX 1: Access Model 0 first, THEN access Chain 'A'
    model = structure[0]
    chain = model[chain_id]
    
    ligand_atoms = []
    protein_residues = []
    
    # 1. Separate the target ligand atoms from the standard amino acid residues
    for residue in chain:
        # FIX 2: residue.id is a tuple: (hetfield, resseq, icode)
        flag = residue.id[0]
        if flag == " ":
            protein_residues.append(residue)
        elif flag == f"H_{ligand_name}":
            ligand_atoms.extend(residue.get_atoms())
            
    if not ligand_atoms:
        print(f"Ligand {ligand_name} not found in Chain {chain_id}.")
        return []

    print(f"Found {len(ligand_atoms)} atoms in ligand {ligand_name}. Scanning for nearby residues (< {cutoff} Å)...")
    print("-" * 60)
    
    # Use a dictionary to store the minimum distance per residue
    pocket_residues = {}
    
    # 2. Measure distance between every protein atom and every ligand atom
    for res in protein_residues:
        res_key = (res.get_resname(), res.id[1])  # (Amino Acid Name, Sequence Number)
        for p_atom in res:
            for l_atom in ligand_atoms:
                distance = p_atom - l_atom
                if distance <= cutoff:
                    # Keep track of the closest contact distance for this residue
                    if res_key not in pocket_residues or distance < pocket_residues[res_key]:
                        pocket_residues[res_key] = round(distance, 2)
                    
    # Sort results by residue sequence position
    sorted_pocket = [(name, pos, dist) for (name, pos), dist in sorted(pocket_residues.items(), key=lambda x: x[0][1])]
    return sorted_pocket

if __name__ == "__main__":
    try:
        nearby = find_binding_pocket_residues("pdb1a3n.ent", chain_id="A", ligand_name="HEM", cutoff=4.0)
        print(f"Active Pocket Residues interacting with Heme in Chain A ({len(nearby)} found):")
        for res_name, pos, dist in nearby:
            print(f"  -> {res_name} (Position {pos}) | Closest contact: {dist} Å")
    except FileNotFoundError:
        print("Please ensure 'pdb1a3n.ent' from Day 14 is in your current directory.")
    except KeyError as e:
        print(f"Could not find key {e}. Available chains in Model 0: {[c.id for c in PDBParser().get_structure('p', 'pdb1a3n.ent')[0]]}")
