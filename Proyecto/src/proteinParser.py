from Bio.PDB import PDBParser
from protein import Protein

parser = PDBParser()

structure = parser.get_structure("AF-A", "../../data/AF-A0A2K6V5L6-F1-model_v6.pdb")

for model in structure:
    print("Modelo:", model.id)
    
    for chain in model:
        print("  Cadena:", chain.id)
        
        for residue in chain:
            print("    Residuo:", residue.get_resname(), residue.id)
            
            for atom in residue:
                print("      Átomo:", atom.get_name(), atom.coord)


