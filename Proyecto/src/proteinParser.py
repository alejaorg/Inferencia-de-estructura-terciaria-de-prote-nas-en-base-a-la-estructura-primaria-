import os
from Bio.PDB import PDBParser
from protein import Protein

parser = PDBParser()

cur_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(cur_dir, "..", "data", "AF-A0A2K6V5L6-F1-model_v6.pdb")
structure = parser.get_structure("AF-A", data_dir)

for model in structure:
    print("Modelo:", model.id)
    
    for chain in model:
        print("  Cadena:", chain.id)
        
        for residue in chain:
            print("    Residuo:", residue.get_resname(), residue.id)
            
            for atom in residue:
                print("      Átomo:", atom.get_name(), atom.coord)
    lista_residuos = list(chain.get_residues())
    print(f"  Total de residuos leídos: {len(lista_residuos)}")
    total_atomos = sum(1 for atom in chain.get_atoms())
    print(f"  Total de átomos extraídos: {total_atomos}")


