import numpy as np
from aminoacid_translation import ( Atoms_composition, Aminoacid_map )
from experimental_data import (
    BOND_LENGTHS,
    BOND_ANGLES,
    RAMACHANDRAN,
    OMEGA
)
import random

class Atom:
    def __init__(self, name, coord):
        self.name = name
        self.coord = np.array(coord, dtype=float)
        
    def set_coor(self, new_coord):
        self.coord = np.array(new_coord, dtype=float)
        
    def __repr__(self):
        return f"{self.name}: {self.coord}"
    
    
class Residue:
    
    def __init__(self, name, index):
        self.name = name
        self.index = index
        self.atoms = {}
        
    def add_atom(self, atom):
        self.atoms[atom.name] = atom
        
    def get_atom(self, name):
        return self.atoms.get(name)
    
    def get_atoms(self):
        return list(self.atoms.values())
    
    def __repr__(self):
        return f"{self.name} {self.index}"

class Protein: 
    
    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.residues = []
        self.build_residue()
        
    def build_residue(self):
        for i, aa in enumerate(self.sequence, start=1):
            res_name = Aminoacid_map[aa]
            residue = Residue(res_name, i)
            
            atom_names = Atoms_composition[res_name]
            
            for atom_name in atom_names:
                residue.add_atom(Atom(atom_name, [0.0, 0.0, 0.0]))
            
            self.residues.append(residue)
            
    def get_residue(self, i):
        return self.residues[i]
    
    def __len__(self):
        return len(self.residues)
    
    def __repr__(self):
        return f"Protein(n_residues={len(self)})"
    
    def init_random_structure(self):
        N = np.array([0.0, 0.0, 0.0])
        CA = np.array([BOND_LENGTHS[("N", "CA")], 0.0, 0.0])
        
        angle = np.radians(BOND_ANGLES[("N", "CA", "C")])
        C = np.array([CA[0] + BOND_LENGTHS[("CA", "C")] * np.cos(angle), BOND_LENGTHS[("CA", "C")] * np.sin(angle), 0.0])
        
        self.residues[0].get_atom("N").set_coor(N)
        self.residues[0].get_atom("CA").set_coor(CA)
        self.residues[0].get_atom("C").set_coor(C)
        
        O = place_atom(N, CA, C, BOND_LENGTHS[("C", "O")], BOND_ANGLES[("CA","C","O")], 0.0)
        self.residues[0].get_atom("O").set_coor(O)
        prev = [N, CA, C]
        
        for i in range(1, len(self.residues)):
            res = self.residues[i]
            
            phi, psi = sample_phi_psi()
            omega = OMEGA["trans"]["angle"]
            
            new_N = place_atom(prev[0], prev[1], prev[2], BOND_LENGTHS[("C", "N")], BOND_ANGLES[("CA", "C", "N")], omega)
            new_CA = place_atom(prev[1], prev[2], new_N, BOND_LENGTHS[("N", "CA")], BOND_ANGLES[("C", "N", "CA")], phi)
            new_C = place_atom(prev[2], new_N, new_CA, BOND_LENGTHS[("CA", "C")], BOND_ANGLES[("N", "CA", "C")], psi)
            new_O = place_atom(new_N, new_CA, new_C,  BOND_LENGTHS[("C", "O")], BOND_ANGLES[("CA","C","O")], 0.0)
            
            res.get_atom("N").set_coor(new_N)
            res.get_atom("CA").set_coor(new_CA)
            res.get_atom("C").set_coor(new_C)
            res.get_atom("O").set_coor(new_O)
            
            prev = [new_N, new_CA, new_C]
            


def normalize(vector):
    return vector / np.linalg.norm(vector)

def place_atom(a, b, c, length, angle_deg, dihedral_deg):
    angle = np.radians(angle_deg)
    dihedral = np.radians(dihedral_deg)
    
    bc = normalize(c - b)
    n = normalize(np.cross(b - a, bc))
    m = np.cross(n, bc)
    
    d = (c + length * (-np.cos(angle) * bc + np.sin(angle) * (np.cos(dihedral) * m + np.sin(dihedral) * n)))
    
    return d

def sample_phi_psi():
    keys = list(RAMACHANDRAN.keys())
    weights = [RAMACHANDRAN[k]["weight"] for k in keys]
    
    region = random.choices(keys, weights=weights)[0]
    r = RAMACHANDRAN[region]
    
    phi = random.uniform(*r["phi"])
    psi = random.uniform(*r["psi"])
    
    return phi, psi 



p = Protein("GASV")
p.init_random_structure()   

for res in p.residues:
    print(res)
    for atom in res.get_atoms():
        print(atom)
        
        
        
    
        
    