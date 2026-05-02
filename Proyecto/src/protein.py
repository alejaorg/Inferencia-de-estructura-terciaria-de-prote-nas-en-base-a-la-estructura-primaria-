import numpy as np
from aminoacid_translation import ( Atoms_composition, Aminoacid_map )
from experimental_data import (
    BOND_LENGTHS,
    BOND_ANGLES,
    RAMACHANDRAN,
    OMEGA
)

from experimental_data import GEOMETRY
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

            
        for res in self.residues:
            build_sidechain(res)
            


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

def sample_rotamer():
    angles = list(GEOMETRY["rotamer_probs"].keys())
    probabilities = list(GEOMETRY["rotamer_probs"].values())
    
    return random.choices(angles, weights=probabilities)[0]

def place_CB(N, CA, C):
    return place_atom(N, CA, C, BOND_LENGTHS[("CA","CB")], 109.5, 122.5)

def build_sidechain(residue):
    
    name = residue.name
    
    if name == "GLY":
        return
    
    N = residue.get_atom("N").coord
    CA = residue.get_atom("CA").coord
    C = residue.get_atom("C").coord
    
    CB = place_CB(N, CA, C)
    residue.get_atom("CB").set_coor(CB)
    
    atom_coords = {
        "N": N,
        "CA": CA,
        "CB": CB
    }
    
    chi_defs = GEOMETRY["chi_definitions"].get(name, [])
    
    for chi in chi_defs:
        
        a_name, b_name, c_name, d_name = chi
        
        a = atom_coords[a_name]
        b = atom_coords[b_name]
        c = atom_coords[c_name]
        
        chi_angle = sample_rotamer()
        
        length = BOND_LENGTHS.get((c_name, d_name), 1.53)
        
        d = place_atom(a, b, c, length, 109.5, chi_angle)
        
        residue.get_atom(d_name).set_coor(d)
        atom_coords[d_name] = d
        
    atoms = residue.get_atoms()
    
    for atom in atoms:
        if np.all(atom.coord == 0.0):
            
            parent = infer_parent(atom.name)
            
            if parent is None or parent not in atom_coords:
                continue
            
            ref = choose_reference_atoms(atom_coords, parent)
            
            if ref is None:
                continue
            
            a, b, c = ref
            
            length = BOND_LENGTHS.get((parent, atom.name), 1.53)
            
            d = place_atom(a, b, c, length, 109.5, sample_rotamer())
            
            atom.set_coor(d)
            atom_coords[atom.name] = d
            
                   
        
def infer_parent(atom_name):

    hierarchy = {
        "CG2": "CB",
        "CD1": "CG",
        "CD2": "CG",
        "CE1": "CD1",
        "CE2": "CD2",
        "CZ": "CE1",
        "NZ": "CE",
        "OE1": "CD",
        "OE2": "CD",
        "OD1": "CG",
        "OD2": "CG",
        "NE": "CD",
        "NH1": "CZ",
        "NH2": "CZ"
    }

    return hierarchy.get(atom_name, None)

def choose_reference_atoms(atom_coords, parent):

    if parent == "CB":
        if "N" in atom_coords and "CA" in atom_coords:
            return atom_coords["N"], atom_coords["CA"], atom_coords["CB"]

    if parent == "CG":
        if "CA" in atom_coords and "CB" in atom_coords:
            return atom_coords["CA"], atom_coords["CB"], atom_coords["CG"]

    if parent == "CD":
        if "CB" in atom_coords and "CG" in atom_coords:
            return atom_coords["CB"], atom_coords["CG"], atom_coords["CD"]

    # fallback seguro
    keys = list(atom_coords.keys())
    if len(keys) >= 3:
        return (
            atom_coords[keys[-3]],
            atom_coords[keys[-2]],
            atom_coords[keys[-1]]
        )

    return None
    
    



p = Protein("GASV")
p.init_random_structure()   

for res in p.residues:
    print(res)
    for atom in res.get_atoms():
        print(atom)
        
        
        
    
        
    