import numpy as np
from aminoacid_translation import ( Atoms_composition, Aminoacid_map )
from experimental_data import (
    BOND_LENGTHS,
    BOND_ANGLES,
    RAMACHANDRAN,
    OMEGA,
    VAN_DER_WAALS_RADIUS,
    VAN_DER_WAALS_TOLERANCE
)

from experimental_data import GEOMETRY

from itertools import product
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
        
        try_build_sidechain(self.residues[0], [])
        built_residues = [self.residues[0]]
        prev = [N, CA, C]

        for i in range(1, len(self.residues)):
            res = self.residues[i]

            prev, success, level = build_residue_clashes(res, built_residues, prev)
            
            if not success:
                print(f"{res}: no se encontró conformación sin clashes") #TODO
            else:
                print(f"{res}: aceptado via {level}")
             
            built_residues.append(res)    

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
    return place_atom(C, N, CA, BOND_LENGTHS[("CA","CB")], 109.5, 122.5)

def place_aromatic_ring(residue, atoms, chi2):
    CB = atoms["CB"]
    CG = atoms["CG"]
    CA = atoms["CA"]
        
    if residue.name == "HIS":
        CD2 = place_atom(CA, CB, CG,
                     BOND_LENGTHS.get(("CG", "CD2")), 122.0, 
                     chi2 + 180.0)
        residue.get_atom("CD2").set_coor(CD2)
        atoms["CD2"] = CD2
        
        ND1 = atoms["ND1"]
        CE1 = place_atom(CB, CG, ND1,
                        BOND_LENGTHS.get(("CD1", "CE1")), BOND_ANGLES.get(("CG", "CD1", "CE1")),
                        180.0)
        residue.get_atom("CE1").set_coor(CE1)
        atoms["CE1"] = CE1
        NE2 = place_atom(CG, ND1, CE1,
                        BOND_LENGTHS.get(("CE1", "NE2")), BOND_ANGLES.get(("ND1", "CE1", "NE2")),
                        0.0)
        residue.get_atom("NE2").set_coor(NE2)
        atoms["NE2"] = NE2
    
    else:    
        CD2 = place_atom(CA, CB, CG,
                     BOND_LENGTHS.get(("CG", "CD2")), BOND_ANGLES.get(("CB", "CG", "CD2")), 
                     chi2 + 180.0)
        residue.get_atom("CD2").set_coor(CD2)
        atoms["CD2"] = CD2
        
        CD1 = atoms["CD1"]
        CE1 = place_atom(CB, CG, CD1,
                        BOND_LENGTHS.get(("CD1", "CE1")), BOND_ANGLES.get(("CG", "CD1", "CE1")),
                        180.0)
        residue.get_atom("CE1").set_coor(CE1)
        atoms["CE1"] = CE1
        
        CE2 = place_atom(CB, CG, CD2,
                        BOND_LENGTHS.get(("CD2", "CE2")), BOND_ANGLES.get(("CG", "CD2", "CE2")),
                        180.0)
        residue.get_atom("CE2").set_coor(CE2)
        atoms["CE2"] = CE2
        
        CZ = place_atom(CG, CD1, CE1,
                    BOND_LENGTHS.get(("CE1", "CZ")), BOND_ANGLES.get(("CD1", "CE1", "CZ")),
                    0.0)
        residue.get_atom("CZ").set_coor(CZ)
        atoms["CZ"] = CZ

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
        "CB": CB,
        "C": C
    }
    
    chi_defs = GEOMETRY["chi_definitions"].get(name, [])
    chi_angles = []
    
    for chi in chi_defs:
        
        a_name, b_name, c_name, d_name = chi
        
        a = atom_coords[a_name]
        b = atom_coords[b_name]
        c = atom_coords[c_name]
        
        chi_angle = sample_rotamer()
        chi_angles.append(chi_angle)
        
        length = BOND_LENGTHS.get((c_name, d_name), 1.53)
        angle = BOND_ANGLES.get((b_name, c_name, d_name), 109.5)
        
        d = place_atom(a, b, c, length, angle, chi_angle)
        
        residue.get_atom(d_name).set_coor(d)
        atom_coords[d_name] = d
        
    # To place CG2 when the residue has, to not place CG1 & CG2 in the same point 
    if name in ("VAL", "ILE", "THR"):
        cg2_atom = residue.get_atom("CG2")
        if cg2_atom is not None:
            length = BOND_LENGTHS.get(("CB", "CG2"), 1.521)
            cg2 = place_atom(N, CA, CB, length, 110.1, chi_angles[0] + 120.0)
            cg2_atom.set_coor(cg2)
            atom_coords["CG2"] = cg2
            
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
            
def try_build_sidechain(residue, built_residues):
    chi_defs = GEOMETRY["chi_definitions"].get(residue.name, [])
    
    if not chi_defs:
        build_sidechain(residue)
        return not detect_clashes_between(residue, built_residues)
    
    rotamers = GEOMETRY["rotamers"]["default"] 
    
    combinations = list(product(rotamers, repeat=len(chi_defs)))
    
    for comb in combinations:
        N = residue.get_atom("N").coord
        CA = residue.get_atom("CA").coord
        C = residue.get_atom("C").coord
        
        CB = place_CB(N, CA, C)
        residue.get_atom("CB").set_coor(CB)
        
        atoms = {"N": N, "CA": CA, "C": C, "CB": CB}
        
        chi_angles = []
        for chi, angle in zip(chi_defs, comb):
            a_name, b_name, c_name, d_name = chi
            
            a = atoms[a_name]
            b = atoms[b_name]
            c = atoms[c_name]
            
            length = BOND_LENGTHS.get((c_name, d_name), 1.53)
            bond_angle = BOND_ANGLES.get((b_name, c_name, d_name), 109.5)
            
            d = place_atom(a, b, c, length, bond_angle, angle)
            residue.get_atom(d_name).set_coor(d)
            atoms[d_name] = d
            chi_angles.append(angle)
        if residue.name in ("VAL", "ILE", "THR"):
            cg2_atom = residue.get_atom("CG2")
            if cg2_atom is not None:
                length = BOND_LENGTHS.get(("CB", "CG2"), 1.521)
                cg2 = place_atom(N, CA, CB, length, 110.1, chi_angles[0] + 120.0)
                cg2_atom.set_coor(cg2)
                atoms["CG2"] = cg2
        
        if residue.name in ("PHE", "TYR", "HIS"):
            place_aromatic_ring(residue, atoms, chi_angles[1])
        if not detect_clashes_between(residue, built_residues):
            return True
    return False

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
    
def detect_clashes(protein, tolerance=VAN_DER_WAALS_TOLERANCE):

    clashes = []

    atoms = []
    for res in protein.residues:
        for atom in res.get_atoms():
            if not np.all(atom.coord == 0.0):
                atoms.append((res, atom))

    # Compara los pares de atomos
    for i in range(len(atoms)):
        res_i, atom_i = atoms[i]
        for j in range(i + 1, len(atoms)):
            res_j, atom_j = atoms[j]

            # Ignore pairs 1-2 & 1-3 in the same residue
            # By definition, as we build the residuesthis is always a clash
            if res_i.index == res_j.index:
                continue

            # Ignore pairs between consecutive residues as the distance between them is define by the angles omega/phi/psi
            backbone = {"N", "CA", "C", "O"}
            if abs(res_i.index - res_j.index) == 1:
                if atom_i.name in backbone and atom_j.name in backbone:
                    continue

            dist = np.linalg.norm(atom_i.coord - atom_j.coord)

            r_i = VAN_DER_WAALS_RADIUS.get(atom_i.name[0], 1.70)
            r_j = VAN_DER_WAALS_RADIUS.get(atom_j.name[0], 1.70)
            threshold = (r_i + r_j) * tolerance
            #Rosetta use a dinamic algorithm to be permisive with clashes, but it's aprox a 0.8 tolerance with the real Van Der Waals radius
            if dist < threshold:
                clashes.append({
                    "res1":  res_i,
                    "atom1": atom_i.name,
                    "res2":  res_j,
                    "atom2": atom_j.name,
                    "dist":  round(dist, 3),
                    "min_allowed": round(threshold, 3)
                })

    return clashes

def detect_clashes_between(res, other_res, tolerance=VAN_DER_WAALS_TOLERANCE):

    clashes = []
    backbone = {"N", "CA", "C", "O"}
    
    for a in other_res:
        for atom_i in res.get_atoms():
            if np.all(atom_i.coord == 0.0):
                continue
            for atom_j in a.get_atoms():
                if np.all(atom_j.coord == 0.0):
                    continue
                
                if abs(res.index - a.index) == 1:
                    if atom_i.name in backbone and atom_j.name in backbone:
                        continue
    
                dist = np.linalg.norm(atom_i.coord - atom_j.coord)

                r_i = VAN_DER_WAALS_RADIUS.get(atom_i.name[0], 1.70)
                r_j = VAN_DER_WAALS_RADIUS.get(atom_j.name[0], 1.70)
                threshold = (r_i + r_j) * tolerance
                #Rosetta use a dinamic algorithm to be permisive with clashes, but we found an aprox 0.8 tolerance with the real Van Der Waals radius is admisible
                if dist < threshold:
                    clashes.append((atom_i.name, atom_j.name, dist))

    return clashes

def build_backbone(res, prev, phi, psi):
    omega = OMEGA["trans"]["angle"]

    new_N = place_atom(prev[0], prev[1], prev[2], BOND_LENGTHS[("C", "N")], BOND_ANGLES[("CA", "C", "N")], omega)
    new_CA = place_atom(prev[1], prev[2], new_N, BOND_LENGTHS[("N", "CA")], BOND_ANGLES[("C", "N", "CA")], phi)
    new_C = place_atom(prev[2], new_N, new_CA, BOND_LENGTHS[("CA", "C")], BOND_ANGLES[("N", "CA", "C")], psi)
    new_O = place_atom(new_N, new_CA, new_C,  BOND_LENGTHS[("C", "O")], BOND_ANGLES[("CA","C","O")], 0.0)
    
    res.get_atom("N").set_coor(new_N)
    res.get_atom("CA").set_coor(new_CA)
    res.get_atom("C").set_coor(new_C)
    res.get_atom("O").set_coor(new_O)

    return [new_N, new_CA, new_C]

def build_residue_clashes(res, built, prev, max_backbone=50):
    phi, psi = sample_phi_psi()
    new_prev = build_backbone(res, prev, phi, psi)
    
    if try_build_sidechain(res, built):
        return new_prev, True, "chi"
    
    for attempt in range(max_backbone):
        phi, psi = sample_phi_psi()
        new_prev = build_backbone(res, prev, phi, psi)
        
        if try_build_sidechain(res, built):
            return new_prev, True, "backbone"
        
    return new_prev, False, "failed"    


        
        
        
    
        
    