import numpy as np
from aminoacid_translation import ( Atoms_composition, Aminoacid_map )
from experimental_data import (
    BOND_LENGTHS,
    BOND_ANGLES,
    RAMACHANDRAN,
    OMEGA,
    VAN_DER_WAALS_RADIUS,
    VAN_DER_WAALS_TOLERANCE,
    BIFURCATIONS,
    ARG_BIFURCATION
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
        
        if residue.name == "TYR":
            OH = place_atom(
                CE1, CE2, CZ,
                BOND_LENGTHS.get(("CZ", "OH"), 1.36),
                BOND_ANGLES.get(("CE1", "CZ", "OH"), 120.0),
                180.0
            )
            residue.get_atom("OH").set_coor(OH)
            atoms["OH"] = OH

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
            
        if residue.name in BIFURCATIONS:
            place_bifurcation(residue, atoms)
            
        if residue.name == "TRP":
            place_TRP(residue, atoms, chi_angles[1])
            
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
        "NH2": "CZ",
        "OH": "CZ",
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

    
    for i in range(len(atoms)):
        res_i, atom_i = atoms[i]
        for j in range(i + 1, len(atoms)):
            res_j, atom_j = atoms[j]

            
            if res_i.index == res_j.index:
                continue

            backbone = {"N", "CA", "C", "O"}
            if abs(res_i.index - res_j.index) == 1:
                if atom_i.name in backbone and atom_j.name in backbone:
                    continue

            dist = np.linalg.norm(atom_i.coord - atom_j.coord)

            r_i = VAN_DER_WAALS_RADIUS.get(atom_i.name[0], 1.70)
            r_j = VAN_DER_WAALS_RADIUS.get(atom_j.name[0], 1.70)
            threshold = (r_i + r_j) * tolerance
            
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
    
    backbone = {"N", "CA", "C", "O"}

    new_atoms  = [a for a in res.get_atoms() if not np.all(a.coord == 0.0)]
    if not new_atoms:
        return []

    coords_new  = np.array([a.coord for a in new_atoms])      
    names_new   = [a.name for a in new_atoms]
    # Radio VdW de cada átomo nuevo
    radii_new   = np.array([VAN_DER_WAALS_RADIUS.get(a.name[0], 1.70)
                            for a in new_atoms])                 

    old_atoms   = []
    old_indices = []   
    for r in other_res:
        for a in r.get_atoms():
            if not np.all(a.coord == 0.0):
                old_atoms.append(a)
                old_indices.append(r.index)

    if not old_atoms:
        return []

    coords_old  = np.array([a.coord for a in old_atoms])        
    names_old   = [a.name for a in old_atoms]
    radii_old   = np.array([VAN_DER_WAALS_RADIUS.get(a.name[0], 1.70)
                            for a in old_atoms])                 
    old_indices = np.array(old_indices)                          

    diff  = coords_new[:, np.newaxis, :] - coords_old[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=2))                  

    thresholds = (radii_new[:, np.newaxis] + radii_old[np.newaxis, :]) * tolerance

    clash_mask = dists < thresholds                             

    res_index = res.index
    consecutive = np.abs(res_index - old_indices) == 1         

    for mi, name_i in enumerate(names_new):
        if name_i in backbone:
            for ki, name_k in enumerate(names_old):
                if consecutive[ki] and name_k in backbone:
                    clash_mask[mi, ki] = False

    if not np.any(clash_mask):
        return []

    clashes = []
    mi_arr, ki_arr = np.where(clash_mask)
    for mi, ki in zip(mi_arr, ki_arr):
        clashes.append((names_new[mi], names_old[ki], float(dists[mi, ki])))

    return clashes

def place_bifurcation(residue, atoms):
    
    if residue.name == "ARG":
        b = ARG_BIFURCATION
        
        pre = grand_grandparent(atoms, b["grandparent"])
        if pre is None:
            return
        
        for atom, offset in [(b["atom1"], b["offset1"]),
                             (b["atom2"], b["offset2"])]:
            length = BOND_LENGTHS.get(("CZ", atom)) 
            b_angle = BOND_ANGLES.get(("NE", "CZ", atom))
            
            coord = place_atom(pre, 
                               atoms[b["grandparent"]],
                               atoms[b["parent"]],
                               length,
                               b_angle,
                               offset)  
            residue.get_atom(atom).set_coor(coord)
            atoms[atom] = coord
        return 
    
    atom1, atom2, parent, grandparent, offset = BIFURCATIONS[residue.name]
    
    pre = grand_grandparent(atoms, grandparent)
    if pre is None:
        return
    
    chi = calculate_dihedral(pre,
                            atoms[grandparent],
                            atoms[parent],
                            atoms[atom1])
    
    length = BOND_LENGTHS.get((parent, atom2))
    b_angle = BOND_ANGLES.get((grandparent, parent, atom2))
    
    coord = place_atom(pre,
                       atoms[grandparent],
                       atoms[parent],
                       length,
                       b_angle,
                       chi + offset)
    residue.get_atom(atom2).set_coor(coord)
    atoms[atom2] = coord
    
def grand_grandparent(atoms, grandparent):
    chain = ["N", "CA", "CB", "CG", "CD", "NE", "CZ"]
    if grandparent in chain:
        idx = chain.index(grandparent)
        if idx > 0:
            grandGrandpa = chain[idx-1]
            return atoms.get(grandGrandpa, None)
    return None
    
def calculate_dihedral(a, b, c , d):
    b1 = b - a
    b2 = c - b
    b3 = d - c

    n1 = normalize(np.cross(b1, b2))
    n2 = normalize(np.cross(b2, b3))
    m = np.cross(n1, normalize(b2))
    
    x = np.dot(n1, n2)
    y = np.dot(m, n2)
    
    return np.degrees(np.arctan2(y, x))

def place_TRP(residue, atoms, chi2):
    
    CA = atoms["CA"]
    CB = atoms["CB"]
    CG = atoms["CG"]
    CD1 = atoms["CD1"]
    
    CD2 = place_atom(CA, CB, CG,
                     BOND_LENGTHS[("CG", "CD2")],
                     BOND_ANGLES[("CB", "CG", "CD2")],
                     chi2 + 180.0)
    residue.get_atom("CD2").set_coor(CD2)
    atoms["CD2"] = CD2
    
    NE1 = place_atom(CB, CG, CD1,
                     BOND_LENGTHS[("CD1", "NE1")],
                     BOND_ANGLES[("CG", "CD1", "NE1")],
                     180.0)
    residue.get_atom("NE1").set_coor(NE1)
    atoms["NE1"] = NE1

    CE2 = place_atom(CG, CD1, NE1,
                     BOND_LENGTHS[("NE1", "CE2")],
                     BOND_ANGLES[("CD1", "NE1", "CE2")],
                     180.0)
    residue.get_atom("CE2").set_coor(CE2)
    atoms["CE2"] = CE2
    
    CE3 = place_atom(CB, CG, CD2,
                     BOND_LENGTHS[("CD2", "CE3")],
                     BOND_ANGLES[("CG", "CD2", "CE3")],
                     180.0)
    residue.get_atom("CE3").set_coor(CE3)
    atoms["CE3"] = CE3
    
    CZ2 = place_atom(CD1, NE1, CE2,
                     BOND_LENGTHS[("CE2", "CZ2")],
                     BOND_ANGLES[("NE1", "CE2", "CZ2")],
                     180.0)
    residue.get_atom("CZ2").set_coor(CZ2)
    atoms["CZ2"] = CZ2

    CZ3 = place_atom(CG, CD2, CE3,
                     BOND_LENGTHS[("CE3", "CZ3")],
                     BOND_ANGLES[("CD2", "CE3", "CZ3")],
                     180.0)
    residue.get_atom("CZ3").set_coor(CZ3)
    atoms["CZ3"] = CZ3
    
    CH2 = place_atom(NE1, CE2, CZ2,
                     BOND_LENGTHS[("CZ2", "CH2")],
                     BOND_ANGLES[("CE2", "CZ2", "CH2")],
                     180.0)
    residue.get_atom("CH2").set_coor(CH2)
    atoms["CH2"] = CH2

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


        
        
        
    
        
    