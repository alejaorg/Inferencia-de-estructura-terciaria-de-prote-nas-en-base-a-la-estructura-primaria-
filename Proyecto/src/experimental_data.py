import numpy as np


BOND_LENGTHS = {

# ===== BACKBONE =====
("N", "CA"): 1.458,
("CA", "C"): 1.522,
("C", "N"): 1.329,
("C", "O"): 1.231,

# ===== ALA =====
("CA", "CB"): 1.530,

# ===== ARG =====
("CB", "CG"): 1.520,
("CG", "CD"): 1.520,
("CD", "NE"): 1.470,
("NE", "CZ"): 1.330,
("CZ", "NH1"): 1.320,
("CZ", "NH2"): 1.320,

# ===== ASN =====
# ("CB", "CG"): 1.520,
("CG", "OD1"): 1.250,
("CG", "ND2"): 1.340,

# ===== ASP =====
# ("CB", "CG"): 1.520,
# ("CG", "OD1"): 1.250,
("CG", "OD2"): 1.250,

# ===== CYS =====
("CB", "SG"): 1.810,

# ===== GLN =====
# ("CB", "CG"): 1.520,
# ("CG", "CD"): 1.520,
("CD", "OE1"): 1.250,
("CD", "NE2"): 1.340,

# ===== GLU =====
# ("CB", "CG"): 1.520,
# ("CG", "CD"): 1.520,
# ("CD", "OE1"): 1.250,
("CD", "OE2"): 1.250,

# ===== GLY =====
# (sin cadena lateral)

# ===== HIS =====
# ("CB", "CG"): 1.370,
("CG", "ND1"): 1.320,
("CG", "CD2"): 1.370,
("ND1", "CE1"): 1.350,
("CD2", "NE2"): 1.320,
("CE1", "NE2"): 1.350,

# ===== ILE =====
("CB", "CG1"): 1.520,
("CB", "CG2"): 1.520,
("CG1", "CD1"): 1.520,

# ===== LEU =====
# ("CB", "CG"): 1.520,
("CG", "CD1"): 1.520,
# ("CG", "CD2"): 1.520,

# ===== LYS =====
# ("CB", "CG"): 1.520,
# ("CG", "CD"): 1.520,
("CD", "CE"): 1.520,
("CE", "NZ"): 1.480,

# ===== MET =====
# ("CB", "CG"): 1.520,
("CG", "SD"): 1.780,
("SD", "CE"): 1.780,

# ===== PHE =====
# ("CB", "CG"): 1.510,
# ("CG", "CD1"): 1.390,
# ("CG", "CD2"): 1.390,
("CD1", "CE1"): 1.390,
("CD2", "CE2"): 1.390,
("CE1", "CZ"): 1.390,
("CE2", "CZ"): 1.390,

# ===== PRO =====
# ("CA", "CB"): 1.530,
# ("CB", "CG"): 1.500,
# ("CG", "CD"): 1.500,
("CD", "N"): 1.470,

# ===== SER =====
("CB", "OG"): 1.360,

# ===== THR =====
("CB", "OG1"): 1.360,
# ("CB", "CG2"): 1.520,

# ===== TRP =====
# ("CB", "CG"): 1.500,
# ("CG", "CD1"): 1.365,
# ("CG", "CD2"): 1.430,
("CD1", "NE1"): 1.375,
# ("CD2", "CE2"): 1.400,
("CE2", "CZ2"): 1.390,
("CZ2", "CH2"): 1.390,

# ===== TYR =====
# ("CB", "CG"): 1.510,
# ("CG", "CD1"): 1.390,
# ("CG", "CD2"): 1.390,
# ("CD1", "CE1"): 1.390,
# ("CD2", "CE2"): 1.390,
# ("CE1", "CZ"): 1.390,
("CZ", "OH"): 1.360,

# ===== VAL =====
# ("CB", "CG1"): 1.520,
# ("CB", "CG2"): 1.520,
}

BOND_ANGLES = {

# ===== BACKBONE (común a todos) =====
("N", "CA", "C"): 111.2,
("CA", "C", "N"): 116.2,
("C", "N", "CA"): 121.7,
("CA", "C", "O"): 120.8,

# ===== ALA =====
("N", "CA", "CB"): 110.5,
("CB", "CA", "C"): 110.5,

# ===== ARG =====
# ("N", "CA", "CB"): 110.5,
("CA", "CB", "CG"): 113.8,
("CB", "CG", "CD"): 113.8,
("CG", "CD", "NE"): 113.8,
("CD", "NE", "CZ"): 120.0,
("NE", "CZ", "NH1"): 120.0,
("NE", "CZ", "NH2"): 120.0,

# ===== ASN =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
("CB", "CG", "OD1"): 120.0,
("CB", "CG", "ND2"): 120.0,

# ===== ASP =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
# ("CB", "CG", "OD1"): 120.0,
("CB", "CG", "OD2"): 120.0,

# ===== CYS =====
# ("N", "CA", "CB"): 110.5,
("CA", "CB", "SG"): 113.8,

# ===== GLN =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
# ("CB", "CG", "CD"): 113.8,
("CG", "CD", "OE1"): 120.0,
("CG", "CD", "NE2"): 120.0,

# ===== GLU =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
# ("CB", "CG", "CD"): 113.8,
# ("CG", "CD", "OE1"): 120.0,
("CG", "CD", "OE2"): 120.0,

# ===== GLY =====
#

# ===== HIS =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 122.0,
("CB", "CG", "ND1"): 108.0,
("CB", "CG", "CD2"): 108.0,
("CG", "ND1", "CE1"): 108.0,
("CG", "CD2", "NE2"): 108.0,
("ND1", "CE1", "NE2"): 108.0,

# ===== ILE =====
# ("N", "CA", "CB"): 110.5,
("CA", "CB", "CG1"): 113.8,
("CA", "CB", "CG2"): 113.8,
("CB", "CG1", "CD1"): 113.8,

# ===== LEU =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
# ("CB", "CG", "CD1"): 113.8,
("CB", "CG", "CD2"): 113.8,

# ===== LYS =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
# ("CB", "CG", "CD"): 113.8,
("CG", "CD", "CE"): 113.8,
("CD", "CE", "NZ"): 111.0,

# ===== MET =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
("CB", "CG", "SD"): 100.0,
("CG", "SD", "CE"): 100.0,

# ===== PHE =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
("CB", "CG", "CD1"): 120.0,
# ("CB", "CG", "CD2"): 120.0,
("CG", "CD1", "CE1"): 120.0,
("CG", "CD2", "CE2"): 120.0,
("CD1", "CE1", "CZ"): 120.0,
("CD2", "CE2", "CZ"): 120.0,

# ===== PRO =====
# ("N", "CA", "CB"): 103.0,
# ("CA", "CB", "CG"): 104.5,
# ("CB", "CG", "CD"): 104.5,
("CG", "CD", "N"): 103.0,

# ===== SER =====
# ("N", "CA", "CB"): 110.5,
("CA", "CB", "OG"): 113.0,

# ===== THR =====
# ("N", "CA", "CB"): 110.5,
("CA", "CB", "OG1"): 113.0,
# ("CA", "CB", "CG2"): 113.8,

# ===== TRP =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
# ("CB", "CG", "CD1"): 126.0,
# ("CB", "CG", "CD2"): 126.0,
("CG", "CD1", "NE1"): 108.0,
("CD2", "CE2", "CZ2"): 120.0,

# ===== TYR =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG"): 113.8,
# ("CB", "CG", "CD1"): 120.0,
# ("CB", "CG", "CD2"): 120.0,
# ("CG", "CD1", "CE1"): 120.0,
# ("CG", "CD2", "CE2"): 120.0,
# ("CD1", "CE1", "CZ"): 120.0,
("CE1", "CZ", "OH"): 120.0,

# ===== VAL =====
# ("N", "CA", "CB"): 110.5,
# ("CA", "CB", "CG1"): 113.8,
# ("CA", "CB", "CG2"): 113.8,
}

BACKBONE_TORSIONS = {
    "phi":  ("C_prev", "N", "CA", "C"),
    "psi":  ("N", "CA", "C", "N_next"),
    "omega": ("CA", "C", "N_next", "CA_next")
}

RAMACHANDRAN = {
    "alpha_helix": {
        "phi": (-80, -40),
        "psi": (-70, -20),
        "weight": 0.4
    },
    "beta_sheet": {
        "phi": (-160, -100),
        "psi": (90, 180),
        "weight": 0.4
    },
    "left_handed": {
        "phi": (30, 90),
        "psi": (0, 90),
        "weight": 0.05
    },
    "other": {
        "phi": (-180, 180),
        "psi": (-180, 180),
        "weight": 0.15
    }
}

OMEGA = {
    "trans": {"angle": 180.0, "prob": 0.995},
    "cis": {"angle": 0.0, "prob": 0.005}
}

CA_DISTANCE = {
    "consecutive": 3.8,
    "contact_threshold": 8.0,
    "max": 20.0,
    "bins": np.linspace(2.0, 20.0, 36)
}

CA_DISTANCE_STATS = {
    "mean": 10.0,
    "std": 4.0,
    "min": 3.8,
    "max": 20.0
}

CHI_DEFINITIONS = {

    # Glicina no tiene CB y Alaina solo tiene un CB por lo que no es rotable
    "GLY": [],
    "ALA": [],

    "SER": [("N", "CA", "CB", "OG")],
    "CYS": [("N", "CA", "CB", "SG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "VAL": [("N", "CA", "CB", "CG1")],

    "ASP": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "OD1")
    ],
    "ASN": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "OD1")
    ],
    "HIS": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "ND1")
    ],
    "ILE": [
        ("N", "CA", "CB", "CG1"),
        ("CA", "CB", "CG1", "CD1")
    ],
    "LEU": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1")
    ],
    "PHE": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1")
    ],
    "TYR": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1")
    ],
    "TRP": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1")
    ],

    "GLU": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "OE1")
    ],
    "GLN": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "OE1")
    ],
    "MET": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "SD"),
        ("CB", "CG", "SD", "CE")
    ],

    "LYS": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "CE"),
        ("CG", "CD", "CE", "NZ")
    ],
    "ARG": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "NE"),
        ("CG", "CD", "NE", "CZ")
    ],


    "PRO": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD")
    ]
}

ROTAMERS = {
    "default": [-60, 60, 180]
}

ROTAMER_PROBS = {
    -60: 0.4,
    60: 0.4,
    180: 0.2
}

GEOMETRY = {
    "bond_lengths": BOND_LENGTHS,
    "bond_angles": BOND_ANGLES,
    "torsions": BACKBONE_TORSIONS,
    "ramachandran": RAMACHANDRAN,
    "omega": OMEGA,
    "ca_distance": CA_DISTANCE,
    "ca_distance_stats": CA_DISTANCE_STATS,
    "chi_definitions": CHI_DEFINITIONS,
    "rotamers": ROTAMERS,
    "rotamer_probs": ROTAMER_PROBS
}

VAN_DER_WAALS_RADIUS = {
    "C":  1.70,
    "N":  1.55,
    "O":  1.52,
    "S":  1.80,
    "H":  1.20,
}

VAN_DER_WAALS_TOLERANCE = 0.8