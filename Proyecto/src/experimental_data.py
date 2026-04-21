import numpy as np


BOND_LENGTHS = { # covalend bond
    ("N", "CA"): 1.458,
    ("CA", "C"): 1.525,
    ("C", "N"): 1.329,   # peptide bond
    ("C", "O"): 1.229,
    ("CA", "CB"): 1.525
}

BOND_ANGLES = {
    ("N", "CA", "C"): 111.2,
    ("CA", "C", "N"): 116.2,
    ("C", "N", "CA"): 121.7,
    ("CA", "C", "O"): 120.8
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
    "SER": [("N", "CA", "CB", "OG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "VAL": [("N", "CA", "CB", "CG1")],
    "LEU": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1")
    ],
    "ILE": [
        ("N", "CA", "CB", "CG1"),
        ("CA", "CB", "CG1", "CD1")
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