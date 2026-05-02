
class EnergyModel:
    
    def __init__(self):
        
        self.bond_k = 100.0
        self.ca_distance = 3.8
        
        self.local_terms = {
            2: {"d0": 5.5, "k": 5.0},
            3: {"d0": 6.8, "k": 5.0},
            4: {"d0": 6.4, "k": 3.0} 
        }
        
        self.sigma = 6.5
        self.epsilon = 2.0
        
        
        self.hydrophobic_k = 1.0

        self.hydrophobic_scale = {
            # highly hydrophobic
            "VAL": 1.0,
            "ILE": 1.0,
            "LEU": 1.0,
            "MET": 0.9,
            "PHE": 0.9,
            "TRP": 0.9,
            "ALA": 0.7,

            # mean
            "TYR": 0.5,
            "HIS": 0.4,

            # polars
            "SER": 0.0,
            "THR": 0.0,
            "ASN": 0.0,
            "GLN": 0.0,

            # charged (hydrophobic)
            "ASP": -0.5,
            "GLU": -0.5,
            "LYS": -0.5,
            "ARG": -0.5,

            # specials
            "GLY": 0.0,
            "PRO": 0.0,
            "CYS": 0.2
        }
        
        #weights
        self.weights = {
            "bond": 1.0,
            "local": 0.8,
            "non_local": 1.2,
            "hydrophobic": 0.5
        }
