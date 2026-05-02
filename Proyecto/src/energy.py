import numpy as np

class Energy:
    
    def __init__(self, model):
        self.model = model
                
    def getEnergies(self, protein):
        
        terms = {
            "bond" : self.bond(protein) * self.model.weights["bond"],
            "local" : self.local(protein) * self.model.weights["local"],
            "nonlocal" : self.non_local(protein) * self.model.weights["non_local"],
            "hydrophobic" : self.hydrophobic(protein) * self.model.weights["hydrophobic"] 
        }
        
        total = sum(terms.values())
        
        return total, terms

    def bond(self, protein):
        
        E = 0.0
        k = self.model.bond_k
        d0 = self.model.ca_distance
        
        for i in range(len(protein.residues) - 1):
            
            ca1 = protein.residues[i].get_atom("CA").coord
            ca2 = protein.residues[i + 1].get_atom("CA").coord
            
            d = np.linalg.norm(ca1 - ca2)
            
            E += k * (d - d0)**2
            
        return E

    def local(self, protein):
        
        E = 0.0
        
        for i in range(len(protein.residues)):
            for offset, params in self.model.local_terms.items():
            
                j = i + offset
                if j >= len(protein.residues):
                    continue
                
                ca1 = protein.residues[i].get_atom("CA").coord
                ca2 = protein.residues[j].get_atom("CA").coord          
                
                d = np.linalg.norm(ca1 - ca2)
                
                E += params["k"] * (d - params["d0"])**2
        return E
    
    def non_local(self, protein):
        
        E = 0.0
        epsilon = self.model.epsilon
        sigma = self.model.sigma
        
        for i in range(len(protein.residues)):
            for j in range(i + 4, len(protein.residues)):
                ca1 = protein.residues[i].get_atom("CA").coord
                ca2 = protein.residues[j].get_atom("CA").coord
                
                r = np.linalg.norm(ca1 - ca2)
                
                if r == 0:
                    continue
                
                sr6 = (sigma / r) ** 6
                sr12 = sr6 ** 2
                
                E += 4 * epsilon * (sr12 - sr6)
        
        return E
    
    def hydrophobic(self, protein):
        
        E = 0.0 
        
        center = np.mean([r.get_atom("CA").coord for r in protein.residues], axis = 0)
        
        for r in protein.residues:
            
            ca = r.get_atom("CA").coord
            
            d = np.linalg.norm(ca - center)
            
            h = self.model.hydrophobic_scale.get(r.name, 0)
            
            E += self.model.hydrophobic_k * h * d**2
            
        return E 
        