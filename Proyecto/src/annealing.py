import numpy as np
import copy 
import random
from protein import *

class Annealing:
    def __init__(self, energy_model, T0=10.0, alpha=0.995, steps=10000, stucked_threshold=150):
        
        self.energy = energy_model
        self.T = T0
        self.alpha = alpha
        self.steps = steps
        self.stucked_threshold = stucked_threshold
        
    
    def run(self, protein):
        
        current = protein
        current_E, part_E = self.energy.getEnergies(current)
        
        best = copy.deepcopy(current)
        best_E = current_E
        
        history = []
        
        for step in range(self.steps):
            
            new_protein = copy.deepcopy(current)
            self.perturb(new_protein)
            
            new_E, new_parts = self.energy.getEnergies(new_protein)
            
            dE = new_E - current_E
            
            if self.accept(dE):
                current = new_protein
                current_E = new_E
                
                if new_E < best_E:
                    best = copy.deepcopy(new_protein)
                    best_E = new_E
                    self.stucked_threshold = 150
            else:
                self.stucked_threshold -= 1
            
            if self.stucked_threshold != 0:
                self.cool()
            else:
                self.heat()
                self.stucked_threshold = 150
            history.append(current_E)
            if step % 20 == 0:
                print(f"Step {step:4d} | E={current_E:.3f} | best={best_E:.3f} | T={self.T:.4f}")

        return best, best_E, history
    
    
    def accept(self, dE):
        
        if dE < 0:
            return True
        
        return random.random() < np.exp(-dE / self.T)
    
    def cool(self):
        self.T *= self.alpha
        
    def heat(self):
        self.T *= 3.0
        
    def perturb(self, protein):
        
        i = random.randint(1, len(protein.residues) - 2)
        
        move = random.choice(["phi", "psi", "sidechain"])
        
        if move == "phi":
            self.rotate_phi(protein, i)
        elif move == "psi":
            self.rotate_psi(protein, i)
        else:
            self.rotate_sidechain(protein, i)
            
    def rotate_phi(self, protein, i):
        
        delta = np.radians(random.uniform(-20, 20)) * (self.T / 10)
        
        self.rotate_backbone(protein, i, delta_phi=delta, delta_psi=0)
    
    
    def rotate_psi(self, protein, i):
        
        delta = np.radians(random.uniform(-20, 20)) * (self.T / 10)
        
        self.rotate_backbone(protein, i, delta_phi=0, delta_psi=delta)
            
    
    def rotate_backbone(self, protein, origin, delta_phi=0, delta_psi=0):
        
        prev = [
            protein.residues[origin - 1].get_atom("N").coord,
            protein.residues[origin - 1].get_atom("CA").coord,
            protein.residues[origin - 1].get_atom("C").coord
        ]
        
        built = protein.residues[:origin]
        
        for i in range(origin, len(protein.residues)):
            
            res = protein.residues[i]
            
            phi, psi =  sample_phi_psi()
            
            if i == origin:
                phi += delta_phi
                psi += delta_psi
            
            prev = build_backbone(res, prev, phi, psi)
            
            try_build_sidechain(res, built)
            
            built.append(res)
        
            
    
        
    def rotate_sidechain(self, protein, i):
        
        res = protein.residues[i]
        
        success = try_build_sidechain(res, protein.residues[:i])
        
        if not success:
            for atom in res.get_atoms():
                if atom.name not in ["N", "CA", "C", "O"]:
                    atom.coord += np.random.normal(0, 0.1, 3)   
        
        
        
        
        
        