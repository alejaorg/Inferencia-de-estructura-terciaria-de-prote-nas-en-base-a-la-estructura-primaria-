import numpy as np
import copy 
import random

class Annealing:
    def __init__(self, energy_model, T0=10.0, alpha=0.995, steps=10000):
        
        self.energy = energy_model
        self.T = T0
        self.alpha = alpha
        self.steps = steps
        
    
    def run(self):
        