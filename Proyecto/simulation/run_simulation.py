import sys
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
src = os.path.abspath(os.path.join(BASE_DIR, '..', 'src'))
sys.path.insert(0, src)

from protein import *
from energy import *
from energyModel import *
from annealing import *
from proteinParser import export_pdb

import time

start = time.time()

p = Protein("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
p.init_random_structure()
        
model = EnergyModel()   
e = Energy(model)
energy = e.getEnergies(p)
optimization = Annealing(e, T0=10.0, alpha=0.995, steps=1500)
best_p, best_E, history = optimization.run(p)

print("\n=== FINAL RESULT ===")

PROYECTO_DIR = os.path.abspath(os.path.join(BASE_DIR, '..'))  
DATA_DIR     = os.path.join(PROYECTO_DIR, 'generated_data')            

export_pdb(best_p, "prueba")
print(f"Total time: {time.time() - start:.2f}s")


