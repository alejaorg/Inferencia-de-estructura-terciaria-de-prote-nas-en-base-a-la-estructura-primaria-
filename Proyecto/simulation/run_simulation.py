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
from proteinViewer import ProteinViewer
import time

start = time.time()

p = Protein("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
p.init_random_structure()
        
model = EnergyModel()   
e = Energy(model)
energy = e.getEnergies(p)

print(energy)
        
optimization = Annealing(e, T0=10.0, alpha=0.995, steps=1500)
best_p, best_E, history = optimization.run(p)

print("\n=== FINAL RESULT ===")

PROYECTO_DIR = os.path.abspath(os.path.join(BASE_DIR, '..'))  
DATA_DIR     = os.path.join(PROYECTO_DIR, 'data')            

export_pdb(best_p, "prueba2")
viewer = ProteinViewer(os.path.join(DATA_DIR, "prueba2"))
print(f"Total time: {time.time() - start:.2f}s")
viewer.show(style="full", color_by="residue")




