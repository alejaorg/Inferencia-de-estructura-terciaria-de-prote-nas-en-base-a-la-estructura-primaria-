import sys
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
src = os.path.abspath(os.path.join(BASE_DIR, '..', 'src'))
sys.path.insert(0, src)

# Now you can import the module normally
from protein import *
from energy import *
from energyModel import *
from annealing import *
from proteinParser import export_pdb
from proteinViewer import ProteinViewer

p = Protein("ACDHI")
p.init_random_structure()


#for res in p.residues:
#    print(res)
#    for atom in res.get_atoms():
#        print(atom)
        
        
model = EnergyModel()   
e = Energy(model)
energy = e.getEnergies(p)

print(energy)
        
optimization = Annealing(e, T0=10.0, alpha=0.995, steps=2000)
best_p, best_E, history = optimization.run(p)

print("\n=== RESULTADO FINAL ===")

print("Energia proteina final:", e.getEnergies(best_p))
print("Energía final:", best_E)

PROYECTO_DIR = os.path.abspath(os.path.join(BASE_DIR, '..'))  
DATA_DIR     = os.path.join(PROYECTO_DIR, 'data')            

export_pdb(best_p, "prueba1")
viewer = ProteinViewer(os.path.join(DATA_DIR, "prueba1"))
viewer.show(style="full", color_by="residue")



#E_final, terms_final = e.getEnergies(best_p)
#print("Desglose final:", terms_final)

#clashes_final = detect_clashes(best_p)
#print("Clashes finales:", len(clashes_final))

#print("\n=== ESTRUCTURA FINAL ===")
#for res in best_p.residues:
#    print(res)
#    for atom in res.get_atoms():
#        print(atom)




#for res in p.residues:
#    print(res)
#    for atom in res.get_atoms():
#        print(atom)


