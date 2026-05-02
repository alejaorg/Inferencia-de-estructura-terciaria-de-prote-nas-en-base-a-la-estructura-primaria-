import sys
import os

src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, os.path.abspath(src))

# Now you can import the module normally
from protein import *
from energy import *
from energyModel import *
from annealing import *

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


