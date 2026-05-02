import sys
import os

src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, os.path.abspath(src))

# Now you can import the module normally
from protein import *
from energy import *
from energyModel import *

p = Protein("ACDHI")
p.init_random_structure()


for res in p.residues:
    print(res)
    for atom in res.get_atoms():
        print(atom)
        
        
model = EnergyModel()   
e = Energy(model)
energy = e.getEnergies(p)

print(energy)
        
clashes = detect_clashes(p)
if clashes:
   print("Clash in some point")


#for res in p.residues:
#    print(res)
#    for atom in res.get_atoms():
#        print(atom)


