import sys
import os


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
src = os.path.abspath(os.path.join(BASE_DIR, '..', 'src'))
sys.path.insert(0, src)


from comparatorPDB import compare_structures


PROYECTO_DIR = os.path.abspath(os.path.join(BASE_DIR, '..'))  
GENERATED_DATA_DIR     = os.path.join(PROYECTO_DIR, 'generated_data')
GROUND_DATA_DIR     = os.path.join(PROYECTO_DIR, 'ground_truth')
     

compare_structures(os.path.join(GENERATED_DATA_DIR, "Crambin"), os.path.join(GROUND_DATA_DIR, "1CRN.pdb"))