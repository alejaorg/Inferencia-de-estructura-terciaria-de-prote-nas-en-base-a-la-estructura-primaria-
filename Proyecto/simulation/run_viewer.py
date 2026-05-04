import sys
import os


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
src = os.path.abspath(os.path.join(BASE_DIR, '..', 'src'))
sys.path.insert(0, src)


from proteinViewer import ProteinViewer


PROYECTO_DIR = os.path.abspath(os.path.join(BASE_DIR, '..'))  
DATA_DIR     = os.path.join(PROYECTO_DIR, 'generated_data')     

viewer = ProteinViewer(os.path.join(DATA_DIR, "Crambin"))
viewer.show(style="full", color_by="residue")