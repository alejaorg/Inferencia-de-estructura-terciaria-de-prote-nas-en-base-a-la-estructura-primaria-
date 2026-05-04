# Inferencia-de-estructura-terciaria-de-prote-nas-en-base-a-la-estructura-primaria-
Un código realizado en python como práctica final de la asignatura de Bioinformática en la universidad complutense de Madrid con la finalidad de ser una versión simplificada de herramientas como Rosetta.



Dependencias del proyecto: 
    numpy, matplotlib, mpl_toolkits, os, itertools, random, copy, sys

Estructura de carpetas esperada en el proyecto:

    Proyecto/
    ├── simulation/
    │   └── test_pipeline.py , run_simulation.py, ...   
    ├── src/
    │   ├── protein.py, annealing.py, energy.py, ...
    ├── generated_data/         ← PDBs generados por el algoritmo
    └── ground_truth/           ← PDBs descargados del RCSB
        └── proteinas_reales.txt

En caso de querer modificar parametros del modelo se deberá modificar la clase energyModel.py


Uso del pipeline usado para los testeos:

test_pipeline.py
-----------------
1. Lee proteinas_reales.txt
2. Para cada proteína: genera estructura con el annealing
3. Descarga el PDB de referencia del RCSB
4. Compara con el comparator (Kabsch RMSD)
5. Genera tabla de estadísticas en CSV y TXT

Uso:
    python test_pipeline.py

Intrucciones de uso para los run_simlation:

run_simulation.py
-----------------
1. Crea una proteina introducion la cadena de aminoacidos en el constructor de la misma
2. Selecciona el nombre del fichero que quieres que genere poniendolo en el segundo argumento de exportPDB

Uso:
    python run_simulation.py

Intrucciones de uso para los run_viewer:

run_viewer.py
-----------------
1. Selecciona el nombre del fichero .pdb que quieres que sea representado poniendolo en el segundo argumento de ProteinViewer (NOTA: este debe encontrarse en la carpeta generated_data en caso de desear otro origen deberá modificarse en la ruta de DATA_DIR)

Uso:
    python run_viewer.py
    
Intrucciones de uso para los run_comparator:

run_comparator.py
-----------------
1. Se deberá indicar el nombre de los ficheros .pdb a comparar en los argumentos de la funcion compare_structures el primero ha de encontrarse en la carpeta generated_data y el segundo en ground_truth (NOTA: en caso de querer cambiar las carpetas de origen se deberá cambiar las rutas definidas dentro dle fichero como se indica en run_viewer)

Uso:
    python run_comparator.py
        

