from Bio.PDB import PDBParser
import numpy as np

def parse_pdb_atoms(filepath, atom_name=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", filepath)

    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom_name is None or atom.name == atom_name:
                        atoms.append((
                            residue.get_id()[1],  # índice del residuo
                            atom.name,
                            atom.get_vector().get_array()
                        ))
        break  # solo primer modelo

    return atoms


def align_coords(coords_mobile, coords_ref):
    """
    Aplica el algoritmo de Kabsch para alinear coords_mobile sobre coords_ref.
    Devuelve las coordenadas alineadas.
    """
    # Centrar ambas nubes de puntos
    center_mobile = coords_mobile.mean(axis=0)
    center_ref    = coords_ref.mean(axis=0)

    mobile_c = coords_mobile - center_mobile
    ref_c    = coords_ref    - center_ref

    # Matriz de covarianza
    H = mobile_c.T @ ref_c

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Corrección para evitar reflexiones
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1, 1, d])

    # Matriz de rotación óptima
    R = Vt.T @ D @ U.T

    # Aplicar rotación y traslación
    aligned = (mobile_c @ R.T) + center_ref

    return aligned


def calculate_rmsd(coords1, coords2):
    diff = coords1 - coords2
    return np.sqrt((diff ** 2).sum(axis=1).mean())


def rmsd_to_similarity(rmsd, max_rmsd=20.0):
    similarity = max(0.0, 1.0 - (rmsd / max_rmsd)) * 100
    return round(similarity, 2)


def compare_structures(generated_path, reference_path):
    gen_atoms = parse_pdb_atoms(generated_path)
    ref_atoms = parse_pdb_atoms(reference_path)

    # Emparejar átomos por índice de residuo y nombre de átomo
    gen_dict = {(r, a): c for r, a, c in gen_atoms}
    ref_dict = {(r, a): c for r, a, c in ref_atoms}

    common_keys = sorted(set(gen_dict.keys()) & set(ref_dict.keys()))
    coords_gen = np.array([gen_dict[k] for k in common_keys])
    coords_ref = np.array([ref_dict[k] for k in common_keys])

    print("Alineando estructuras (Kabsch)...")
    coords_aligned = align_coords(coords_gen, coords_ref)

    rmsd = calculate_rmsd(coords_aligned, coords_ref)
    similarity = rmsd_to_similarity(rmsd)

    print(f"\n📊 Resultados:")
    print(f"   Átomos comparados : {len(common_keys)}")
    print(f"   RMSD              : {round(rmsd, 3)} Å")
    print(f"   Similitud         : {similarity}%")

    return rmsd, similarity
