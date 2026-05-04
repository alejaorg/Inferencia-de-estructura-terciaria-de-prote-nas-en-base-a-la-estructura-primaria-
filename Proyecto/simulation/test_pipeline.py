

import sys
import os
import urllib.request
import time
import csv

# ── Rutas ─────────────────────────────────────────────────────────────────────
BASE_DIR     = os.path.dirname(os.path.abspath(__file__))
PROYECTO_DIR = os.path.abspath(os.path.join(BASE_DIR, '..'))
SRC_DIR      = os.path.join(PROYECTO_DIR, 'src')
GEN_DIR      = os.path.join(PROYECTO_DIR, 'generated_data')
GT_DIR       = os.path.join(PROYECTO_DIR, 'ground_truth')
PROTEINS_TXT = os.path.join(GT_DIR, 'proteinas_reales.txt')

sys.path.insert(0, SRC_DIR)

os.makedirs(GEN_DIR, exist_ok=True)
os.makedirs(GT_DIR,  exist_ok=True)

# ── Imports del proyecto ──────────────────────────────────────────────────────
from protein      import Protein
from energy       import Energy
from energyModel  import EnergyModel
from annealing    import Annealing
from proteinParser import export_pdb
from comparatorPDB   import compare_structures

# Configuración del annealing
T0    = 10.0
ALPHA = 0.995
STEPS = 1000

# Leer proteínas del txt
def load_proteins(path):
    proteins = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            proteins.append({
                "name": parts[0],
                "pdb":  parts[1],
                "aa":   int(parts[2]),
                "seq":  parts[3].replace(" ", ""),
                "desc": parts[4] if len(parts) > 4 else "",
            })
    return proteins

# Descargar PDB de referencia
def download_pdb(pdb_code, dest_path):
    if os.path.exists(dest_path):
        print(f"  [GT] {pdb_code} ya descargado.")
        return True
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    try:
        urllib.request.urlretrieve(url, dest_path)
        print(f"  [GT] {pdb_code} descargado.")
        return True
    except Exception as e:
        print(f"  [GT] Error descargando {pdb_code}: {e}")
        return False

# Generar estructura 
def generate_structure(seq, out_path_noext):
    model        = EnergyModel()
    e            = Energy(model)
    optimization = Annealing(e, T0=T0, alpha=ALPHA, steps=STEPS)

    p = Protein(seq)
    p.init_random_structure()

    best_p, best_E, _ = optimization.run(p)

    export_pdb(best_p, out_path_noext)
    return best_E

# Pipeline principal
def run_pipeline():
    proteins = load_proteins(PROTEINS_TXT)
    print(f"Proteínas cargadas: {len(proteins)}\n")

    results = []

    for i, prot in enumerate(proteins):
        name     = prot["name"]
        pdb_code = prot["pdb"]
        seq      = prot["seq"]
        aa       = prot["aa"]

        print(f"{'='*60}")
        print(f"[{i+1}/{len(proteins)}] {name} ({pdb_code}, {aa} AA)")
        print(f"{'='*60}")

        # Nombre de archivo seguro 
        safe_name = name.replace(" ", "_").replace("/", "-").replace("(", "").replace(")", "")

        gen_path_noext = os.path.join(GEN_DIR, safe_name)
        gt_path        = os.path.join(GT_DIR,  f"{pdb_code}.pdb")

        # 1. Generar estructura 
        print(f"  [GEN] Generando estructura...")
        t0 = time.time()
        try:
            best_E = generate_structure(seq, gen_path_noext)
            gen_time = round(time.time() - t0, 1)
            print(f"  [GEN] Energia final: {best_E:.3f} | Tiempo: {gen_time}s")
        except Exception as ex:
            print(f"  [GEN] ERROR: {ex}")
            results.append({
                "name": name, "pdb": pdb_code, "aa": aa,
                "status": "ERROR_GEN", "rmsd": None,
                "similarity": None, "energy": None, "time_s": None
            })
            continue

        # Detectar extensión que usó export_pdb
        if os.path.exists(gen_path_noext + ".pdb"):
            gen_path = gen_path_noext + ".pdb"
        elif os.path.exists(gen_path_noext):
            gen_path = gen_path_noext
        else:
            print(f"  [GEN] No se encontró el archivo generado.")
            results.append({
                "name": name, "pdb": pdb_code, "aa": aa,
                "status": "ERROR_FILE", "rmsd": None,
                "similarity": None, "energy": best_E, "time_s": gen_time
            })
            continue

        # 2. Descargar referencia
        ok = download_pdb(pdb_code, gt_path)
        if not ok:
            results.append({
                "name": name, "pdb": pdb_code, "aa": aa,
                "status": "ERROR_DL", "rmsd": None,
                "similarity": None, "energy": best_E, "time_s": gen_time
            })
            continue

        # 3. Comparar 
        print(f"  [CMP] Comparando estructuras...")
        try:
            rmsd, similarity = compare_structures(gen_path, gt_path)
            print(f"  [CMP] RMSD: {rmsd:.3f} Å | Similitud: {similarity}%")
            results.append({
                "name": name, "pdb": pdb_code, "aa": aa,
                "status": "OK", "rmsd": round(rmsd, 3),
                "similarity": similarity, "energy": round(best_E, 3),
                "time_s": gen_time
            })
        except Exception as ex:
            print(f"  [CMP] ERROR: {ex}")
            results.append({
                "name": name, "pdb": pdb_code, "aa": aa,
                "status": "ERROR_CMP", "rmsd": None,
                "similarity": None, "energy": round(best_E, 3),
                "time_s": gen_time
            })

        print()

    # 4. Guardar resultados 
    save_results(results)

# Guardar tabla de resultados
def save_results(results):
    # CSV
    csv_path = os.path.join(PROYECTO_DIR, 'resultados_batch.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["name","pdb","aa","status","rmsd","similarity","energy","time_s"])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nCSV guardado en: {csv_path}")

    # TXT legible
    txt_path = os.path.join(PROYECTO_DIR, 'resultados_batch.txt')
    ok_results = [r for r in results if r["status"] == "OK"]

    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("=" * 90 + "\n")
        f.write("  RESULTADOS BATCH — INFERENCIA DE ESTRUCTURA TERCIARIA\n")
        f.write("=" * 90 + "\n\n")

        # Tabla
        header = f"{'Proteína':<35} {'PDB':<6} {'AA':<5} {'RMSD(Å)':<10} {'Similitud%':<12} {'Energía':<14} {'Tiempo(s)':<10} {'Estado'}"
        f.write(header + "\n")
        f.write("-" * 90 + "\n")

        for r in results:
            rmsd_s  = f"{r['rmsd']:.3f}"      if r['rmsd']       is not None else "—"
            sim_s   = f"{r['similarity']:.2f}" if r['similarity'] is not None else "—"
            e_s     = f"{r['energy']:.1f}"     if r['energy']     is not None else "—"
            time_s  = f"{r['time_s']}"         if r['time_s']     is not None else "—"
            line = f"{r['name']:<35} {r['pdb']:<6} {r['aa']:<5} {rmsd_s:<10} {sim_s:<12} {e_s:<14} {time_s:<10} {r['status']}"
            f.write(line + "\n")

        f.write("-" * 90 + "\n\n")

        # Estadísticas globales
        if ok_results:
            rmsds = [r['rmsd'] for r in ok_results]
            sims  = [r['similarity'] for r in ok_results]
            times = [r['time_s'] for r in ok_results]

            f.write("ESTADÍSTICAS GLOBALES\n")
            f.write(f"  Proteínas procesadas : {len(results)}\n")
            f.write(f"  Completadas con éxito: {len(ok_results)}\n")
            f.write(f"  Errores              : {len(results) - len(ok_results)}\n\n")
            f.write(f"  RMSD medio           : {sum(rmsds)/len(rmsds):.3f} Å\n")
            f.write(f"  RMSD mínimo          : {min(rmsds):.3f} Å  ({ok_results[rmsds.index(min(rmsds))]['name']})\n")
            f.write(f"  RMSD máximo          : {max(rmsds):.3f} Å  ({ok_results[rmsds.index(max(rmsds))]['name']})\n\n")
            f.write(f"  Similitud media      : {sum(sims)/len(sims):.2f}%\n")
            f.write(f"  Mejor similitud      : {max(sims):.2f}%  ({ok_results[sims.index(max(sims))]['name']})\n")
            f.write(f"  Peor similitud       : {min(sims):.2f}%  ({ok_results[sims.index(min(sims))]['name']})\n\n")
            f.write(f"  Tiempo medio/proteína: {sum(times)/len(times):.1f}s\n")
            f.write(f"  Tiempo total         : {sum(times):.1f}s\n")

    print(f"TXT guardado en:  {txt_path}")

# Entry point
if __name__ == "__main__":
    run_pipeline()