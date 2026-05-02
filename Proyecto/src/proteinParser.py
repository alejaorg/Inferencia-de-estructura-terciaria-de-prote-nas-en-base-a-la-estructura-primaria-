import os
from aminoacid_translation import Atoms_composition 

# ── Define your output folder here ───────────────────────────────────────────
# This file lives in:      project/
#                              your_code/   ← este archivo está aquí
#                              outputs/     ← carpeta hermana donde se guarda el PDB
#
# Adjust the folder name to match yours:
OUTPUT_FOLDER = os.path.join(os.path.dirname(__file__), "..", "data")


def export_pdb(protein, filename):

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    filepath = os.path.join(OUTPUT_FOLDER, filename)

    with open(filepath, "w") as f:
        atom_serial = 1

        for residue in protein.residues:
            res_name = residue.name[:3].upper()
            res_seq  = residue.index
            chain_id = "A"

            atom_lookup   = {atom.name: atom for atom in residue.get_atoms()}
            ordered_names = Atoms_composition.get(res_name, list(atom_lookup.keys()))

            for atom_name in ordered_names:
                atom = atom_lookup.get(atom_name)
                if atom is None:
                    continue

                x, y, z = atom.coord
                element  = atom_name[0]

                if len(atom_name) >= 4:
                    name_field = atom_name[:4].ljust(4)
                elif len(atom_name) == 1:
                    name_field = f" {atom_name}  "
                else:
                    name_field = f" {atom_name:<3}"

                line = (
                    f"{'ATOM':<6}"
                    f"{atom_serial:>5} "
                    f"{name_field:<4}"
                    f" "
                    f"{res_name:<3} "
                    f"{chain_id}"
                    f"{res_seq:>4}    "
                    f"{x:>8.3f}"
                    f"{y:>8.3f}"
                    f"{z:>8.3f}"
                    f"{'1.00':>6}"
                    f"{'0.00':>6}          "
                    f"{element:>2}"
                )

                f.write(line + "\n")
                atom_serial += 1

        f.write("END\n")