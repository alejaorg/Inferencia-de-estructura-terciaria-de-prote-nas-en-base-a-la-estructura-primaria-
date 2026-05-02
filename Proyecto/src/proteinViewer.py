"""
protein_viewer.py
-----------------
Clase para visualizar proteínas en 3D a partir de archivos .pdb.

Dependencias:
    pip install matplotlib numpy biopython

Uso básico:
    viewer = ProteinViewer("mi_proteina.pdb")
    viewer.show()

Uso avanzado:
    viewer = ProteinViewer("mi_proteina.pdb")
    viewer.show(style="ribbon", color_by="chain", background="dark")
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from typing import Optional
import warnings

# Paleta de colores para cadenas y tipos de átomos
CHAIN_COLORS = [
    "#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
    "#F4A261", "#9B2335", "#6A4C93", "#1982C4",
]

ATOM_COLORS = {
    "C": "#A0A0A0",
    "N": "#3A86FF",
    "O": "#FF595E",
    "S": "#FFCA3A",
    "H": "#FFFFFF",
    "P": "#FF7700",
    "default": "#CCCCCC",
}

RESIDUE_COLORS = {
    # Hidrofóbicos
    "ALA": "#FF6B6B", "VAL": "#FF8E53", "LEU": "#FFA07A",
    "ILE": "#FFB347", "MET": "#FFD700", "PHE": "#FFC0CB",
    "TRP": "#DA70D6", "PRO": "#EE82EE",
    # Hidrofílicos polares
    "SER": "#87CEEB", "THR": "#ADD8E6", "CYS": "#98FB98",
    "TYR": "#90EE90", "ASN": "#00CED1", "GLN": "#20B2AA",
    # Cargados
    "ASP": "#FF4500", "GLU": "#FF6347",
    "LYS": "#4169E1", "ARG": "#0000CD", "HIS": "#6495ED",
    # Glicina
    "GLY": "#D3D3D3",
    "default": "#AAAAAA",
}


class ProteinViewer:
    """
    Visualizador 3D de proteínas a partir de archivos PDB.

    Parámetros
    ----------
    pdb_path : str
        Ruta al archivo .pdb.

    Ejemplo
    -------
    >>> viewer = ProteinViewer("1crn.pdb")
    >>> viewer.show()
    >>> viewer.show(style="spheres", color_by="element")
    >>> viewer.show(style="ribbon", color_by="residue")
    """

    def __init__(self, pdb_path: str):
        self.pdb_path = pdb_path
        self.atoms = []          # Lista de dicts con info de cada átomo
        self.residues = {}       # Dict: (chain, res_seq) -> list of atoms
        self.chains = {}         # Dict: chain_id -> list of CA atoms
        self._load_pdb()

    # ------------------------------------------------------------------ #
    #  Parseo del archivo PDB (sin dependencias externas)                  #
    # ------------------------------------------------------------------ #

    def _load_pdb(self):
        """Lee el archivo PDB y almacena la información de los átomos."""
        self.atoms.clear()
        self.residues.clear()
        self.chains.clear()

        try:
            with open(self.pdb_path, "r") as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"No se encontró el archivo: {self.pdb_path}")

        for line in lines:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue

            try:
                atom = {
                    "record":   record,
                    "serial":   int(line[6:11]),
                    "name":     line[12:16].strip(),
                    "res_name": line[17:20].strip(),
                    "chain":    line[21].strip() or "A",
                    "res_seq":  int(line[22:26]),
                    "x":        float(line[30:38]),
                    "y":        float(line[38:46]),
                    "z":        float(line[46:54]),
                    "element":  line[76:78].strip() if len(line) > 76 else line[12:14].strip().lstrip("0123456789"),
                }
                self.atoms.append(atom)

                # Agrupar por residuo
                key = (atom["chain"], atom["res_seq"])
                self.residues.setdefault(key, []).append(atom)

                # Guardar alfa-carbonos para el ribbon
                if atom["name"] == "CA":
                    self.chains.setdefault(atom["chain"], []).append(atom)

            except (ValueError, IndexError):
                continue  # Ignorar líneas malformadas

        if not self.atoms:
            raise ValueError("No se encontraron átomos ATOM/HETATM en el archivo PDB.")

        print(f"✓ Cargados {len(self.atoms)} átomos | "
              f"{len(self.residues)} residuos | "
              f"{len(self.chains)} cadena(s): {list(self.chains.keys())}")

    # ------------------------------------------------------------------ #
    #  Helpers de color                                                    #
    # ------------------------------------------------------------------ #

    def _color_by_chain(self, atom: dict) -> str:
        chains_sorted = sorted(self.chains.keys())
        idx = chains_sorted.index(atom["chain"]) if atom["chain"] in chains_sorted else 0
        return CHAIN_COLORS[idx % len(CHAIN_COLORS)]

    def _color_by_element(self, atom: dict) -> str:
        el = atom["element"].upper()[:1]
        return ATOM_COLORS.get(el, ATOM_COLORS["default"])

    def _color_by_residue(self, atom: dict) -> str:
        return RESIDUE_COLORS.get(atom["res_name"], RESIDUE_COLORS["default"])

    def _color_by_bfactor(self, atom: dict) -> str:
        """Placeholder – requeriría leer B-factor; usa degradado azul-rojo."""
        return "#8888FF"

    def _get_color(self, atom: dict, color_by: str) -> str:
        dispatch = {
            "chain":   self._color_by_chain,
            "element": self._color_by_element,
            "residue": self._color_by_residue,
        }
        fn = dispatch.get(color_by, self._color_by_chain)
        return fn(atom)

    # ------------------------------------------------------------------ #
    #  Métodos de renderizado                                              #
    # ------------------------------------------------------------------ #

    def _render_backbone(self, ax, color_by: str, alpha: float = 0.8):
        """Dibuja el backbone conectando alfa-carbonos por cadena."""
        for chain_id, ca_atoms in self.chains.items():
            if len(ca_atoms) < 2:
                continue
            xs = [a["x"] for a in ca_atoms]
            ys = [a["y"] for a in ca_atoms]
            zs = [a["z"] for a in ca_atoms]

            # Líneas coloreadas por segmento
            points = np.array([xs, ys, zs]).T.reshape(-1, 1, 3)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            colors = [self._get_color(ca_atoms[i], color_by) for i in range(len(ca_atoms) - 1)]

            lc = Line3DCollection(segments, colors=colors, linewidths=1.5, alpha=alpha)
            ax.add_collection3d(lc)

    def _render_ribbon(self, ax, color_by: str):
        """
        Pseudo-ribbon: dibuja el backbone con líneas gruesas suavizadas
        (spline interpolado de alfa-carbonos).
        """
        try:
            from scipy.interpolate import splprep, splev
            use_scipy = True
        except ImportError:
            use_scipy = False
            warnings.warn("scipy no disponible; el ribbon usará líneas rectas.")

        for chain_id, ca_atoms in self.chains.items():
            if len(ca_atoms) < 4:
                continue

            xs = np.array([a["x"] for a in ca_atoms])
            ys = np.array([a["y"] for a in ca_atoms])
            zs = np.array([a["z"] for a in ca_atoms])

            if use_scipy:
                tck, u = splprep([xs, ys, zs], s=0, k=min(3, len(xs) - 1))
                u_fine = np.linspace(0, 1, len(ca_atoms) * 10)
                xs, ys, zs = splev(u_fine, tck)

            color = CHAIN_COLORS[sorted(self.chains.keys()).index(chain_id) % len(CHAIN_COLORS)]
            ax.plot(xs, ys, zs, color=color, linewidth=3, alpha=0.85)

    def _render_spheres(self, ax, color_by: str, atom_filter: Optional[str] = None,
                        max_atoms: int = 2000):
        """Dibuja átomos como esferas (scatter 3D)."""
        atoms = self.atoms
        if atom_filter:
            atoms = [a for a in atoms if a["name"] == atom_filter]

        # Submuestrear si hay demasiados átomos
        if len(atoms) > max_atoms:
            step = len(atoms) // max_atoms
            atoms = atoms[::step]

        xs = [a["x"] for a in atoms]
        ys = [a["y"] for a in atoms]
        zs = [a["z"] for a in atoms]
        colors = [self._get_color(a, color_by) for a in atoms]

        ax.scatter(xs, ys, zs, c=colors, s=20, alpha=0.7, edgecolors="none", depthshade=True)

    def _render_backbone_and_sidechains(self, ax, color_by: str):
        """
        Estilo similar a AlphaFold:
        - Backbone grueso conectando N -> CA -> C de cada residuo
        - Sidechains como líneas finas saliendo del CA
        """
        # Átomos del backbone principal
        BACKBONE_ATOMS = {"N", "CA", "C", "O"}

        # ── Backbone ──────────────────────────────────────────────────
        for chain_id, ca_atoms in self.chains.items():
            # Reconstruir orden N->CA->C por residuo
            backbone_coords = []
            backbone_colors = []

            residue_keys = sorted(
                [k for k in self.residues if k[0] == chain_id],
                key=lambda k: k[1]
            )

            for key in residue_keys:
                res_atoms = self.residues[key]
                atom_map = {a["name"]: a for a in res_atoms}
                for bname in ("N", "CA", "C"):
                    a = atom_map.get(bname)
                    if a:
                        backbone_coords.append([a["x"], a["y"], a["z"]])
                        backbone_colors.append(self._get_color(a, color_by))

            if len(backbone_coords) < 2:
                continue

            pts = np.array(backbone_coords).reshape(-1, 1, 3)
            segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
            lc = Line3DCollection(segs, colors=backbone_colors[:-1], linewidths=2.5, alpha=0.9)
            ax.add_collection3d(lc)

        # ── Sidechains ────────────────────────────────────────────────
        sidechain_segs = []
        sidechain_colors = []

        for key, res_atoms in self.residues.items():
            atom_map = {a["name"]: a for a in res_atoms}
            ca = atom_map.get("CA")
            if ca is None:
                continue

            for a in res_atoms:
                if a["name"] in BACKBONE_ATOMS:
                    continue  # solo sidechain
                sidechain_segs.append([
                    [ca["x"], ca["y"], ca["z"]],
                    [a["x"],  a["y"],  a["z"]],
                ])
                sidechain_colors.append(self._get_color(a, color_by))

        if sidechain_segs:
            lc2 = Line3DCollection(
                sidechain_segs, colors=sidechain_colors,
                linewidths=0.9, alpha=0.6
            )
            ax.add_collection3d(lc2)

        # ── Esferas en CA ─────────────────────────────────────────────
        for chain_id, ca_atoms in self.chains.items():
            xs = [a["x"] for a in ca_atoms]
            ys = [a["y"] for a in ca_atoms]
            zs = [a["z"] for a in ca_atoms]
            colors = [self._get_color(a, color_by) for a in ca_atoms]
            ax.scatter(xs, ys, zs, c=colors, s=30, alpha=0.95,
                       edgecolors="none", depthshade=True, zorder=5)

    def _render_sticks(self, ax, color_by: str, max_bonds: int = 5000):
        """Dibuja enlaces entre átomos del mismo residuo."""
        bond_segments = []
        bond_colors = []
        count = 0

        for (chain, res_seq), res_atoms in self.residues.items():
            if count >= max_bonds:
                break
            for i, a1 in enumerate(res_atoms):
                for a2 in res_atoms[i + 1:]:
                    dist = np.sqrt((a1["x"] - a2["x"])**2 +
                                   (a1["y"] - a2["y"])**2 +
                                   (a1["z"] - a2["z"])**2)
                    if dist < 2.0:  # Umbral de enlace covalente (Å)
                        bond_segments.append([
                            [a1["x"], a1["y"], a1["z"]],
                            [a2["x"], a2["y"], a2["z"]],
                        ])
                        bond_colors.append(self._get_color(a1, color_by))
                        count += 1

        if bond_segments:
            lc = Line3DCollection(bond_segments, colors=bond_colors, linewidths=0.8, alpha=0.6)
            ax.add_collection3d(lc)

    # ------------------------------------------------------------------ #
    #  API pública                                                         #
    # ------------------------------------------------------------------ #

    def show(
        self,
        style: str = "backbone",
        color_by: str = "chain",
        background: str = "dark",
        title: Optional[str] = None,
        figsize: tuple = (10, 8),
        save_path: Optional[str] = None,
    ):
        """
        Muestra la proteína en 3D.

        Parámetros
        ----------
        style : str
            Estilo de renderizado:
            - "backbone"  : línea por alfa-carbonos (rápido)
            - "ribbon"    : cinta suavizada (requiere scipy)
            - "spheres"   : esfera por átomo
            - "sticks"    : enlaces covalentes
            - "ca_spheres": solo alfa-carbonos como esferas
        color_by : str
            Esquema de colores: "chain" | "element" | "residue"
        background : str
            "dark" o "light"
        title : str, opcional
            Título del gráfico. Por defecto usa el nombre del archivo.
        figsize : tuple
            Tamaño de la figura en pulgadas.
        save_path : str, opcional
            Ruta donde guardar la imagen (p.ej. "output.png").
        """

        valid_styles = ("backbone", "ribbon", "spheres", "sticks", "ca_spheres", "full")
        if style not in valid_styles:
            raise ValueError(f"style debe ser uno de: {valid_styles}")

        bg_color = "#0D0D0D" if background == "dark" else "#F5F5F5"
        text_color = "#EEEEEE" if background == "dark" else "#111111"

        fig = plt.figure(figsize=figsize, facecolor=bg_color)
        ax = fig.add_subplot(111, projection="3d", facecolor=bg_color)

        # --- Renderizar según estilo ---
        if style == "backbone":
            self._render_backbone(ax, color_by)
        elif style == "ribbon":
            self._render_ribbon(ax, color_by)
        elif style == "spheres":
            self._render_spheres(ax, color_by)
        elif style == "ca_spheres":
            self._render_spheres(ax, color_by, atom_filter="CA", max_atoms=99999)
        elif style == "sticks":
            self._render_sticks(ax, color_by)
        elif style == "full":
            self._render_backbone_and_sidechains(ax, color_by)

        # --- Estética de los ejes ---
        ax.set_xlabel("X (Å)", color=text_color, labelpad=5)
        ax.set_ylabel("Y (Å)", color=text_color, labelpad=5)
        ax.set_zlabel("Z (Å)", color=text_color, labelpad=5)
        ax.tick_params(colors=text_color)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor("#333333" if background == "dark" else "#CCCCCC")
        ax.yaxis.pane.set_edgecolor("#333333" if background == "dark" else "#CCCCCC")
        ax.zaxis.pane.set_edgecolor("#333333" if background == "dark" else "#CCCCCC")
        ax.grid(True, color="#222222" if background == "dark" else "#DDDDDD", alpha=0.4)

        # --- Título ---
        if title is None:
            import os
            title = os.path.basename(self.pdb_path)
        fig.suptitle(
            title,
            color=text_color,
            fontsize=14,
            fontweight="bold",
            y=0.97,
        )

        # --- Leyenda de cadenas ---
        if color_by == "chain":
            chains_sorted = sorted(self.chains.keys())
            handles = [
                plt.Line2D([0], [0], color=CHAIN_COLORS[i % len(CHAIN_COLORS)],
                           linewidth=3, label=f"Cadena {c}")
                for i, c in enumerate(chains_sorted)
            ]
            ax.legend(handles=handles, loc="upper left",
                      facecolor=bg_color, labelcolor=text_color, framealpha=0.5)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight", facecolor=bg_color)
            print(f"✓ Imagen guardada en: {save_path}")

        plt.show()

    def info(self):
        """Imprime un resumen de la proteína cargada."""
        print(f"\n{'='*50}")
        print(f"  Archivo  : {self.pdb_path}")
        print(f"  Átomos   : {len(self.atoms)}")
        print(f"  Residuos : {len(self.residues)}")
        print(f"  Cadenas  : {list(self.chains.keys())}")
        xs = [a["x"] for a in self.atoms]
        ys = [a["y"] for a in self.atoms]
        zs = [a["z"] for a in self.atoms]
        print(f"  Rango X  : [{min(xs):.1f}, {max(xs):.1f}] Å")
        print(f"  Rango Y  : [{min(ys):.1f}, {max(ys):.1f}] Å")
        print(f"  Rango Z  : [{min(zs):.1f}, {max(zs):.1f}] Å")
        print(f"{'='*50}\n")


# ------------------------------------------------------------------ #
#  Ejecución directa                                                   #
# ------------------------------------------------------------------ #

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Uso: python protein_viewer.py <archivo.pdb> [style] [color_by]")
        print("  style   : backbone | ribbon | spheres | sticks | ca_spheres  (default: backbone)")
        print("  color_by: chain | element | residue                          (default: chain)")
        sys.exit(1)

    pdb_file = sys.argv[1]
    style    = sys.argv[2] if len(sys.argv) > 2 else "backbone"
    color_by = sys.argv[3] if len(sys.argv) > 3 else "chain"

    viewer = ProteinViewer(pdb_file)
    viewer.info()
    viewer.show(style=style, color_by=color_by)