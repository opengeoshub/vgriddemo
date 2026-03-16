import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def draw_rhom_grid(
    n: int = 4,
    side: float = 1.0,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
):
    """
    Subdivide a rhombus made from two regular (equilateral) triangles
    into n^2 smaller, congruent rhombi, mirroring the logic of
    `draw_square_grid` in `square.py` but using a triangular lattice
    like `rhombi_octa.py`.

    - The outer parent rhombus is built from two equilateral triangles
      sharing a horizontal edge.
    - Each child cell is also a rhombus that can be decomposed into two
      equilateral triangles of side length ``side / n``.
    """
    fig, ax = plt.subplots(figsize=(4, 4))

    # Height of an equilateral triangle with side length "side"
    h = np.sqrt(3) / 2.0 * side

    # Edge vectors of the rhombus in the equilateral-triangle lattice.
    # These have length "side" and meet at 60 degrees.
    u = np.array([side / 2.0, h])   # up-right
    v = np.array([side / 2.0, -h])  # down-right

    # Place the parent rhombus so that its center is at the origin.
    # Parallelogram (rhombus) with base point A and edge vectors u, v:
    # vertices: A, A+u, A+u+v, A+v
    # Choose A so that centroid is at (0,0).
    A = -0.5 * (u + v)
    B = A + u
    C = A + u + v
    D = A + v

    big_rhombus = np.vstack([A, B, C, D])

    # Draw and fill outer parent rhombus once (background)
    ax.add_patch(
        Polygon(big_rhombus, closed=True, facecolor=facecolor, edgecolor=edgecolor, linewidth=2)
    )

    # Draw the center (centroid) of the outer rhombus
    centroid = big_rhombus.mean(axis=0)
    ax.scatter(centroid[0], centroid[1], s=20, c="red", zorder=10)

    # Fine-grid edge vectors: divide u and v into n segments
    u_f = u / n
    v_f = v / n

    # Draw all n^2 child rhombi as outlines only (small cells with dash-dot lines)
    for i in range(n):
        for j in range(n):
            P = A + i * u_f + j * v_f
            cell = np.vstack([P, P + u_f, P + u_f + v_f, P + v_f])
            ax.add_patch(
                Polygon(
                    cell,
                    closed=True,
                    fill=False,
                    edgecolor=edgecolor,
                    linewidth=1
                )
            )

    # Helper to draw coarser "parent" grids with thicker lines
    def draw_coarse_grid(n_coarse: int, lw: float):
        u_c = u / n_coarse
        v_c = v / n_coarse
        for i in range(n_coarse):
            for j in range(n_coarse):
                P_c = A + i * u_c + j * v_c
                cell_c = np.vstack([P_c, P_c + u_c, P_c + u_c + v_c, P_c + v_c])
                ax.add_patch(
                    Polygon(cell_c, closed=True, fill=False, edgecolor=edgecolor, linewidth=lw)
                )

    # Emphasize parent grids analogous to triangle.py and square.py behaviour:
    if n % 2 == 0:
        draw_coarse_grid(n // 2, lw=3)
    if n % 3 == 0:
        draw_coarse_grid(n // 3, lw=3)

    # Draw the center point again in edgecolor for visibility
    centroid = big_rhombus.mean(axis=0)
    ax.scatter(centroid[0], centroid[1], s=20, c=edgecolor, zorder=5)

    # Camera / canvas settings similar to square.py
    ax.set_aspect("equal")
    margin = 0.1 * side
    xs = big_rhombus[:, 0]
    ys = big_rhombus[:, 1]
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Example: n=9 → aperture = n^2 = 81
    draw_rhom_grid(n=9, side=1.0)

