import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def draw_square_grid(n: int = 4, side: float = 1.0, facecolor: str = "#d8e9f8", edgecolor: str = "black"):
    """
    Subdivide a square into n^2 smaller squares, analogous to the triangle
    subdivision in `triangle.py`:

    - The outer parent square is filled once with `facecolor`.
    - All child cells (n x n grid of small squares) are drawn as outlines
      on top of the filled parent, so the whole interior has the same color.
    - When n is a multiple of 2 or 3, a coarser "parent" grid (n/2 or n/3)
      is overdrawn with thicker lines to emphasize the direct parents.
    """
    fig, ax = plt.subplots(figsize=(4, 4))

    # Parent square vertices, centered at origin for symmetry
    half = side / 2.0
    A = np.array([-half, -half])  # bottom left
    B = np.array([half, -half])   # bottom right
    C = np.array([half, half])    # top right
    D = np.array([-half, half])   # top left

    # Draw and fill outer parent square once (background)
    big_square = np.vstack([A, B, C, D])
    ax.add_patch(
        Polygon(big_square, closed=True, facecolor=facecolor, edgecolor=edgecolor, linewidth=2)
    )

    # Draw the center (centroid) of the outer square
    centroid = big_square.mean(axis=0)
    ax.scatter(centroid[0], centroid[1], s=20, c="red", zorder=10)


    # Edge subdivision vectors for the fine grid (n segments per side)
    u = (B - A) / n  # along bottom edge (x direction)
    v = (D - A) / n  # along left edge (y direction)

    # Draw all n^2 child squares as outlines only
    for i in range(n):
        for j in range(n):
            P = A + i * u + j * v
            sq = np.vstack([P, P + u, P + u + v, P + v])
            ax.add_patch(
                Polygon(sq, closed=True, fill=False, edgecolor=edgecolor, linewidth=1)
            )

    # Helper to draw coarser "parent" grids with thicker lines
    def draw_coarse_grid(n_coarse: int, lw: float):
        step_u = (B - A) / n_coarse
        step_v = (D - A) / n_coarse
        for i in range(n_coarse):
            for j in range(n_coarse):
                P_c = A + i * step_u + j * step_v
                sq_c = np.vstack([P_c, P_c + step_u, P_c + step_u + step_v, P_c + step_v])
                ax.add_patch(
                    Polygon(sq_c, closed=True, fill=False, edgecolor=edgecolor, linewidth=lw)
                )

    # Emphasize parent grids analogous to triangle.py behaviour:
    # - if n is divisible by 2, draw grid with n/2
    # - if n is divisible by 3, also draw grid with n/3
    if n % 2 == 0:
        draw_coarse_grid(n // 2, lw=3)
    if n % 3 == 0:
        draw_coarse_grid(n // 3, lw=3)

    # Draw the center point of the outer square
    centroid = big_square.mean(axis=0)
    ax.scatter(centroid[0], centroid[1], s=20, c=edgecolor, zorder=5)

    # Camera / canvas settings similar to triangle.py
    ax.set_aspect("equal")
    margin = 0.1 * side
    xs = big_square[:, 0]
    ys = big_square[:, 1]
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Example: n=4 → aperture = n^2 = 16
    draw_square_grid(n=4, side=1.0)

