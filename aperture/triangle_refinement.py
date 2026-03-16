import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
import numpy as np


def draw_triangle_grid(n=100, side=1.0, facecolor="#d8e9f8", edgecolor="black"):
    """
    Subdivide an equilateral triangle into n^2 small triangles.

    For an equilateral triangle, subdivision with aperture n^2 (n > 0) is
    achieved by dividing each edge into n equal segments and connecting
    division points so that all small triangle edges are parallel to the
    corresponding edges of the parent triangle.
    """
    fig, ax = plt.subplots(figsize=(4, 4))

    # Height of an equilateral triangle with side length "side"
    h = np.sqrt(3) / 2 * side

    # Parent triangle vertices (centered)
    A = np.array([-side / 2.0, -h / 3.0])  # bottom left
    B = np.array([side / 2.0, -h / 3.0])   # bottom right
    C = np.array([0.0, 2.0 * h / 3.0])     # top

    # Edge subdivision vectors (each edge is divided into n equal segments)
    u = (B - A) / n  # along AB
    v = (C - A) / n  # along AC

    # Draw and fill outer parent triangle once (underneath)
    big_triangle = np.vstack([A, B, C])
    ax.add_patch(
        Polygon(big_triangle, closed=True, facecolor=facecolor, edgecolor=edgecolor, linewidth=2)
    )

    # Draw the center (centroid) of the outer triangle
    centroid = big_triangle.mean(axis=0)
    ax.scatter(centroid[0], centroid[1], s=20, c="red", zorder=10)

    # Generate all small triangles: there are exactly n^2 of them in theory.
    # Draw all child triangles as outlines only on top of the filled parent.

    for i in range(n):
        for j in range(n - i):
            # Base lattice point for this "cell"
            P = A + i * u + j * v

            # Upward-pointing small triangle (edges only)
            up = np.vstack([P, P + u, P + v])
            ax.add_patch(
                Polygon(up, closed=True, fill=False, edgecolor=edgecolor, linewidth=1)
            )

            # Downward-pointing small triangle (exists except on the left edge)
            if i > 0:
                Q = A + i * u + j * v
                down = np.vstack([Q, Q - u, Q + v - u])
                ax.add_patch(
                    Polygon(down, closed=True, fill=False, edgecolor=edgecolor, linewidth=1)
                )

    # Emphasize "direct parent" grids for certain n by overdrawing
    # coarser subdivisions with thicker lines.
    def draw_coarse_grid(n_coarse: int, lw: float):
        u_c = (B - A) / n_coarse
        v_c = (C - A) / n_coarse
        for i in range(n_coarse):
            for j in range(n_coarse - i):
                P_c = A + i * u_c + j * v_c
                up_c = np.vstack([P_c, P_c + u_c, P_c + v_c])
                ax.add_patch(
                    Polygon(up_c, closed=True, fill=False, edgecolor=edgecolor, linewidth=lw)
                )
                if i > 0:
                    Q_c = A + i * u_c + j * v_c
                    down_c = np.vstack([Q_c, Q_c - u_c, Q_c + v_c - u_c])
                    ax.add_patch(
                        Polygon(down_c, closed=True, fill=False, edgecolor=edgecolor, linewidth=lw)
                    )

    # Case n = 4: emphasize parent grid n=2
    if n % 2 == 0:
        draw_coarse_grid(n // 2, lw=3)

    # Case n = 9: emphasize parent grid n=3
    if n % 3 == 0:
        draw_coarse_grid(n // 3, lw=3)

    # Set equal aspect and a tight view around the parent triangle
    ax.set_aspect("equal")
    margin = 0.1 * side
    xs = big_triangle[:, 0]
    ys = big_triangle[:, 1]
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Example: n=4 → aperture = n^2 = 16
    draw_triangle_grid(n=9, side=1.0)

