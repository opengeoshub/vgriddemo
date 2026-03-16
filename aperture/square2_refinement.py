import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def draw_square_with_rhombus(
    side: float = 1.0,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
    grid_lw: float = 1.5,
    diag_lw: float = 3,
):
    fig, ax = plt.subplots(figsize=(4, 4))
    half = side / 2.0
    A = np.array([-half, -half])  # bottom left
    B = np.array([half, -half])   # bottom right
    C = np.array([half, half])    # top right
    D = np.array([-half, half])   # top left
    outer = np.vstack([A, B, C, D])
    ax.add_patch(
        Polygon(outer, closed=True, facecolor=facecolor, edgecolor=edgecolor, linewidth=grid_lw)
    )
    # Central vertical and horizontal lines (2×2 grid)
    ax.plot([0, 0], [-half, half], color=edgecolor, linewidth=grid_lw)
    ax.plot([-half, half], [0, 0], color=edgecolor, linewidth=grid_lw)
    # Rhombus in each of the 4 small squares
    step = side / 2.0
    for i in range(2):
        for j in range(2):
            x0 = -half + i * step
            y0 = -half + j * step
            cx = x0 + step / 2.0
            cy = y0 + step / 2.0
            mid_left   = (x0,        cy)
            mid_right  = (x0+step,   cy)
            mid_bottom = (cx,        y0)
            mid_top    = (cx,        y0+step)
            rhombus = np.array([mid_left, mid_top, mid_right, mid_bottom])
            ax.add_patch(
                Polygon(rhombus, closed=True, fill=False, edgecolor=edgecolor, linewidth=diag_lw)
            )
    # Centroid of the outer square (red dot)
    # centroid = outer.mean(axis=0)
    # ax.scatter(centroid[0], centroid[1], s=40, c="red", zorder=5)
    ax.set_aspect("equal")
    margin = 0.05 * side
    ax.set_xlim(-half - margin, half + margin)
    ax.set_ylim(-half - margin, half + margin)
    ax.axis("off")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    draw_square_with_rhombus(side=1.0)

