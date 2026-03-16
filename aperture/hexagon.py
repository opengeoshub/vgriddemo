import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def regular_hexagon(center_x: float, center_y: float, side: float) -> np.ndarray:
    """
    Return vertices of a regular hexagon with a vertical long axis
    (pointy top) and given side length.
    """
    angles = np.deg2rad(np.array([90, 150, 210, 270, 330, 30]))
    xs = center_x + side * np.cos(angles)
    ys = center_y + side * np.sin(angles)
    return np.column_stack([xs, ys])


def draw_hexagon_neighbour_distance(
    side: float = 1.0,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
) -> None:
    """
    Draw a small cluster of regular hexagons (one central cell and its
    six neighbours), shade the central hexagon, and draw a single
    distance ray d_h from the centre of the shaded cell to one
    neighbouring hexagon centre. All neighbours are at the same
    distance, so only one ray is shown.
    """
    fig, ax = plt.subplots(figsize=(4, 4))

    # Hexagon side length and derived vertical spacing for pointy‑top layout
    s = side
    dx = np.sqrt(3.0) * s       # horizontal spacing between centres in same row
    dy = 1.5 * s                # vertical spacing between rows

    # Central hexagon at the origin
    cx, cy = 0.0, 0.0

    # Neighbour centre: choose the one to the right (same y)
    cx_n = cx + dx
    cy_n = cy

    # List of hexagons to draw: central + its six neighbours
    centres = [
        (0, 0),          # centre
        (dx, 0),         # right
        (-dx, 0),        # left
        (0.5 * dx, dy),  # upper-right
        (-0.5 * dx, dy), # upper-left
        (0.5 * dx, -dy), # lower-right
        (-0.5 * dx, -dy) # lower-left
    ]

    for (hx, hy) in centres:
        verts = regular_hexagon(hx, hy, s)
        is_central = (hx == 0 and hy == 0)
        ax.add_patch(
            Polygon(
                verts,
                closed=True,
                fill=is_central,
                facecolor=facecolor if is_central else "white",
                edgecolor=edgecolor,
                linewidth=1.5 if is_central else 1.0,
            )
        )

    # Draw central point
    ax.plot(cx, cy, "ro", markersize=5, zorder=3)

    # Draw centroids and rays from the centre to all neighbouring hexagons.
    # Use a solid line for the right neighbour (labelled d_h) and dashed
    # lines for all other neighbours to indicate equal distance.
    for hx, hy in centres:
        if hx == 0 and hy == 0:
            continue
        ax.plot(hx, hy, "ro", markersize=4, zorder=3)
        linestyle = "-" if (hx == cx_n and hy == cy_n) else "--"
        ax.plot([cx, hx], [cy, hy], color=edgecolor, linewidth=1.2, linestyle=linestyle)

    # Label d_h slightly above the solid ray, with the same offset behaviour
    # used in square.py (positioned along the ray and nudged
    # perpendicularly away from it).
    along_offset = 3.0 / 4.0
    px = cx + along_offset * (cx_n - cx)
    py = cy + along_offset * (cy_n - cy)
    dx_ray = cx_n - cx
    dy_ray = cy_n - cy
    length = np.hypot(dx_ray, dy_ray)
    if length > 0:
        nx = -dy_ray / length
        ny = dx_ray / length
        normal_offset = 0.16 * side
        px += normal_offset * nx
        py += normal_offset * ny

    ax.text(px, py, r"$d_h$", fontsize=20, ha="center", va="center")

    ax.set_aspect("equal")
    margin = 0.5 * side
    all_x = [c[0] for c in centres]
    all_y = [c[1] for c in centres]
    ax.set_xlim(min(all_x) - 2 * side, max(all_x) + 2 * side)
    ax.set_ylim(min(all_y) - 2 * side, max(all_y) + 2 * side)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    draw_hexagon_neighbour_distance(side=1.0)

