import matplotlib.pyplot as plt
import numpy as np


def draw_square_neighbour_distances(
    side: float = 1.0,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
) -> None:
    """
    Draw a 3x3 grid of squares with the center cell shaded and two
    distance rays d_{s1}, d_{s2} from the center of the shaded cell
    to the centers of its right and upper-right neighbours, matching
    the reference figure.
    """
    fig, ax = plt.subplots(figsize=(4, 4))

    # Outer 3x3 grid geometry
    n = 3
    half = side * 1.5  # total width = 3 * cell_size = 3 * side
    cell = side

    # Coordinates for grid lines
    xs = np.linspace(-half, half, n + 1)
    ys = np.linspace(-half, half, n + 1)

    # Draw vertical grid lines
    for x in xs:
        ax.plot([x, x], [ys[0], ys[-1]], color=edgecolor, linewidth=3)

    # Draw horizontal grid lines
    for y in ys:
        ax.plot([xs[0], xs[-1]], [y, y], color=edgecolor, linewidth=3)

    # Shade the center cell
    # Index (1,1) in the 3x3 grid, using (i,j) row/column indices from 0..2
    i_c = 1
    j_c = 1
    x0 = xs[j_c]
    x1 = xs[j_c + 1]
    y0 = ys[i_c]
    y1 = ys[i_c + 1]
    ax.fill_between([x0, x1], y0, y1, color=facecolor, zorder=0)

    # Centroid of the shaded (central) cell
    cx = (x0 + x1) / 2.0
    cy = (y0 + y1) / 2.0
    # Red dot at the center of the central square
    ax.plot(cx, cy, "ro", markersize=5, zorder=3)

    # Centers of the right neighbour and upper-right neighbour
    # Right neighbour (same row, one column to the right)
    x0_r = xs[j_c + 1]
    x1_r = xs[j_c + 2]
    cx_r = (x0_r + x1_r) / 2.0
    cy_r = cy

    # Upper-right neighbour (one row up, one column to the right)
    y0_ur = ys[i_c + 1]
    y1_ur = ys[i_c + 2]
    cx_ur = cx_r
    cy_ur = (y0_ur + y1_ur) / 2.0

    # Draw the two distance rays
    def draw_ray(
        target_x: float,
        target_y: float,
        label: str,
        along_offset: float = 0.6,
        normal_offset: float = 0.1,
    ):
        ax.plot([cx, target_x], [cy, target_y], color=edgecolor, linewidth=1.5)
        # Red dot at the end of the ray
        ax.plot(target_x, target_y, "ro", markersize=4)

        # Position along the ray
        px = cx + along_offset * (target_x - cx)
        py = cy + along_offset * (target_y - cy)

        # Small offset perpendicular to the ray so the label sits "above" it
        dx = target_x - cx
        dy = target_y - cy
        length = (dx**2 + dy**2) ** 0.5
        if length > 0:
            nx = -dy / length  # rotate (dx,dy) by +90 degrees and normalize
            ny = dx / length
            px += normal_offset * nx
            py += normal_offset * ny

        ax.text(px, py, label, fontsize=20, ha="center", va="center")

    # Slightly different along-offsets so labels don't overlap rays or dots
    draw_ray(cx_r, cy_r, r"$d_{s1}$", along_offset=3.0 / 4.0)
    # Place d_{s2} at about 2/3 along its ray
    draw_ray(cx_ur, cy_ur, r"$d_{s2}$", along_offset=3.0 / 4.0)

    ax.set_aspect("equal")
    margin = 0.1 * side
    ax.set_xlim(xs[0] - margin, xs[-1] + margin)
    ax.set_ylim(ys[0] - margin, ys[-1] + margin)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    draw_square_neighbour_distances(side=1.0)

