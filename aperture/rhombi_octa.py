import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def draw_rhombi(
    count: int = 5,
    side: float = 1.0,
    filled_index: int = 1,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
):
    """
    Draw a horizontal row of identical rhombuses built from
    two equilateral triangles, in vertical (point-up / point-down)
    orientation, with one of them filled.

    - ``count``: how many rhombuses in the row.
    - ``side``: edge length of the equilateral triangles (also the
      edge length of each rhombus).
    - ``filled_index``: zero-based index of the rhombus to fill; if out of
      range, none are filled.
    """
    fig, ax = plt.subplots(figsize=(6, 2))

    # Height of an equilateral triangle with side length "side"
    h = np.sqrt(3) / 2.0 * side

    # Rhombus made of two equilateral triangles, centered at origin,
    # with long diagonal vertical:
    # left, top, right, bottom
    base_vertices = np.array(
        [
            (-side / 2.0, 0.0),
            (0.0, h),
            (side / 2.0, 0.0),
            (0.0, -h),
        ]
    )

    # Horizontal spacing so rhombi touch at the left/right vertices
    spacing = side

    for i in range(count):
        cx = (i - (count - 1) / 2.0) * spacing
        cy = 0.0
        verts = base_vertices + np.array([cx, cy])

        is_filled = (i == filled_index)
        ax.add_patch(
            Polygon(
                verts,
                closed=True,
                facecolor=facecolor if is_filled else "white",
                edgecolor=edgecolor,
                linewidth=3,
            )
        )

        # Thin dash-dot line separating the two equilateral triangles that form
        # each rhombus (shared horizontal edge through the center).
        ax.plot(
            [cx - side / 2.0, cx + side / 2.0],
            [cy, cy],
            color=edgecolor,
            linewidth=0.75,
            linestyle="-.",
        )

    ax.set_aspect("equal")
    margin = side * 0.5
    total_width = (count - 1) * spacing + side
    ax.set_xlim(-total_width / 2.0 - margin, total_width / 2.0 + margin)
    ax.set_ylim(-h - margin, h + margin)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Draw 4 vertical rhombi in a row, fill the second one (index 1)
    draw_rhombi(count=4, side=1.0, filled_index=1)

