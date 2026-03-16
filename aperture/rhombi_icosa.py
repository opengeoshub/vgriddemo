import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def draw_rhombi_strip(
    count: int = 10,
    side: float = 1.0,
    filled_index: int = 3,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
):
    """
    Draw `count` rhombi built from pairs of equilateral triangles,
    **interleaved and vertically stacked** (two staggered rows)
    instead of lying in a single horizontal row.

    - Rhombi are diamonds with vertical long diagonal, made from two
      equilateral triangles (one up, one down).
    - Even indices sit on a lower row, odd indices on an upper row,
      horizontally interleaved so neighboring rhombi nest.
    - One rhombus (``filled_index``) is shaded.
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    h = np.sqrt(3) / 2.0 * side  # height of an equilateral triangle

    # Base rhombus centered at origin: left, top, right, bottom
    base_vertices = np.array(
        [
            (-side / 2.0, 0.0),
            (0.0, h),
            (side / 2.0, 0.0),
            (0.0, -h),
        ]
    )

    # Horizontal spacing chosen so interleaved rhombi visually touch
    # (no visible gaps between their outlines).
    spacing_x = 0.5 * side
    # Vertical offset between the two staggered rows
    offset_y = h

    # Precompute centers for all rhombi, then recenter the whole layout
    centers = []
    for i in range(count):
        cx = i * spacing_x
        cy = 0.0 if (i % 2 == 0) else offset_y
        centers.append((cx, cy))

    centers = np.array(centers)
    cx_mean = centers[:, 0].mean()
    cy_mean = centers[:, 1].mean()
    centers[:, 0] -= cx_mean
    centers[:, 1] -= cy_mean

    # Draw rhombi
    for i, (cx, cy) in enumerate(centers):
        verts = base_vertices + np.array([cx, cy])

        is_filled = (i == filled_index)
        fc = facecolor if is_filled else "white"

        ax.add_patch(
            Polygon(
                verts,
                closed=True,
                facecolor=fc,
                edgecolor=edgecolor,
                linewidth=3,
            )
        )

        # Thin dash-dot internal horizontal line separating the two triangles
        ax.plot(
            [cx - side / 2.0, cx + side / 2.0],
            [cy, cy],
            color=edgecolor,
            linewidth=0.75,
            linestyle="-.",
        )

    ax.set_aspect("equal")
    margin_x = side
    margin_y = h
    xs = centers[:, 0]
    ys = centers[:, 1]
    total_width = xs.max() - xs.min() + side
    ax.set_xlim(xs.min() - margin_x - side / 2.0, xs.max() + margin_x + side / 2.0)
    ax.set_ylim(ys.min() - margin_y - h, ys.max() + margin_y + h)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Draw 10 rhombi arranged in the ICS-like strip, fill one in the middle.
    draw_rhombi_strip(count=10, side=1.0, filled_index=3)

