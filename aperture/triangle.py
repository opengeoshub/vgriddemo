import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def draw_triangle_with_neighbour_distances(
    n: int = 4,
    side: float = 1.0,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
) -> None:
    """
    Recreate the refined equilateral triangle subdivision from
    ``triangle_refinement.py`` with refinement n (default 4), then
    choose the triangle whose centroid is closest to the centroid of
    the big triangle and draw three distances from it to three
    different "types" of neighbours:

    - d1: neighbour sharing an entire edge,
    - d2: neighbour sharing only a vertex (no common edge),
    - d3: neighbour not touching the central triangle at all.
    """
    fig, ax = plt.subplots(figsize=(4, 4))

    # Height of an equilateral triangle with side length "side"
    h = np.sqrt(3) / 2 * side

    # Parent triangle vertices (centered), same as triangle_refinement.py
    A = np.array([-side / 2.0, -h / 3.0])  # bottom left
    B = np.array([side / 2.0, -h / 3.0])   # bottom right
    C = np.array([0.0, 2.0 * h / 3.0])     # top

    # Edge subdivision vectors
    u = (B - A) / n  # along AB
    v = (C - A) / n  # along AC

    # Outer parent triangle (used only to set view limits; keep unfilled)
    big_triangle = np.vstack([A, B, C])
    ax.add_patch(
        Polygon(
            big_triangle,
            closed=True,
            fill=False,
            facecolor="none",
            edgecolor=edgecolor,
            linewidth=2,
        )
    )

    # Centroid of the outer triangle
    big_centroid = big_triangle.mean(axis=0)

    # Collect all small triangles (both orientations) along with their centroids
    small_tris = []  # list of (vertices, centroid)

    for i in range(n):
        for j in range(n - i):
            P = A + i * u + j * v

            # Upward-pointing small triangle
            up = np.vstack([P, P + u, P + v])
            small_tris.append((up, up.mean(axis=0)))

            # Downward-pointing small triangle (exists except on the left edge)
            if i > 0:
                Q = A + i * u + j * v
                down = np.vstack([Q, Q - u, Q + v - u])
                small_tris.append((down, down.mean(axis=0)))

    # Choose the "central" triangle: centroid closest to big triangle centroid
    centroids = np.array([c for (_, c) in small_tris])
    d_to_big = np.linalg.norm(centroids - big_centroid, axis=1)
    central_idx = int(np.argmin(d_to_big))
    central_tri, central_centroid = small_tris[central_idx]

    # Draw all small triangles, filling only the central one with facecolor
    for idx, (tri, centroid) in enumerate(small_tris):
        is_central = (idx == central_idx)
        ax.add_patch(
            Polygon(
                tri,
                closed=True,
                fill=is_central,
                facecolor=facecolor if is_central else "white",
                edgecolor=edgecolor,
                linewidth=1.5 if is_central else 1.0,
            )
        )

    ax.set_aspect("equal")
    # Use all triangle vertices and centroids to set a view box that
    # guarantees every centroid is visible.
    all_points = [big_triangle] + [tri for (tri, _) in small_tris] + [centroids]
    all_points = np.vstack(all_points)
    margin = 0.1 * side
    xs = all_points[:, 0]
    ys = all_points[:, 1]
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    draw_triangle_with_neighbour_distances(n=3, side=1.0)

