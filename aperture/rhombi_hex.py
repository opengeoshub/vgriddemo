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


def _clip_polygon_to_convex(poly: np.ndarray, clip_poly: np.ndarray) -> np.ndarray:
    """
    Clip polygon `poly` (Nx2) against a convex polygon `clip_poly` (Mx2)
    using the Sutherland–Hodgman algorithm. Both polygons are expected
    to have vertices ordered counter‑clockwise.
    """
    def inside(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> bool:
        # Point p is inside if it is on the left side of edge a->b
        return np.cross(b - a, p - a) >= 0.0

    def intersect(p1: np.ndarray, p2: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        # Intersection of segment p1->p2 with line a->b
        d = p2 - p1
        e = b - a
        denom = np.cross(d, e)
        if abs(denom) < 1e-12:
            return p1  # Parallel; fallback to p1
        t = np.cross(a - p1, e) / denom
        return p1 + t * d

    # Ensure the clipping polygon is counter‑clockwise; if not, reverse it.
    signed_area = 0.0
    for i in range(len(clip_poly)):
        x1, y1 = clip_poly[i]
        x2, y2 = clip_poly[(i + 1) % len(clip_poly)]
        signed_area += (x1 * y2 - x2 * y1)
    if signed_area < 0:
        clip_poly = clip_poly[::-1]

    output = poly
    n_clip = len(clip_poly)
    for i in range(n_clip):
        a = clip_poly[i]
        b = clip_poly[(i + 1) % n_clip]
        input_list = output
        if len(input_list) == 0:
            break
        output = []
        s = input_list[-1]
        for p in input_list:
            if inside(p, a, b):
                if not inside(s, a, b):
                    output.append(intersect(s, p, a, b))
                output.append(p)
            elif inside(s, a, b):
                output.append(intersect(s, p, a, b))
            s = p
        output = np.array(output)
    return output


def draw_rhombi_hexagon(
    n: int = 2,
    side: float = 1.0,
    facecolor: str = "#d8e9f8",
    edgecolor: str = "black",
) -> None:
    """
    Draw a rhombus (made from two equilateral triangles) subdivided
    with ``n = 4`` as in ``draw_rhom_grid`` and, at each vertex where
    four child rhombi meet, place a regular hexagon whose center is
    exactly that vertex.
    """
    fig, ax = plt.subplots(figsize=(4, 7))

    # --- Outer rhombus geometry, matching `draw_rhom_grid` for a given side ---
    # Height of an equilateral triangle with side length "side"
    h = np.sqrt(3) / 2.0 * side

    # Edge vectors of the rhombus in the equilateral-triangle lattice.
    # These have length "side" and meet at 60 degrees.
    u = np.array([side / 2.0, h])   # up-right
    v = np.array([side / 2.0, -h])  # down-right

    # Place the parent rhombus so that its center is at the origin.
    # Parallelogram (rhombus) with base point A and edge vectors u, v:
    # vertices: A, A+u, A+u+v, A+v.
    A = -0.5 * (u + v)
    B = A + u
    C = A + u + v
    D = A + v
    big_rhombus = np.vstack([A, B, C, D])
    # Geometric center of the outer rhombus (used for the hexagon cluster)
    centroid = big_rhombus.mean(axis=0)
    # Draw centroid in red, analogous to rhombi.py
    ax.scatter(centroid[0], centroid[1], s=20, c="red", zorder=5)

    # --- Subdivide into n^2 child rhombi ---
    u_f = u / n
    v_f = v / n

    # Draw the 4×4 children rhombi (optional but helps to debug alignment)
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
                    linewidth=1,
                    linestyle="dashdot",
                    zorder=2,
                )
            )

    # --- Dual vertices: points where 4 children rhombi meet ---
    # Grid vertices in this subdivision are:
    #   G(i,j) = A + i*u_f + j*v_f, for i,j in {0,…,n}.
    # Interior grid points (1..n-1) x (1..n-1) are where 4 rhombi meet.
    centers = []
    for i in range(1, n):
        for j in range(1, n):
            G = A + i * u_f + j * v_f
            centers.append(G)

    centers = np.array(centers)

    # Choose hexagon side relative to the smallest edge length of the child mesh
    edge_lengths = [
        np.linalg.norm(u_f),
        np.linalg.norm(v_f),
        np.linalg.norm(u_f - v_f),
    ]
    base_len = min(edge_lengths)
    # Make the hexagons fit comfortably inside the outer rhombus.
    hex_side_eff = 0.6 * base_len

    # Generate a full pointy‑top hexagonal grid of centers and draw
    # every hexagon that intersects the outer rhombus.
    #
    # Standard axial‑to‑pixel mapping for pointy‑top hexes
    # (see redblobgames.com): for side length s,
    #   x = s * sqrt(3) * (q + r/2)
    #   y = s * 3/2 * r
    def axial_to_xy(q: int, r: int, s: float) -> np.ndarray:
        x = s * np.sqrt(3.0) * (q + 0.5 * r)
        y = s * 1.5 * r
        return np.array([x, y]) + centroid

    # Choose a range large enough so that all hexagons that could
    # intersect the rhombus are included. Scaling with n keeps this
    # roughly proportional to the figure size.
    radius = 2 * n + 2
    all_hexes = []
    for q in range(-radius, radius + 1):
        for r in range(-radius, radius + 1):
            center = axial_to_xy(q, r, hex_side_eff)
            cx, cy = center
            all_hexes.append(regular_hexagon(cx, cy, hex_side_eff))

    # Clip each hexagon to the outer rhombus boundary and draw the result
    for hex_poly in all_hexes:
        clipped = _clip_polygon_to_convex(hex_poly, big_rhombus)
        if clipped.size >= 6:  # at least 3 vertices
            ax.add_patch(
                Polygon(
                    clipped,
                    closed=True,
                    fill=False,
                    edgecolor=edgecolor,
                    linewidth=3,
                    joinstyle="round",
                    zorder=2,
                )
            )

    # Finally draw the outer rhombus, filled
    ax.add_patch(
        Polygon(
            big_rhombus,
            closed=True,
            fill=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=3,
            joinstyle="round",
            zorder=1,
        )
    )

    ax.set_aspect("equal")
    margin_x = side * 0.4
    margin_y = side * 0.4
    xs = big_rhombus[:, 0]
    ys = big_rhombus[:, 1]
    ax.set_xlim(xs.min() - margin_x, xs.max() + margin_x)
    ax.set_ylim(ys.min() - margin_y, ys.max() + margin_y)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Example consistent with rhombi.py: n controls subdivision density
    draw_rhombi_hexagon(n=2, side=1.0)

