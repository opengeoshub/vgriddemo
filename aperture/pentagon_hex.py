import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon


def get_regular_pentagon(center_x: float, center_y: float, side_len: float, rotation: float = 0.0) -> Polygon:
    """
    Regular pentagon on the plane.
    rotation is in radians, applied to the first vertex on +x.
    """
    angles = np.linspace(0, 2 * np.pi, 6)[:-1] + rotation
    points = [
        (center_x + side_len * np.cos(a), center_y + side_len * np.sin(a))
        for a in angles
    ]
    return Polygon(points)


def get_hex_poly(center_x: float, center_y: float, side_len: float, orientation: float = 0.0) -> Polygon:
    """
    Regular hexagon, same helper as in triangle_hex.
    """
    angles = np.linspace(0, 2 * np.pi, 7)[:-1] + orientation
    points = [
        (center_x + side_len * np.cos(a), center_y + side_len * np.sin(a))
        for a in angles
    ]
    return Polygon(points)


def plot_pentagon_hex_grid(side: float = 1.0, aperture: float = 3.0) -> None:
    """
    Build an outer regular pentagon, an inner opposite‑oriented pentagon,
    and a surrounding ring of hexagons, similar to the screenshot example.
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    # Outer pentagon
    p_side = side
    # Rotate so one edge is roughly horizontal
    parent_poly = get_regular_pentagon(0.0, 0.0, p_side, rotation=np.pi / 2.0)

    # Coordinates of outer pentagon (used for limits and clipping only)
    px, py = parent_poly.exterior.xy

    # Inner pentagon: centered, opposite orientation to the outer pentagon.
    # Slightly smaller so it sits well inside.
    cx0, cy0 = parent_poly.centroid.x, parent_poly.centroid.y
    # For similar polygons, area scales with side^2.
    # To make the inner pentagon `aperture` times smaller in area:
    #   inner_side / p_side = 1 / sqrt(aperture)
    inner_scale = 1.0 / np.sqrt(aperture)
    inner_side = inner_scale * p_side
    # Adjust inner pentagon orientation depending on aperture to better
    # match the desired surrounding‑cell layout.
    if aperture == 7.0:
        # Slight rotation relative to outer for aperture=7
        inner_rotation = np.pi / 2.0 + np.pi / 10.0  # 18° offset
    elif aperture == 5.0:
        # Nearly aligned with outer for aperture=5
        inner_rotation = np.pi / 2.0                 # 0° relative offset
    else:
        # Default: flipped relative to outer pentagon
        inner_rotation = np.pi / 2.0 + np.pi / 5.0   # 36° offset
    inner_pent = get_regular_pentagon(cx0, cy0, inner_side, rotation=inner_rotation)
    ipx, ipy = inner_pent.exterior.xy
    ax.plot(ipx, ipy, color="red", lw=2.5, zorder=4)

    # --- Five adjacent non‑regular "hexagons" that tile between
    # the inner pentagon and an outer, scaled copy with no gaps. ---
    inner_coords = np.array(inner_pent.exterior.coords[:-1])  # 5 vertices
    scale_outer = 3.0
    outer_coords = cx0 + (inner_coords - cx0) * scale_outer

    for i in range(5):
        p0 = inner_coords[i]
        p1 = inner_coords[(i + 1) % 5]
        P0 = outer_coords[i]
        P1 = outer_coords[(i + 1) % 5]

        # Construct a 6‑vertex polygon whose inner edge is p0‑p1,
        # outer edge is P0‑P1, and which shares full edges with its
        # neighbors so there are no gaps.
        m0 = 0.5 * (p0 + P0)
        m1 = 0.5 * (p1 + P1)

        hex_like = np.array([p0, p1, m1, P1, P0, m0])
        ring_poly = Polygon(hex_like)

        # Clip the surrounding cell by the original outer pentagon,
        # so that no part extends beyond that boundary.
        inter = parent_poly.intersection(ring_poly)
        if inter.is_empty:
            continue

        geoms = [inter] if isinstance(inter, Polygon) else list(inter.geoms)
        for g in geoms:
            hx, hy = g.exterior.xy
            ax.plot(hx, hy, color="black", lw=1.5, zorder=3)

    ax.set_aspect("equal")
    xs = np.array(px)
    ys = np.array(py)
    margin = side * 0.8
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin)
    ax.axis("off")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_pentagon_hex_grid(side=1.0, aperture=3)

