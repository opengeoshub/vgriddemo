import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon


def line_intersection(p1, p2, p3, p4):
    """
    Intersection of infinite lines p1–p2 and p3–p4 in 2D.
    Returns the intersection point as a NumPy array.
    """
    p1, p2, p3, p4 = map(np.asarray, (p1, p2, p3, p4))
    d1 = p2 - p1
    d2 = p4 - p3
    denom = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(denom) < 1e-12:
        return None
    t = ((p3[0] - p1[0]) * d2[1] - (p3[1] - p1[1]) * d2[0]) / denom
    return p1 + t * d1

def get_hex_poly(center_x, center_y, side_len, orientation=0):
    angles = np.linspace(0, 2 * np.pi, 7)[:-1] + orientation
    points = [(center_x + side_len * np.cos(a),
               center_y + side_len * np.sin(a)) for a in angles]
    return Polygon(points)

def get_eisenstein_generator(N):
    for a in range(int(np.sqrt(N)) + 1, 0, -1):
        for b in range(a + 1):
            if a**2 + a*b + b**2 == N: return a, b
    return None, None

def plot_substrate_test(apertures, vertex_aligned=False):
    cols = 3
    rows = (len(apertures) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(18, 4.5 * rows))
    axes = axes.flatten()

    for i, N in enumerate(apertures):
        a, b = get_eisenstein_generator(N)
        if a is None: continue

        p_side = 1.0
        parent_poly = get_hex_poly(0, 0, p_side)

        # Scale & Rotation logic (equivalent to previous depth=1 behavior)
        lin_scale = np.sqrt(N)
        c_side = p_side / lin_scale

        total_rot_rad = np.arctan2(b * np.sqrt(3), 2*a + b)

        # Calculate Vertex Offset (HV Mode)
        off_x, off_y = 0, 0
        if vertex_aligned:
            # Shift the origin to a vertex of a child cell
            off_x = c_side * np.cos(total_rot_rad)
            off_y = c_side * np.sin(total_rot_rad)

        ax = axes[i]
        px, py = parent_poly.exterior.xy
        ax.plot(px, py, 'r-', lw=2, zorder=3)
        # Draw centroid of the outer (parent) hexagon
        cx_parent, cy_parent = parent_poly.centroid.x, parent_poly.centroid.y
        ax.plot(cx_parent, cy_parent, "ro", markersize=5, zorder=4)

        # From the red hexagon, construct an outer triangle by
        # extending three alternating edges and intersecting them.
        # Apply this construction for all apertures so each panel
        # has a corresponding outer triangle.
        verts = np.array(parent_poly.exterior.coords[:-1])  # 6 vertices
        A = line_intersection(verts[0], verts[1], verts[2], verts[3])
        B = line_intersection(verts[2], verts[3], verts[4], verts[5])
        C = line_intersection(verts[4], verts[5], verts[0], verts[1])
        tri = np.vstack([A, B, C, A])
        ax.plot(tri[:, 0], tri[:, 1], color="blue", lw=2.5, zorder=2.5)
        # Mark each triangle vertex with a red dot
        for vx, vy in (A, B, C):
            ax.plot(vx, vy, "ro", markersize=4, zorder=3.5)
        # Use only the outer triangle as the clipping region
        tri_poly = Polygon([A, B, C])
        clip_region = tri_poly

        # Render Children
        limit = int(lin_scale + 2)
        for q in range(-limit, limit + 1):
            for r in range(-limit, limit + 1):
                # Basic axial to cartesian
                tx = c_side * (1.5 * q)
                ty = c_side * (np.sqrt(3)/2 * q + np.sqrt(3) * r)

                # Rotate and apply Vertex Offset
                cx = tx * np.cos(total_rot_rad) - ty * np.sin(total_rot_rad) + off_x
                cy = tx * np.sin(total_rot_rad) + ty * np.cos(total_rot_rad) + off_y

                # Intersection check for clean boundaries against
                # both the parent hexagon and the outer triangle.
                child_poly = get_hex_poly(cx, cy, c_side, orientation=total_rot_rad)
                if clip_region.intersects(child_poly):
                    inter = clip_region.intersection(child_poly)
                    if inter.is_empty:
                        continue
                    # Handle possible MultiPolygon from the intersection
                    geoms = [inter] if isinstance(inter, Polygon) else list(inter.geoms)
                    for g in geoms:
                        c_x, c_y = g.exterior.xy
                        ax.plot(c_x, c_y, color='black', lw=0.5, alpha=0.8)

        # ax.set_title(f"aperture={N} ({a},{b})", y=-0.1)
        ax.set_aspect('equal')
        ax.axis('off')

    for j in range(i + 1, len(axes)): axes[j].axis('off')
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05)
    plt.show()

# Run for standard apertures in HV (Vertex-Aligned) mode
valid_apertures = [3, 4, 7]
plot_substrate_test(valid_apertures, vertex_aligned=False)
