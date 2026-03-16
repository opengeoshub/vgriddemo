import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

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

def plot_substrate_test(apertures, depth=1, vertex_aligned=False):
    cols = 4
    rows = (len(apertures) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(18, 4.5 * rows))
    axes = axes.flatten()

    for i, N in enumerate(apertures):
        a, b = get_eisenstein_generator(N)
        if a is None: continue

        p_side = 1.0
        parent_poly = get_hex_poly(0, 0, p_side)

        # Scale & Rotation logic
        lin_scale = np.sqrt(N)**depth
        c_side = p_side / lin_scale

        total_rot_rad = np.arctan2(b * np.sqrt(3), 2*a + b) if (depth % 2 != 0) else 0

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

                # Intersection check for clean boundaries
                child_poly = get_hex_poly(cx, cy, c_side, orientation=total_rot_rad)
                if parent_poly.intersects(child_poly):
                    c_x, c_y = child_poly.exterior.xy
                    ax.plot(c_x, c_y, 'b-', lw=0.5, alpha=0.6)
                    ax.scatter(cx, cy, c='blue', s=2, alpha=0.5)

        ax.set_title(f"aperture={N} ({a},{b})", y=-0.1)
        ax.set_aspect('equal')
        ax.axis('off')

    for j in range(i + 1, len(axes)): axes[j].axis('off')
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05)
    plt.show()

# Run for standard apertures in HV (Vertex-Aligned) mode
valid_apertures = [3, 4, 7, 9, 12, 13, 16, 19, 21, 25, 27, 28, 31, 37, 39, 43]
plot_substrate_test(valid_apertures, depth=1, vertex_aligned=False)
