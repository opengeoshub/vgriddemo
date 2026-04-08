import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString


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


def get_hex_poly(center_x, center_y, side_len, orientation=0.0):
    angles = np.linspace(0, 2 * np.pi, 7)[:-1] + orientation
    points = [
        (center_x + side_len * np.cos(a), center_y + side_len * np.sin(a))
        for a in angles
    ]
    return Polygon(points)


def _point_to_segment_dist(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """Signed distance from point p to line through a, b (segment not extended)."""
    p, a, b = np.asarray(p), np.asarray(a), np.asarray(b)
    ab = b - a
    ap = p - a
    denom = np.dot(ab, ab) + 1e-20
    t = np.clip(np.dot(ap, ab) / denom, 0.0, 1.0)
    proj = a + t * ab
    return np.linalg.norm(p - proj)


def _segment_on_ref(edge_a: np.ndarray, edge_b: np.ndarray,
                    ref_a: np.ndarray, ref_b: np.ndarray,
                    dist_tol: float = 0.08, angle_tol: float = 0.2) -> bool:
    """True if hex edge (edge_a, edge_b) lies on reference segment (ref_a, ref_b)."""
    edge_a, edge_b = np.asarray(edge_a), np.asarray(edge_b)
    ref_a, ref_b = np.asarray(ref_a), np.asarray(ref_b)
    mid = 0.5 * (edge_a + edge_b)
    if _point_to_segment_dist(mid, ref_a, ref_b) > dist_tol:
        return False
    d_edge = edge_b - edge_a
    d_ref = ref_b - ref_a
    le = np.linalg.norm(d_edge)
    lr = np.linalg.norm(d_ref)
    if le < 1e-12 or lr < 1e-12:
        return False
    cos_a = np.dot(d_edge, d_ref) / (le * lr)
    return abs(cos_a) >= 1.0 - angle_tol  # parallel or anti-parallel


def get_eisenstein_generator(N: int):
    """
    Find integers (a, b) such that a^2 + a b + b^2 = N, if possible.
    For N = 4 this returns (2, 0), which gives the standard
    class I (axial) orientation used here.
    """
    for a in range(int(np.sqrt(N)) + 1, 0, -1):
        for b in range(a + 1):
            if a * a + a * b + b * b == N:
                return a, b
    return None, None


def plot_triangle_hex_classI(
    aperture: int = 4,
    step: int = 1,
    draw_hex: bool = False,
    draw_extras: bool = False,
):
    """
    Generate the class I tessellation of an aperture-4 hexagon
    inside a regular triangle, similar to the reference image.
    """
    a, b = get_eisenstein_generator(aperture)
    if a is None:
        raise ValueError(f"No Eisenstein generator found for aperture={aperture}")

    # ------------------------------------------------------------------
    # Step 1: Construct the outer equilateral triangle and its altitudes
    # ------------------------------------------------------------------
    side = 1.0
    h = np.sqrt(3.0) / 2.0 * side

    # Same coordinates as used in triangle.py
    A = np.array([-side / 2.0, -h / 3.0])  # bottom left
    B = np.array([side / 2.0, -h / 3.0])   # bottom right
    C = np.array([0.0, 2.0 * h / 3.0])     # top

    # Triangle polygon used as clipping region
    tri = np.vstack([A, B, C, A])
    tri_poly = Polygon([A, B, C])

    # Centroid (coincides with intersection of altitudes/medians)
    centroid = (A + B + C) / 3.0

    # Midpoints of the three sides
    mid_AB = 0.5 * (A + B)
    mid_BC = 0.5 * (B + C)
    mid_CA = 0.5 * (C + A)

    # Linear scale factor and child cell side length
    lin_scale = np.sqrt(aperture)
    child_side = (side / 2.0) / lin_scale  # fit neatly within the triangle

    # Orientation for class I (derived from the Eisenstein generator)
    total_rot_rad = np.arctan2(b * np.sqrt(3.0), 2.0 * a + b)

    fig, ax = plt.subplots(figsize=(5, 5))

    # Draw the outer triangle
    ax.plot(tri[:, 0], tri[:, 1], color="black", lw=2.0, zorder=2.5)

    # First‑step subdivision: three segments from centroid to side midpoints
    # (thin; thick "tree" is drawn separately when step >= 2).
    mids = (mid_AB, mid_BC, mid_CA)
    ray_dirs = []
    base_lw = 4.0
    thin_lw = 0.8
    base_len = np.linalg.norm(mid_AB - centroid)
    for M in mids:
        ax.plot(
            [centroid[0], M[0]],
            [centroid[1], M[1]],
            color="black",
            lw=thin_lw,
            zorder=3,
        )
        v = M - centroid
        v_len = np.linalg.norm(v)
        ray_dirs.append(v / v_len if v_len > 1e-12 else np.zeros_like(v))

    # Second‑step: extend each line towards its corresponding vertex.
    extended_points = []
    if step >= 2:
        verts = (C, A, B)  # vertices corresponding to mid_AB, mid_BC, mid_CA
        for M, V in zip(mids, verts):
            v = V - centroid
            u = M - centroid
            u_len = np.linalg.norm(u)
            v_len = np.linalg.norm(v)
            if v_len < 1e-12:
                extended_points.append(M)
                continue
            t = u_len / v_len
            P = centroid + t * v  # new point on the ray towards the vertex
            extended_points.append(P)
            ax.plot(
                [centroid[0], P[0]],
                [centroid[1], P[1]],
                color="black",
                lw=thin_lw,
                zorder=3,
            )

        # Thick "tree": when not drawing hex, use geometric segments; when
        # drawing hex, skeleton is drawn as thick hex edges below.
        if len(extended_points) >= 3 and not draw_hex:
            ax.plot(
                [C[0], centroid[0]], [C[1], centroid[1]],
                color="black", lw=base_lw, zorder=3.5,
            )
            P_left = extended_points[1]
            P_right = extended_points[2]
            ax.plot(
                [centroid[0], P_left[0]], [centroid[1], P_left[1]],
                color="black", lw=base_lw, zorder=3.5,
            )
            ax.plot(
                [centroid[0], P_right[0]], [centroid[1], P_right[1]],
                color="black", lw=base_lw, zorder=3.5,
            )
            seg_len = 0.6 * base_len
            down = np.array([0.0, -1.0])
            for P_branch in (P_left, P_right):
                Q_branch = P_branch + seg_len * down
                seg = LineString([tuple(P_branch), tuple(Q_branch)])
                inter = seg.intersection(tri_poly)
                if inter.is_empty:
                    continue
                geoms = [inter] if isinstance(inter, LineString) else list(inter.geoms)
                for g in geoms:
                    coords = np.array(g.coords)
                    if len(coords) >= 2:
                        ax.plot(
                            coords[:, 0], coords[:, 1],
                            color="black", lw=base_lw, zorder=3.5,
                        )

    # Extra edges (Step 2): for each kite, from the endpoint of its ray
    # after extension, draw two short segments parallel to the *inner*
    # kite edges that emanate from the centroid (not the outer triangle
    # sides).
    if draw_extras and step >= 2 and extended_points:
        # Ray directions from the centroid (order corresponds to mids):
        #  - ray_dirs[0]: vertical (centroid -> mid_AB, bottom)
        #  - ray_dirs[1]: up-right  (centroid -> mid_BC)
        #  - ray_dirs[2]: up-left   (centroid -> mid_CA)
        #
        # For each kite we want the new segments parallel to the two
        # centroid‑based kite edges:
        #  - top kite (vertical ray):  two slanted edges → use rays 1 and 2
        #  - lower-right kite:        vertical + up-right → rays 0 and 1
        #  - lower-left kite:         vertical + up-left  → rays 0 and 2

        def n(v):
            L = np.linalg.norm(v)
            return v / L if L > 1e-12 else v

        r0 = n(ray_dirs[0])
        r1 = n(ray_dirs[1])
        r2 = n(ray_dirs[2])

        # For each kite, use the directions of the two *other* rays from
        # the centre (so they are parallel to centre‑based edges but do
        # not coincide with the ray that defines that kite).
        kite_dirs = [
            (r1, r2),  # top kite: use the two slanted rays
            (r0, r2),  # lower-right kite: vertical + opposite slanted ray
            (r0, r1),  # lower-left kite:  vertical + opposite slanted ray
        ]

        # Length of added segments relative to centroid‑to‑midpoint
        base_len = np.linalg.norm(mid_AB - centroid)
        seg_len = 0.6 * base_len

        for P, (d1, d2) in zip(extended_points, kite_dirs):
            for d in (d1, d2):
                d_unit = n(d)
                Q = P + seg_len * d_unit
                seg = LineString([tuple(P), tuple(Q)])
                inter = seg.intersection(tri_poly)
                if inter.is_empty:
                    continue
                # inter can be a LineString or MultiLineString; handle both
                if isinstance(inter, LineString):
                    lines = [inter]
                else:
                    lines = list(inter.geoms)
                for ln in lines:
                    (x0, y0), (x1, y1) = ln.coords[0], ln.coords[-1]
                    ax.plot(
                        [x0, x1],
                        [y0, y1],
                        color="black",
                        lw=thin_lw,
                        zorder=4,
                    )

    # Step 3: repeat the same construction inside three "small kites"
    # nearer the triangle corners. We approximate these by moving further
    # along each main ray and repeating the parallel-edge pattern with
    # shorter segments.
    if draw_extras and step >= 3 and extended_points:
        base_len = np.linalg.norm(mid_AB - centroid)
        if base_len > 1e-12:
            # reuse kite_dirs but with scaled segments. The new extra
            # edges will have 3/4 the length of the old extra edges from
            # step 2, and for each of the three small kites we again draw
            # two parallel edges.
            seg_len = 0.6 * base_len               # length used in step 2
            small_seg_len = (3.0 / 4.0) * seg_len  # 3/4 of the old extra length

            kite_dirs = [
                (r1, r2),  # top small kite
                (r0, r2),  # lower-right small kite
                (r0, r1),  # lower-left small kite
            ]
            for i, (P, (d1, d2)) in enumerate(zip(extended_points, kite_dirs)):
                main_vec = P - centroid
                main_len = np.linalg.norm(main_vec)
                if main_len < 1e-12:
                    continue
                main_dir = main_vec / main_len

                # First, add the "extra" central edge continuing the same ray
                # from P to a new point Q_main further toward the corner.
                Q_main = P + small_seg_len * main_dir
                seg_main = LineString([tuple(P), tuple(Q_main)])
                inter_main = seg_main.intersection(tri_poly)
                if not inter_main.is_empty:
                    if isinstance(inter_main, LineString):
                        lines_main = [inter_main]
                    else:
                        lines_main = list(inter_main.geoms)
                    for ln in lines_main:
                        (x0, y0), (x1, y1) = ln.coords[0], ln.coords[-1]
                        ax.plot(
                            [x0, x1],
                            [y0, y1],
                            color="black",
                            lw=thin_lw,
                            zorder=4,
                        )

                # Then add the two parallel edges inside each small kite,
                # starting at Q_main so they form a small Y‑shape, just like
                # step 2 but scaled down.
                for d in (d1, d2):
                    d_unit = n(d)
                    Q = Q_main + small_seg_len * d_unit
                    seg = LineString([tuple(Q_main), tuple(Q)])
                    inter = seg.intersection(tri_poly)
                    if inter.is_empty:
                        continue
                    if isinstance(inter, LineString):
                        lines = [inter]
                    else:
                        lines = list(inter.geoms)
                    for ln in lines:
                        (x0, y0), (x1, y1) = ln.coords[0], ln.coords[-1]
                        ax.plot(
                            [x0, x1],
                            [y0, y1],
                            color="black",
                            lw=thin_lw,
                            zorder=4,
                        )

    # Optional: render class‑I hex tessellation; skeleton edges (tree) thick, rest thin
    if draw_hex:
        clip_region = tri_poly
        limit = int(lin_scale + 2)

        # Build skeleton reference segments when step >= 2 (apex→center→branches→drops)
        skeleton_refs = []
        if step >= 2 and len(extended_points) >= 3:
            P_left = extended_points[1]
            P_right = extended_points[2]
            seg_len = 0.6 * base_len
            down = np.array([0.0, -1.0])
            skeleton_refs = [
                (np.asarray(C), np.asarray(centroid)),
                (np.asarray(centroid), np.asarray(P_left)),
                (np.asarray(centroid), np.asarray(P_right)),
                (np.asarray(P_left), np.asarray(P_left) + seg_len * down),
                (np.asarray(P_right), np.asarray(P_right) + seg_len * down),
            ]

        # Collect unique hex boundary segments (clipped to triangle)
        def _norm_edge(a, b):
            a = tuple(round(float(x), 10) for x in np.asarray(a).flat)
            b = tuple(round(float(x), 10) for x in np.asarray(b).flat)
            return (a, b) if a <= b else (b, a)

        seen = set()
        all_edges = []
        for q in range(-limit, limit + 1):
            for r in range(-limit, limit + 1):
                tx = child_side * (1.5 * q)
                ty = child_side * (np.sqrt(3.0) / 2.0 * q + np.sqrt(3.0) * r)
                cx = tx * np.cos(total_rot_rad) - ty * np.sin(total_rot_rad)
                cy = tx * np.sin(total_rot_rad) + ty * np.cos(total_rot_rad)
                child_poly = get_hex_poly(cx, cy, child_side, orientation=total_rot_rad)
                if not clip_region.intersects(child_poly):
                    continue
                inter = clip_region.intersection(child_poly)
                if inter.is_empty:
                    continue
                geoms = [inter] if isinstance(inter, Polygon) else list(inter.geoms)
                for g in geoms:
                    coords = list(g.exterior.coords)[:-1]
                    for i in range(len(coords)):
                        a, b = np.array(coords[i]), np.array(coords[(i + 1) % len(coords)])
                        key = _norm_edge(a, b)
                        if key not in seen:
                            seen.add(key)
                            all_edges.append((a, b))

        # Classify and draw: thin first, then thick skeleton edges
        for (a, b) in all_edges:
            on_skeleton = False
            for ref_a, ref_b in skeleton_refs:
                if _segment_on_ref(a, b, ref_a, ref_b):
                    on_skeleton = True
                    break
            lw = base_lw if on_skeleton else thin_lw
            ax.plot([a[0], b[0]], [a[1], b[1]], color="black", lw=lw, alpha=0.9, zorder=3)

    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Step 3: subdivision lines + hex tessellation like the reference image
    plot_triangle_hex_classI(aperture=4, step=3, draw_hex=True, draw_extras=True)

