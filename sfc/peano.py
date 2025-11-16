# Peano curve rotated 90° CW, starting from bottom-left corner
# We first generate the curve starting at (0,0), then re-shift so the minimum y is 0 and minimum x is 0
# so the start becomes bottom-left visually.

import matplotlib.pyplot as plt
import math
import geopandas as gpd

# L-system rules
axiom = "X"
rules = {
    "X": "XFYFX+F+YFXFY-F-XFYFX",
    "Y": "YFXFY-F-XFYFX+F+YFXFY"
}

def peano_lsystem(axiom, rules, iterations):
    seq = axiom
    for _ in range(iterations):
        seq = "".join([rules.get(ch, ch) for ch in seq])
    return seq

def peano(order):
    """Plot the 2D Peano curve for a given order."""
    coords, size = peano_coords(order)
    xs, ys = zip(*coords)

    plt.figure(figsize=(6, 6))
    plt.plot(xs, ys, linewidth=0.6)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlim(-1, size)
    plt.ylim(-1, size)
    plt.xticks([])
    plt.yticks([])
    plt.gca().invert_yaxis()
    # plt.title(f"Peano Space-Filling Curve — order {order} (grid {size}×{size})")
    plt.show()

def webmercator_y_to_lat(y_norm):
    """Convert normalized Web Mercator Y coordinate to latitude."""
    return math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y_norm))))

def peano_to_lonlat(order):
    """Convert Peano curve coordinates to longitude/latitude pairs."""
    coords, size = peano_coords(order)
    lonlat = []
    for x, y in coords:
        x_norm = x / size
        y_norm = y / size
        lon = x_norm * 360.0 - 180.0
        lat = webmercator_y_to_lat(y_norm)
        lonlat.append((lon, lat))
    return lonlat

def peano_coords(order):
    """Generate Peano curve coordinates for a given order."""
    seq = peano_lsystem(axiom, rules, order)
    
    # Use step size of 1 for coordinate generation (this gives the actual curve structure)
    # Start facing down so it's rotated 90° clockwise
    angle = -90
    x, y = 0.0, 0.0
    pts = [(x, y)]

    for cmd in seq:
        if cmd == "F":
            x += math.cos(math.radians(angle)) * 1
            y += math.sin(math.radians(angle)) * 1
            pts.append((x, y))
        elif cmd == "+":
            angle += 90
        elif cmd == "-":
            angle -= 90

    # Shift so the min x and y become zero (bottom-left origin)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    min_x = min(xs)
    min_y = min(ys)
    pts_shifted = [(px - min_x, py - min_y) for px, py in pts]
    
    size = 3 ** order
    return pts_shifted, size


def peano_mollweide(order, projection="mollweide"):
    """
    Generate and plot Peano curve on a global Mollweide map.

    Args:
        order: Resolution order (e.g., 3 for 27x27 grid)
        projection: Map projection ('mollweide' or 'moll')
    """
    pts = peano_to_lonlat(order)

    lons = [math.radians(((lon + 180) % 360) - 180) for lon, lat in pts]
    lats = [math.radians(lat) for lon, lat in pts]

    # Map 'moll' to 'mollweide' for matplotlib compatibility
    if projection == "moll":
        projection = "mollweide"

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection=projection)
    ax.plot(lons, lats, linewidth=0.6)
    ax.grid(True, linestyle=":", linewidth=0.5)

    # Add world countries boundary layer (same approach as hilbert.py)
    try:
        world_countries = gpd.read_file(
            "https://raw.githubusercontent.com/opengeoshub/vopendata/refs/heads/main/shape/world_countries.geojson"
        )
        for geom in world_countries.geometry:
            if geom is None:
                continue
            if geom.geom_type == "Polygon":
                polys = [geom]
            elif geom.geom_type == "MultiPolygon":
                polys = list(geom.geoms)
            else:
                continue
            for poly in polys:
                xs, ys = zip(*poly.exterior.coords)
                lon_wrapped = [((x + 180) % 360) - 180 for x in xs]
                lon_rad = [math.radians(x) for x in lon_wrapped]
                lat_rad = [math.radians(y) for y in ys]
                ax.plot(lon_rad, lat_rad, color="black", linewidth=0.1)
    except Exception:
        pass

    plt.tight_layout()
    # out_path = f"peano_global_order{order}.png"
    # plt.savefig(out_path, dpi=200)
    plt.show()


def peano_cli():
    """
    Command-line interface for Peano curve generation and visualization.

    Usage examples:
        python peano.py --order 3 --projection mollweide
        python peano.py --order 4 --projection moll
        python peano.py --order 3
        python peano.py --mollweide --order 3
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Generate and visualize Peano space-filling curves",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python peano.py --order 3                    # Generate order 3 curve (2D view)
  python peano.py --mollweide --order 3        # Generate order 3 curve on Mollweide projection
  python peano.py --order 4 --projection moll  # Generate order 4 curve with moll projection
  python peano.py --help                       # Show this help message
        """,
    )

    parser.add_argument(
        "--order",
        "-o",
        type=int,
        default=3,
        help="Resolution order (default: 3, creates 27x27 grid)",
    )

    parser.add_argument(
        "--projection",
        "-p",
        choices=["mollweide", "moll"],
        default="mollweide",
        help="Map projection to use for mollweide view (default: mollweide)",
    )

    parser.add_argument(
        "--mollweide",
        "-m",
        action="store_true",
        help="Generate Peano curve on Mollweide projection instead of 2D view",
    )

    args = parser.parse_args()

    # Validate order
    # if args.order < 1 or args.order > 8:
    #     print(
    #         f"Warning: Order {args.order} may create very large grids. Consider using order <= 6 for reasonable performance."
    #     )
    #     response = input("Continue anyway? (y/N): ")
    #     if response.lower() != "y":
    #         print("Operation cancelled.")
    #         sys.exit(0)

    # Non-interactive execution without prompts
    # if args.mollweide:
    #     peano_mollweide(args.order, args.projection)
    # else:
    #     peano(args.order)
    peano_mollweide(3)

if __name__ == "__main__":
    # Run CLI if script is executed directly
    peano_cli()
