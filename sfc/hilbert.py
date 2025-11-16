import matplotlib.pyplot as plt
import math
import geopandas as gpd
from shapely.geometry import LineString


def hilbert_index_to_xy(index, order):
    """Convert Hilbert index to (x, y)."""
    x = y = 0
    n = 1 << order
    idx = index
    s = 1
    while s < n:
        rx = 1 & (idx // 2)
        ry = 1 & (idx ^ rx)
        if ry == 0:
            if rx == 1:
                x, y = s - 1 - x, s - 1 - y
            x, y = y, x
        x += s * rx
        y += s * ry
        idx //= 4
        s *= 2
    return x, y


def hilbert(order):
    """Plot the 2D Hilbert curve for a given order."""
    n_points = 1 << (2 * order)
    pts = [hilbert_index_to_xy(i, order) for i in range(n_points)]
    size = 1 << order

    xs, ys = zip(*pts)
    plt.figure(figsize=(6, 6))
    plt.plot(xs, ys, linewidth=2)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlim(-1, size)
    plt.ylim(-1, size)
    plt.xticks([])
    plt.yticks([])
    # plt.title(f"Hilbert Space-Filling Curve — order {order} (grid {size}×{size})")
    plt.show()


def webmercator_y_to_lat(y_norm):
    """Convert normalized Web Mercator Y coordinate to latitude."""
    return math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y_norm))))


def hilbert_to_lonlat(order):
    """Convert Hilbert curve coordinates to longitude/latitude pairs."""
    coords, size = hilbert_coords(order)
    lonlat = []
    for x, y in coords:
        x_norm = x / size
        y_norm = y / size
        lon = x_norm * 360.0 - 180.0
        lat = webmercator_y_to_lat(y_norm)
        lonlat.append((lon, lat))
    return lonlat


def hilbert_coords(order):
    """Generate Hilbert curve coordinates for a given order."""
    n_points = 1 << (2 * order)
    coords = [hilbert_index_to_xy(i, order) for i in range(n_points)]
    size = 1 << order
    return coords, size


def hilbert_mollweide(order, projection="mollweide"):
    """
    Generate and plot Hilbert curve on a global Mollweide map.

    Args:
        order: Resolution order (e.g., 6 for 64x64 grid)
        projection: Map projection ('mollweide' or 'moll')
    """
    pts = hilbert_to_lonlat(order)

    lons = [math.radians(((lon + 180) % 360) - 180) for lon, lat in pts]
    lats = [math.radians(lat) for lon, lat in pts]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection=projection)
    ax.plot(lons, lats, linewidth=1)
    ax.grid(True, linestyle=":", linewidth=0.5)

    # Add world countries boundary layer
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
    # out_path = f"hilbert_global_order{order}.png"
    # plt.savefig(out_path, dpi=200)
    plt.show()


def hilbert_json(order, output_path=None):
    """
    Generate Hilbert curve and save as GeoJSON in WGS84 CRS.
    
    Args:
        order: Resolution order (e.g., 6 for 64x64 grid)
        output_path: Output file path (default: hilbert_order{order}.geojson)
    
    Returns:
        str: Path to the saved GeoJSON file
    """
    # Generate Hilbert curve coordinates in lon/lat
    pts = hilbert_to_lonlat(order)
    
    # Create LineString geometry from the points
    line = LineString(pts)
    
    # Create GeoDataFrame with WGS84 CRS
    gdf = gpd.GeoDataFrame(
        [{"order": order, "geometry": line}],
        crs="EPSG:4326"  # WGS84
    )
    
    # Set output path if not provided
    if output_path is None:
        output_path = f"hilbert_order{order}.geojson"
    
    # Save as GeoJSON
    gdf.to_file(output_path, driver="GeoJSON")
    
    print(f"Hilbert curve saved to: {output_path}")
    print(f"Order: {order}, Points: {len(pts)}, CRS: WGS84 (EPSG:4326)")
    
    return output_path


def hilbert_cli():
    """
    Command-line interface for Hilbert curve generation and visualization.

    Usage examples:
        python hilbert.py --order 6 --projection mollweide
        python hilbert.py --order 8 --projection moll
        python hilbert.py --order 7
        python hilbert.py --mollweide --order 6
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Generate and visualize Hilbert space-filling curves",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python hilbert.py --order 6                    # Generate order 6 curve (2D view)
  python hilbert.py --mollweide --order 6        # Generate order 6 curve on Mollweide projection
  python hilbert.py --order 8 --projection moll  # Generate order 8 curve with moll projection
  python hilbert.py --json --order 6             # Export order 6 curve as GeoJSON
  python hilbert.py --json --order 6 --output my_curve.geojson  # Custom output path
  python hilbert.py --help                       # Show this help message
        """,
    )

    parser.add_argument(
        "--order",
        "-o",
        type=int,
        default=6,
        help="Resolution order (default: 6, creates 64x64 grid)",
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
        help="Generate Hilbert curve on Mollweide projection instead of 2D view",
    )

    parser.add_argument(
        "--json",
        "-j",
        action="store_true",
        help="Export Hilbert curve as GeoJSON file in WGS84 CRS",
    )

    parser.add_argument(
        "--output",
        "-output",
        type=str,
        help="Output file path for GeoJSON export (default: hilbert_order{order}.geojson)",
    )

    args = parser.parse_args()

    # Validate order
    # if args.order < 1 or args.order > 12:
    #     print(
    #         f"Warning: Order {args.order} may create very large grids. Consider using order <= 10 for reasonable performance."
    #     )
    #     response = input("Continue anyway? (y/N): ")
    #     if response.lower() != "y":
    #         print("Operation cancelled.")
    #         sys.exit(0)

    # Non-interactive execution without prompts
    # if args.json:
    #     hilbert_json(args.order, args.output)
    # elif args.mollweide:
    #     hilbert_mollweide(args.order, args.projection)
    # else:
    #     hilbert(args.order)
    hilbert_json(4)

if __name__ == "__main__":
    # Run CLI if script is executed directly
    hilbert_cli()
