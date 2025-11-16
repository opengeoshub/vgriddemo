# Retry: Morton curve on global Mollweide map with smaller default order (6 -> 64x64 = 4096 points)
import math
import matplotlib.pyplot as plt
import geopandas as gpd


def deinterleave_bits(n: int):
    x = 0
    y = 0
    bit = 0
    while n:
        x |= (n & 1) << bit
        n >>= 1
        y |= (n & 1) << bit
        n >>= 1
        bit += 1
    return x, y


def morton_coords(order: int):
    size = 1 << order
    n_points = size * size
    coords = [deinterleave_bits(i) for i in range(n_points)]
    return coords, size


def webmercator_y_to_lat(y_norm):
    return math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y_norm))))


def morton_to_lonlat(order):
    coords, size = morton_coords(order)
    lonlat = []
    for x, y in coords:
        x_norm = x / size
        y_norm = y / size
        lon = x_norm * 360.0 - 180.0
        lat = webmercator_y_to_lat(y_norm)
        lonlat.append((lon, lat))
    return lonlat


def morton(order: int):
    """
    Plot the 2D Morton (Z-order) curve for a given order.

    Args:
        order: Resolution order (e.g., 6 for 64x64 grid)
    """
    coords, size = morton_coords(order)
    xs, ys = zip(*coords)
    plt.figure(figsize=(6, 6))
    plt.plot(xs, ys, linewidth=2)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlim(-1, size)
    plt.ylim(-1, size)
    plt.xticks([])
    plt.yticks([])
    # plt.title(f"Morton (Z-order) Curve — order {order} (grid {size}×{size})")
    plt.show()


def morton_mollweide(order: int, projection: str = "mollweide"):
    """
    Generate and plot Morton (Z-order) curve on a global Mollweide map.

    Args:
        order: Resolution order (e.g., 6 for 64x64 grid)
        projection: Map projection ('mollweide' or 'moll')
    """
    pts = morton_to_lonlat(order)

    lons = [math.radians(((lon + 180) % 360) - 180) for lon, lat in pts]
    lats = [math.radians(lat) for lon, lat in pts]

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
    # out_path = f"morton_global_order{order}.png"
    # plt.savefig(out_path, dpi=200)
    plt.show()


def morton_cli():
    """
    Command-line interface for Morton curve generation and visualization.

    Usage examples:
        python morton.py --order 6 --projection mollweide
        python morton.py --order 8 --projection moll
        python morton.py --order 7
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Generate and visualize Morton (Z-order) curves",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python morton.py --order 6                    # Generate order 6 curve with default mollweide projection
  python morton.py --order 8 --projection moll  # Generate order 8 curve with moll projection
  python morton.py --order 7                    # Generate order 7 curve with default projection
  python morton.py --help                       # Show this help message
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
        help="Generate Morton curve on Mollweide projection instead of 2D view",
    )

    args = parser.parse_args()

    # # Validate order
    # if args.order < 1 or args.order > 12:
    #     print(
    #         f"Warning: Order {args.order} may create very large grids. Consider using order <= 10 for reasonable performance."
    #     )
    #     response = input("Continue anyway? (y/N): ")
    #     if response.lower() != "y":
    #         print("Operation cancelled.")
    #         sys.exit(0)

    # Non-interactive execution without prompts
    # if args.mollweide:
    #     morton_mollweide(args.order, args.projection)
    # else:
    #     morton(args.order)
    morton_mollweide(5)


if __name__ == "__main__":
    # Run CLI if script is executed directly
    morton_cli()
