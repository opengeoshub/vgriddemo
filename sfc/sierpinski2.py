# Sierpinski arrowhead curve with Mollweide projection support and CLI
import math
import matplotlib.pyplot as plt


def sierpinski_instructions(order):
    """L-system string for the Sierpinski arrowhead curve."""
    s = "A"
    rules = {"A": "B-A-B", "B": "A+B+A"}
    for _ in range(order):
        s = "".join(rules[ch] if ch in rules else ch for ch in s)
    return s


def sierpinski_points(order):
    """Interpret L-system into (x,y) points (turtle graphics)."""
    instr = sierpinski_instructions(order)
    angle = math.radians(60)
    x, y = 0.0, 0.0
    dx, dy = 1.0, 0.0  # initial heading = right
    pts = [(x, y)]
    for ch in instr:
        if ch in ("A", "B"):
            x += dx
            y += dy
            pts.append((x, y))
        elif ch == "+":  # left 60°
            ndx = dx * math.cos(angle) - dy * math.sin(angle)
            ndy = dx * math.sin(angle) + dy * math.cos(angle)
            dx, dy = ndx, ndy
        elif ch == "-":  # right 60°
            ndx = dx * math.cos(-angle) - dy * math.sin(-angle)
            ndy = dx * math.sin(-angle) + dy * math.cos(-angle)
            dx, dy = ndx, ndy
    return pts


def sierpinski_coords(order, grid_size=None, orient="tl"):
    """
    Generate Sierpinski arrowhead curve coordinates.
    - order: recursion depth
    - grid_size: scale target (if None, uses 2**order)
    - orient: 'tl','bl','tr','br' places/rotates the triangle
    - returns: (scaled_points, size)
    """
    pts = sierpinski_points(order)
    xs, ys = zip(*pts)
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    width = xmax - xmin
    height = ymax - ymin

    if grid_size is None:
        grid_size = 2**order

    # scale so largest span fits grid_size-1 (0..grid_size-1)
    span = max(width, height) or 1.0
    s = (grid_size - 1) / span
    tx = -xmin * s
    ty = -ymin * s
    scaled = [(x * s + tx, y * s + ty) for (x, y) in pts]

    # optional orientation (rotate/flip to match expected placement)
    if orient == "tl":  # top-left: rotate and flip
        oriented = [(y, grid_size - 1 - x) for (x, y) in scaled]
    elif orient == "bl":  # bottom-left (default orientation)
        oriented = scaled
    elif orient == "tr":
        oriented = [(grid_size - 1 - x, y) for (x, y) in scaled]
    elif orient == "br":
        oriented = [(grid_size - 1 - y, x) for (x, y) in scaled]
    else:
        oriented = scaled

    return oriented, grid_size


def sierpinski(order, grid_size=None, orient="tl", show=True):
    """
    Generate & plot Sierpinski arrowhead curve.
    - order: recursion depth
    - grid_size: scale target (if None, uses 2**order)
    - orient: 'tl','bl','tr','br' places/rotates the triangle
    - show: whether to display the plot
    - returns: (scaled_points, bbox)
    """
    oriented, size = sierpinski_coords(order, grid_size, orient)
    
    if show:
        xs2, ys2 = zip(*oriented)
        plt.figure(figsize=(6, 6))
        plt.plot(xs2, ys2, linewidth=0.9)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.xlim(-1, size)
        plt.ylim(-1, size)
        plt.xticks([])
        plt.yticks([])
        plt.title(
            f"Sierpiński Arrowhead — order {order} (scaled to {size}×{size})"
        )
        plt.show()

    # Calculate bbox for return
    xs, ys = zip(*oriented)
    bbox = ((min(xs), min(ys)), (max(xs), max(ys)))
    return oriented, bbox


def webmercator_y_to_lat(y_norm):
    """Convert normalized Web Mercator Y coordinate to latitude."""
    return math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y_norm))))


def sierpinski_to_lonlat(order, grid_size=None, orient="tl"):
    """Convert Sierpinski curve coordinates to longitude/latitude pairs."""
    coords, size = sierpinski_coords(order, grid_size, orient)
    lonlat = []
    for x, y in coords:
        x_norm = x / size
        y_norm = y / size
        lon = x_norm * 360.0 - 180.0
        lat = webmercator_y_to_lat(y_norm)
        lonlat.append((lon, lat))
    return lonlat


def sierpinski_mollweide(order, grid_size=None, orient="tl", projection="mollweide"):
    """
    Generate and plot Sierpinski curve on a global Mollweide map.

    Args:
        order: Resolution order (e.g., 6 for 64x64 grid)
        grid_size: scale target (if None, uses 2**order)
        orient: 'tl','bl','tr','br' places/rotates the triangle
        projection: Map projection ('mollweide' or 'moll')
    """
    pts = sierpinski_to_lonlat(order, grid_size, orient)

    lons = [math.radians(((lon + 180) % 360) - 180) for lon, lat in pts]
    lats = [math.radians(lat) for lon, lat in pts]

    # Map 'moll' to 'mollweide' for matplotlib compatibility
    if projection == "moll":
        projection = "mollweide"

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection=projection)
    ax.plot(lons, lats, linewidth=0.6)
    ax.grid(True, linestyle=":", linewidth=0.5)
    ax.set_title(
        f"Sierpiński Arrowhead Curve on a Global Mollweide Map — order={order} (grid={grid_size or 2**order}×{grid_size or 2**order})"
    )

    plt.tight_layout()
    # out_path = f"sierpinski_global_order{order}.png"
    # plt.savefig(out_path, dpi=200)
    plt.show()


def sierpinski_cli():
    """
    Command-line interface for Sierpinski curve generation and visualization.

    Usage examples:
        python sierpinski2.py --order 6 --projection mollweide
        python sierpinski2.py --order 8 --projection moll
        python sierpinski2.py --order 7
        python sierpinski2.py --mollweide --order 6
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Generate and visualize Sierpinski arrowhead curves",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sierpinski2.py --order 6                    # Generate order 6 curve (2D view)
  python sierpinski2.py --mollweide --order 6        # Generate order 6 curve on Mollweide projection
  python sierpinski2.py --order 8 --projection moll  # Generate order 8 curve with moll projection
  python sierpinski2.py --help                       # Show this help message
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
        "--grid-size",
        "-g",
        type=int,
        help="Grid size override (default: 2^order)",
    )

    parser.add_argument(
        "--orientation",
        choices=["tl", "bl", "tr", "br"],
        default="tl",
        help="Triangle orientation (default: tl)",
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
        help="Generate Sierpinski curve on Mollweide projection instead of 2D view",
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

    # try:
    #     if args.mollweide:
    #         print(
    #             f"Generating Sierpinski curve for order {args.order} with {args.projection} projection..."
    #         )
    #         sierpinski_mollweide(args.order, args.grid_size, args.orientation, args.projection)
    #         print(
    #             "Sierpinski curve generation on Mollweide projection completed successfully!"
    #         )
    #     else:
    #         print(f"Generating Sierpinski curve for order {args.order} (2D view)...")
    #         sierpinski(args.order, args.grid_size, args.orientation)
    #         print("Sierpinski curve generation completed successfully!")

    # except Exception as e:
    #     print(f"Error: {e}")
    #     sys.exit(1)
    sierpinski_mollweide(6, None, orient="br", projection="mollweide")

if __name__ == "__main__":
    # Run CLI if script is executed directly
    sierpinski_cli()
