import json
from pathlib import Path
from typing import Optional, Union

from matplotlib.font_manager import FontProperties
import matplotlib.patches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.textpath
import matplotlib.transforms

import cartopy.crs as ccrs
from shapely.geometry import Polygon, mapping, shape
from shapely.ops import unary_union


def _path_to_geometry(path: mpath.Path):
    """Convert a Matplotlib path into a Shapely geometry preserving holes."""
    loops = []
    for polygon in path.to_polygons(closed_only=False):
        if len(polygon) < 3:
            continue
        coords = polygon.tolist()
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        loops.append(coords)

    if not loops:
        return None

    raw_polygons = []
    for coords in loops:
        poly = Polygon(coords)
        if poly.is_empty:
            continue
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty:
            continue
        raw_polygons.append(poly)

    if not raw_polygons:
        return None

    raw_polygons.sort(key=lambda p: p.area, reverse=True)

    EPSILON = 1e-9

    structures = []
    for poly in raw_polygons:
        rep_point = poly.representative_point()
        containers = [
            struct for struct in structures
            if struct["polygon"].buffer(EPSILON).contains(rep_point)
        ]

        if not containers:
            structures.append({
                "polygon": poly,
                "holes": [],
                "hole_polygons": []
            })
            continue

        # choose the smallest containing shell
        containers.sort(key=lambda struct: struct["polygon"].area)
        container = containers[0]

        depth = 1
        for hole_poly in container["hole_polygons"]:
            if hole_poly.buffer(EPSILON).contains(rep_point):
                depth += 1

        if depth % 2 == 1:
            container["holes"].append(poly.exterior.coords[:])
            container["hole_polygons"].append(poly)
        else:
            structures.append({
                "polygon": poly,
                "holes": [],
                "hole_polygons": []
            })

    polygons = []
    for struct in structures:
        shell = struct["polygon"]
        if shell.is_empty:
            continue
        hole_rings = [hole.exterior.coords[:] for hole in struct["hole_polygons"] if not hole.is_empty]
        polygon = Polygon(shell.exterior.coords[:], hole_rings)
        if not polygon.is_valid:
            polygon = polygon.buffer(0)
        if polygon.is_empty:
            continue
        polygons.append(polygon)

    if not polygons:
        return None

    geometry = unary_union(polygons)
    if geometry.is_empty:
        return None

    geometry = geometry.buffer(0)
    if geometry.is_empty:
        return None

    return geometry


def save_text_path_as_geojson(path: mpath.Path, output_path: Path, text: str):
    """Serialize a Matplotlib text path to a GeoJSON file."""
    geometry = _path_to_geometry(path)
    if geometry is None:
        raise ValueError("No valid polygons were generated from the text path.")

    if geometry.geom_type == "Polygon":
        geometries = [geometry]
    elif geometry.geom_type == "MultiPolygon":
        geometries = list(geometry.geoms)
    else:
        raise ValueError(f"Unsupported geometry type: {geometry.geom_type}")

    # Sort polygons left-to-right based on centroid x coordinate
    geometries.sort(key=lambda g: (g.centroid.x, g.centroid.y))

    features = []
    for idx, geom in enumerate(geometries):
        char = text[idx] if idx < len(text) else None
        properties = {"index": idx}
        if char is not None:
            properties["character"] = char
        features.append(
            {
                "type": "Feature",
                "properties": properties,
                "geometry": mapping(geom),
            }
        )

    feature_collection = {"type": "FeatureCollection", "features": features}

    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(feature_collection, fh, indent=2)


def main(text: str = "DGGAL", output_path: Optional[Union[str, Path]] = None):
    fig = plt.figure(figsize=[20, 10])
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    ax.coastlines()
    ax.gridlines(color='gray', linestyle='--', linewidth=0.5)
    ax.tissot(facecolor='#573681', alpha=0.2)
    ax.set_global()

    # add dggal GeoJSON overlay if available
    # geojson_path = Path(__file__).with_name("dggal.geojson")
    # if geojson_path.exists():
    #     with geojson_path.open("r", encoding="utf-8") as fh:
    #         dggal_geojson = json.load(fh)
    #     dggal_geoms = [
    #         shape(feature["geometry"])
    #         for feature in dggal_geojson.get("features", [])
    #         if "geometry" in feature
    #     ]
    #     if dggal_geoms:
    #         ax.add_geometries(
    #             dggal_geoms,
    #             ccrs.PlateCarree(),
    #             facecolor="none",
    #             edgecolor="#00008b",
    #             linewidth=1,
    #         )

    # generate a matplotlib path representing the word "cartopy"
    fp = FontProperties(family='monospace', weight='bold')
    dggs_lat = -47
    dggs_lon = -128
    bhs_lat = -47
    bhs_lon = -95
    vgrid_lat =  -47
    vgrid_lon = -108
    gitc_lat = -47
    gitc_lon = -101
    thangquach_lat = -32
    thangquach_lon = -149.6
    hanquach_lat = -32
    hanquach_lon = -123
    nhungnguyen_lat = -32
    nhungnguyen_lon = -167
    tuannguyen_lat = -32
    tuannguyen_lon = -145
    # hcmgis_lat = -47
    # hcmgis_lon = -155
    techcraft_lat = -32
    techcraft_lon = -138
    dggal_lat = -47
    dggal_lon = -148
    a5_lat = -39
    a5_lon = -150
    ngannguyen_lat = -32
    ngannguyen_lon = -151
    kimloinguyen_lat = -32
    kimloinguyen_lon = -175
    dggrid_lat = -47  
    dggrid_lon = -178
    logo_path = matplotlib.textpath.TextPath(
        (dggal_lon, dggal_lat),
        text,
        size=80,
        prop=fp,
    )

    # scale the letters up to sensible longitude and latitude sizes
    transform = matplotlib.transforms.Affine2D().scale(1, 2).translate(0, 35)

    # add a background image
    im = ax.stock_img()
    # Apply the scale transform and then the map coordinate transform
    plate_carree_transform = (transform +
                              ccrs.PlateCarree()._as_mpl_transform(ax))

    # Optionally save the transformed text path as GeoJSON
    transformed_logo_path = transform.transform_path(logo_path)
    if output_path is None:
        filename = f"{text.lower()}_logo.geojson"
        output_path = Path(__file__).with_name(filename)
    else:
        output_path = Path(output_path)
    save_text_path_as_geojson(transformed_logo_path, output_path, text=text)

    # add the path as a patch, drawing black outlines around the text
    patch = matplotlib.patches.PathPatch(logo_path,
                                         facecolor='none', edgecolor='black',
                                         transform=plate_carree_transform)
    im.set_clip_path(patch)
    ax.add_patch(patch)

    plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Generate a GeoJSON logo from text.")
    parser.add_argument(
        "text",
        nargs="?",
        default="DGGAL",
        help="Text to render (defaults to 'DGGRID').",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Optional output GeoJSON file path.",
    )
    args = parser.parse_args()
    main(text=args.text, output_path=args.output)