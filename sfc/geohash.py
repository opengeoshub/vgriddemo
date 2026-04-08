"""
Plot a Plate Carree-style world map (lon/lat) and draw split lines with labels.

`plot_world` builds the base map. `plot_split0` draws split_lon + 0/1 labels.
`plot_split1` draws the same split_lon meridian, the equator on the western half,
and top/bottom labels (default 01 north, 00 south). `plot_split2` repeats the meridian
and western equator from split0/split1, then splits the northern western cell into 010 / 011.
`plot_split3` reuses split2 lines (no 010/011 text) plus a horizontal cut through 011: 0111 (north) / 0110 (south).
`plot_split4` adds a vertical cut through 0111 with labels 01110 / 01111 only (no plot_split3 0110/0111 text).

"""

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import cartopy.crs as ccrs

DIVIDER_COLOR = "red"
DEFAULT_FIGSIZE = (5.8, 3.1)
STOCK_IMG_ALPHA = 0.2
# Coastlines use the same alpha as the stock image for a consistent faded base map.
COASTLINE_ALPHA = STOCK_IMG_ALPHA
# Same default width for map outline and all horizontal/vertical split lines.
SPLIT_LINEWIDTH = 4.0


def _style_map_boundary(
    ax,
    *,
    edgecolor: str = DIVIDER_COLOR,
    linewidth: float = SPLIT_LINEWIDTH,
) -> None:
    """Cartopy map outline; use the same linewidth as split lines."""
    geo = ax.spines.get("geo")
    if geo is not None:
        geo.set_edgecolor(edgecolor)
        geo.set_linewidth(linewidth)
        geo.set_visible(True)
    outline = getattr(ax, "outline_patch", None)
    if outline is not None:
        outline.set_edgecolor(edgecolor)
        outline.set_linewidth(linewidth)


def _plot_split_lon_meridian(
    ax,
    split_lon: float,
    *,
    linewidth: float,
) -> None:
    """Full-height vertical meridian in Plate Carrée (world edges)."""
    ax.plot(
        [split_lon, split_lon],
        [-90.0, 90.0],
        color=DIVIDER_COLOR,
        linewidth=linewidth,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )


def _draw_binary_label_colored(
    ax,
    x: float,
    y: float,
    label: str,
    *,
    fontsize: float,
    zorder: int = 6,
    digit_colors: dict[str, str] | None = None,
    colors_by_index: list[str] | None = None,
) -> None:
    """
    Draw each character of `label` separately, so digits can be color-coded.
    Intended for strings like "01", "00", etc.
    """
    s = str(label)
    if not s:
        return

    if digit_colors is None:
        # Fallback if a caller doesn't provide colors_by_index.
        digit_colors = {"0": "red", "1": "blue"}

    # Requested behavior: interleave colors by character position.
    # Index 0 -> red, index 1 -> blue, index 2 -> red, ...
    if colors_by_index is None:
        colors_by_index = ["red" if (i % 2 == 0) else "blue" for i in range(len(s))]

    # Draw each character using "offset points" so spacing doesn't depend on map scaling.
    # This also avoids the previous transform-composition issues where text could disappear.
    xycoords = ccrs.PlateCarree()._as_mpl_transform(ax)

    # Spacing between digit centers in points.
    # For typical fontsize values (e.g. 70), this prevents overlap for 2-character labels.
    dx_points = 0.65 * float(fontsize)
    start_points = -dx_points * (len(s) - 1) / 2.0

    for i, ch in enumerate(s):
        xoff = start_points + i * dx_points
        color = (
            colors_by_index[i]
            if colors_by_index is not None and i < len(colors_by_index)
            else digit_colors.get(ch, DIVIDER_COLOR)
        )
        ax.annotate(
            ch,
            xy=(x, y),
            xycoords=xycoords,
            textcoords="offset points",
            xytext=(xoff, 0.0),
            ha="center",
            va="center",
            fontsize=fontsize,
            color=color,
            zorder=zorder,
        )


def plot_world(figsize: tuple[float, float] = DEFAULT_FIGSIZE) -> tuple[plt.Figure, Axes]:
    """Create figure/axes and draw the global base map (no split lines)."""
    fig, ax = plt.subplots(
        figsize=figsize,
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    ax.set_global()
    ax.coastlines(alpha=COASTLINE_ALPHA)
    # ax.gridlines(color="gray", linestyle="--", linewidth=0.5)
    # ax.set_xticks([])
    # ax.set_yticks([])
    ax.stock_img()
    for im in ax.get_images():
        im.set_alpha(STOCK_IMG_ALPHA)
    # Match split line width unless a plot_split* overrides with its own linewidth.
    _style_map_boundary(ax, linewidth=SPLIT_LINEWIDTH)
    return fig, ax


def plot_split0(
    ax,
    split_lon: float = 0.0,
    left_label: str = "0",
    right_label: str = "1",
    *,
    linewidth: float = SPLIT_LINEWIDTH,
    fontsize: float = 80,
) -> None:
    """Vertical meridian at split_lon; label left/right halves."""
    _style_map_boundary(ax, linewidth=linewidth)
    split_lon = float(split_lon)
    _plot_split_lon_meridian(ax, split_lon, linewidth=linewidth)
    left_x = (-180.0 + split_lon) / 2.0
    right_x = (split_lon + 180.0) / 2.0
    _draw_binary_label_colored(
        ax,
        left_x,
        0.0,
        left_label,
        fontsize=fontsize,
        zorder=6,
    )
    _draw_binary_label_colored(
        ax,
        right_x,
        0.0,
        right_label,
        fontsize=fontsize,
        zorder=6,
    )


def plot_split1(
    ax,
    split_lon: float = 0.0,
    top_label: str = "01",
    bottom_label: str = "00",
    *,
    linewidth: float = SPLIT_LINEWIDTH,
    fontsize: float = 70,
) -> None:
    """
    Same vertical split_lon meridian as plot_split0, plus horizontal split on the left half
    at the equator. Default labels: top (north) \"01\", bottom (south) \"00\".
    """
    _style_map_boundary(ax, linewidth=linewidth)
    split_lon = float(split_lon)
    _plot_split_lon_meridian(ax, split_lon, linewidth=linewidth)
    ax.plot(
        [-180.0, split_lon],
        [0.0, 0.0],
        color=DIVIDER_COLOR,
        linewidth=linewidth,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )
    left_x = (-180.0 + split_lon) / 2.0
    # Centers of northern and southern sub-rectangles in the left half
    y_top = 45.0
    y_bottom = -45.0
    _draw_binary_label_colored(
        ax,
        left_x,
        y_top,
        top_label,
        fontsize=fontsize,
        zorder=6,
    )
    _draw_binary_label_colored(
        ax,
        left_x,
        y_bottom,
        bottom_label,
        fontsize=fontsize,
        zorder=6,
    )


def plot_split2(
    ax,
    split_lon: float = 0.0,
    left_label: str = "010",
    right_label: str = "011",
    *,
    linewidth: float = SPLIT_LINEWIDTH,
    fontsize: float = 60,
    draw_labels: bool = True,
) -> None:
    """
    Draws the same geometry as plot_split0 + plot_split1 (full-height meridian at split_lon,
    equator segment on the western half), then subdivides the northern western cell (01)
    with a vertical at the longitude midpoint. Optionally labels 010 (west) / 011 (east).
    """
    _style_map_boundary(ax, linewidth=linewidth)
    split_lon = float(split_lon)
    # As in plot_split0 / plot_split1: world meridian at split_lon
    _plot_split_lon_meridian(ax, split_lon, linewidth=linewidth)
    # As in plot_split1: equator from dateline side of the western half to split_lon
    ax.plot(
        [-180.0, split_lon],
        [0.0, 0.0],
        color=DIVIDER_COLOR,
        linewidth=linewidth,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )
    mid_lon = (-180.0 + split_lon) / 2.0
    ax.plot(
        [mid_lon, mid_lon],
        [0.0, 90.0],
        color=DIVIDER_COLOR,
        linewidth=linewidth,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )
    if draw_labels:
        x_w = (-180.0 + mid_lon) / 2.0
        x_e = (mid_lon + split_lon) / 2.0
        y_c = 45.0
        _draw_binary_label_colored(
            ax,
            x_w,
            y_c,
            left_label,
            fontsize=fontsize,
            zorder=6,
        )
        _draw_binary_label_colored(
            ax,
            x_e,
            y_c,
            right_label,
            fontsize=fontsize,
            zorder=6,
        )


def plot_split3(
    ax,
    split_lon: float = 0.0,
    top_label: str = "0111",
    bottom_label: str = "0110",
    *,
    linewidth: float = SPLIT_LINEWIDTH,
    fontsize_011_cells: float = 50,
    draw_labels: bool = True,
) -> None:
    """
    Same lines as plot_split2 but without 010/011 labels, then subdivide cell 011 at lat 45°:
    top 0111, bottom 0110 (longitude span from mid_lon to split_lon). Set draw_labels=False
    to omit 0110/0111 text (e.g. for plot_split4).
    """
    plot_split2(
        ax,
        split_lon=split_lon,
        linewidth=linewidth,
        draw_labels=False,
    )
    split_lon = float(split_lon)
    mid_lon = (-180.0 + split_lon) / 2.0
    lat_split = 45.0
    ax.plot(
        [mid_lon, split_lon],
        [lat_split, lat_split],
        color=DIVIDER_COLOR,
        linewidth=linewidth,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )
    x_011 = (mid_lon + split_lon) / 2.0
    y_top = 67.5
    y_bottom = 22.5
    if draw_labels:
        _draw_binary_label_colored(
            ax,
            x_011,
            y_top,
            top_label,
            fontsize=fontsize_011_cells,
            zorder=6,
        )
        _draw_binary_label_colored(
            ax,
            x_011,
            y_bottom,
            bottom_label,
            fontsize=fontsize_011_cells,
            zorder=6,
        )


def plot_split4(
    ax,
    split_lon: float = 0.0,
    left_label: str = "01110",
    right_label: str = "01111",
    *,
    linewidth: float = SPLIT_LINEWIDTH,
    fontsize_0111_children: float = 35,
) -> None:
    """
    Lines through plot_split3 with plot_split3 labels omitted (draw_labels=False), then split the
    northern 011 half at lon x_011. Only 01110 / 01111 are drawn — no 0110/0111 text from split3.
    """
    plot_split3(
        ax,
        split_lon=split_lon,
        linewidth=linewidth,
        draw_labels=False,
    )
    split_lon = float(split_lon)
    mid_lon = (-180.0 + split_lon) / 2.0
    x_011 = (mid_lon + split_lon) / 2.0
    ax.plot(
        [x_011, x_011],
        [45.0, 90.0],
        color=DIVIDER_COLOR,
        linewidth=linewidth,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )
    x_w = (mid_lon + x_011) / 2.0
    x_e = (x_011 + split_lon) / 2.0
    y_c = 67.5
    _draw_binary_label_colored(
        ax,
        x_w,
        y_c,
        left_label,
        fontsize=fontsize_0111_children,
        zorder=6,
    )
    _draw_binary_label_colored(
        ax,
        x_e,
        y_c,
        right_label,
        fontsize=fontsize_0111_children,
        zorder=6,
    )


def main() -> None:
    fig, ax = plot_world()
    # plot_split0(ax)
    # plot_split1(ax)
    # plot_split2(ax)
    # plot_split3(ax)
    plot_split4(ax)
    plt.tight_layout(pad=0.2)
    plt.show()


if __name__ == "__main__":
    main()
