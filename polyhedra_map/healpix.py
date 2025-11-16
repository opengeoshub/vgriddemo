from __future__ import division

import matplotlib.pyplot as plt

import json
import os
from cartopy.crs import Projection
import shapely.geometry as sgeom
import cartopy.crs as ccrs
import geopandas as gpd
from shapely.ops import transform as shp_transform
import cartopy.feature as cfeature


class HEALPix(Projection):
    def __init__(self, central_longitude=0):
        proj4_params = [('proj', 'healpix'),
                        ('lon_0', central_longitude)]
        super(HEALPix, self).__init__(proj4_params)

        # Boundary is based on units of m, with a standard spherical ellipse.
        width = 2e7
        h = width/2
        box_h = width/4
 
        points = [[width, -box_h],
                  [width, box_h],
                  [width - 1*h + h/2, h],
                  [width - 1*h, box_h],
                  [width - 2*h + h / 2, h],
                  [width - 2*h, box_h],
                  [width - 3*h + h / 2, h],
                  [width - 3*h, box_h],
                  [width - 4*h + h / 2, h],
                  [width - 4*h, box_h],
                  [width - 4*h, -box_h],
                  [width - 4*h + h / 2, -h],
                  [width - 3*h, -box_h],
                  [width - 3*h + h / 2, -h],
                  [width - 2*h, -box_h],
                  [width - 2*h + h / 2, -h],
                  [width - 1*h, -box_h],
                  [width - 1*h + h / 2, -h]]
        self._boundary = sgeom.LinearRing(points)

        xs, ys = zip(*points)
        self._x_limits = min(xs), max(xs)
        self._y_limits = min(ys), max(ys)
        self._threshold = (self.x_limits[1] - self.x_limits[0]) / 1e4

    @property
    def boundary(self):
        return self._boundary

    @property
    def threshold(self):
        return self._threshold

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class rHEALPix(Projection):
    def __init__(self, central_longitude=0, north_square=0, south_square=0):
        proj4_params = [('proj', 'rhealpix'),
                        ('north_square', north_square),
                        ('south_square', south_square),
                        ('lon_0', central_longitude)]
        super(rHEALPix, self).__init__(proj4_params)

        # Boundary is based on units of m, with a standard spherical ellipse.
        nrth_x_pos = (north_square - 2) * 1e7
        sth_x_pos = (south_square - 2) * 1e7
        top = 5.05e6
        points = []
        points.extend([
                  [2e7, -5e6],
                  [2e7, top],
                  [nrth_x_pos + 1e7, top],
                  [nrth_x_pos + 1e7, 1.5e7],
                  [nrth_x_pos, 1.5e7],
                  [nrth_x_pos, top],
                  [-2e7, top]])
        if south_square != 0:
            points.append([-2e7, -top])
        points.extend([
                  [sth_x_pos, -5e6],
                  [sth_x_pos, -1.5e7],
                  [sth_x_pos + 1e7, -1.5e7],
                  [sth_x_pos + 1e7, -5e6],
                  ])
        self._boundary = sgeom.LineString(points[::-1])

        xs, ys = zip(*points)
        self._x_limits = min(xs), max(xs)
        self._y_limits = min(ys), max(ys)
        self._threshold = (self.x_limits[1] - self.x_limits[0]) / 1e4

    @property
    def boundary(self):
        return self._boundary

    @property
    def threshold(self):
        return self._threshold

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


crs = ccrs.Mollweide()
# crs = HEALPix(central_longitude=50)
# crs = rHEALPix(central_longitude=50, north_square=0, south_square=0)
ax = plt.axes(projection=crs)
ax.stock_img()
ax.coastlines()
ax.gridlines(color='gray', linestyle='--', linewidth=1)
ax.tissot(color='red', alpha=0.5)
ax.set_global()
# ax.gridlines(color='gray', linestyle='--', linewidth=0.5)
# geojson_path = os.path.join(os.path.dirname(__file__), 'rhealpix_0.geojson')
# if os.path.exists(geojson_path):
#     gdf = gpd.read_file(geojson_path)
#     gdf.boundary.to_crs(crs).plot(
#         color=None, edgecolor="red", linewidth=1, ax=ax
#     )

# # Optionally add centroid labels from centroids.geojson (label only, no points)
# centroids_path = os.path.join(os.path.dirname(__file__), 'centroids.geojson')
# if os.path.exists(centroids_path):
#     centroids_gdf = gpd.read_file(centroids_path)
#     for _, row in centroids_gdf.iterrows():
#         geom = row.get('geometry')
#         if geom is None or geom.geom_type != 'Point':
#             continue
#         label = row.get('dggal_rhealpix')
#         if not label:
#             continue
#         lon, lat = geom.x, geom.y
#         ax.text(
#             lon,
#             lat,
#             label,
#             transform=ccrs.PlateCarree(),
#             color='red',
#             fontsize=20,
#             ha='center',
#             va='center',
#             zorder=10
#         )

plt.show()

## Removed standalone fallback figure; overlay is drawn on the Healpix axes above
