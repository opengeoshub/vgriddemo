import matplotlib.pyplot as plt
from cartopy.crs import Projection
import cartopy.crs as ccrs
import shapely.geometry as sgeom


class rHEALPix(Projection):
    def __init__(self, central_longitude=50, north_square=0, south_square=0):
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

import matplotlib.pyplot as plt
# ax = plt.axes(projection=rHEALPix())
ax = plt.axes(projection=ccrs.TransverseMercator())
ax.stock_img(name='ne_shaded')
ax.coastlines()
ax.gridlines()
ax.set_global()
plt.show()