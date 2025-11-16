import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.crs import Projection
import shapely.geometry as sgeom

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

# # Set image extent in degrees
img_extent = (-180, 180, -90, 90)

# Load image
img = plt.imread('./vgrid/utils/polyhedra_map/blue_marble.jpg')
# plt.imshow(img)

# Set dimensions corresponding to the Rubiks cube size
fig = plt.figure(figsize=(8.973, 6.73))

# Add subplot with the rHEALPix projection
ax = fig.add_subplot(1, 1, 1, projection=rHEALPix(), frameon=False)

ax.set_global()

# Add image
ax.imshow(img, origin='upper', extent=img_extent, transform=ccrs.PlateCarree())

# Add vertical cutting lines
x_min, x_max = ax.get_xlim()
for x in np.arange(x_min, x_max, (x_max-x_min) / 12):
    ax.axvline(x, c='w', lw=0.1)

# Add horizontal cutting lines
y_min, y_max = ax.get_ylim()
for y in np.arange(y_min, y_max, (y_max-y_min) / 9):
    ax.axhline(y, c='w', lw=0.1)

# Save figure
# plt.savefig('./blue_marble_cube.png', dpi=1800, bbox_inches='tight', pad_inches=0)
plt.show()