""" This code generates Voronoi polygons on a sphere under the consideration of a curved planetary surface. It automatically considers 
polar and Date Line intersections. The generated polygons are saved in the GeoJSON format. 

Global Voronoi polygons can be constructed from randomly generated points or a list of geographic lat/lon point coordinates. If you want to use a set of 
point coordinates to generate Voronoi polygons, you modify the input_samples list (Format: [[lon1, lat1],[lon_2, lat_2],[lon_n, lat_n]]) and 
set use_random_set_of_input_points = False.

Dependencies:
- numpy
- scipy (SphericalVoronoi)
- shapely
- pyproj

The following variables need to be changed (lines 700 ff):

use_random_set_of_input_points - True if using randomly distributed points to generate Voronoi polygons, False if using your own point coordinates (20 minimum)
number_of_random_points - Number of random points that should be generated, if using randomly distributed points (>= 20 points required)
input_samples - Point coordinates to generate Voronoi polygons from, if using own point coordinates (>= 20 points required)
input_samples_point_shapefile_output_path - Path to save the input points GeoJSON: 'C:\path\to\point_file.shp' (will be saved as .geojson)
voronoi_polygons_shapefile_output_path - Path to save the Voronoi polygon GeoJSON: 'C:\path\to\polygon_file.shp' (will be saved as .geojson)
geogr_sr_text - Geographic coordinate system WKT, defining the reference body the data is located on 

Minimum of 20 input points is required because I encoutered difficulties when using the SphericalVoronoi functions from the scipy library with 
less data points. """

import numpy, math, re, random, json
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, mapping
from scipy.spatial import SphericalVoronoi
from pyproj import Proj, Transformer, CRS

def direct_vincenty(flattening, major_axis, vertices_angle_list, dist_buffer, bufferfactor): 
	global buffer_vertices_list
	buffer_vertices_list = []
	for pair in vertices_angle_list:
		
		""" Calculation of Point 1 (phi) on auxiliary sphere, azimuth of geodesic at equator and length of geodesic between 
		equator and Point 1. """
		
		f = flattening
		a = major_axis
		phi1 = pair[1]
		lambda1 = pair[0]		
		alpha12 = pair[2]
		s = bufferfactor * dist_buffer
		piD4 = math.atan(1.0) 
		two_pi = piD4 * 8.0 
		phi1 = phi1 * piD4 / 45.0 
		lambda1 = lambda1 * piD4 / 45.0 
		alpha12 = alpha12 * piD4 / 45.0 
		if alpha12 < 0.0: 
			alpha12 = alpha12 + two_pi 
		if alpha12 > two_pi: 
			alpha12 = alpha12 - two_pi
		b = a * (1.0 - f) 
		tanU1 = (1-f) * math.tan(phi1) 
		U1 = math.atan(tanU1) 
		sigma1 = math.atan2(tanU1, math.cos(alpha12)) 
		Sinalpha0 = math.cos(U1) * math.sin(alpha12) 
		cosalpha0_sq = 1.0 - Sinalpha0 * Sinalpha0 
		u_sq = cosalpha0_sq * (a * a - b * b ) / (b * b) 
		A = 1.0 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * \
			(320 - 175 * u_sq))) 
		B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))) 
		sigma = (s / (b * A)) 
		last_sigma = -9999999.9
		
		""" Approximation of sigma (distance Point1-Point2 on auxiliary sphere). """
		
		while abs(last_sigma - sigma) > 1.0e-12:
			two_sigma_m = 2 * sigma1 + sigma 
			delta_sigma = B * math.sin(sigma) * (math.cos(two_sigma_m) + (B/4) * (math.cos(sigma) * \
				(-1 + 2 * math.pow(math.cos(two_sigma_m), 2) - (B/6) * math.cos(two_sigma_m) * \
				(-3 + 4 * math.pow(math.sin(sigma), 2 )) * (-3 + 4 * math.pow(math.cos(two_sigma_m), 2)))))
			last_sigma = sigma 
			sigma = (s / (b * A)) + delta_sigma 
		
		""" Calculation of Point 2 coordinates on ellipsoid. """
		
		phi2 = math.atan2 ((math.sin(U1) * math.cos(sigma) + math.cos(U1) * math.sin(sigma) * math.cos(alpha12)), \
			((1-f) * math.sqrt(math.pow(Sinalpha0, 2) + pow(math.sin(U1) * math.sin(sigma) - math.cos(U1) * \
			math.cos(sigma) * math.cos(alpha12), 2))))
		lambda_new = math.atan2((math.sin(sigma) * math.sin(alpha12)), (math.cos(U1) * math.cos(sigma) -  \
			math.sin(U1) *  math.sin(sigma) * math.cos(alpha12))) 
		C = (f/16) * cosalpha0_sq * (4 + f * (4 - 3 * cosalpha0_sq)) 
		L = lambda_new - (1-C) * f * Sinalpha0 * (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) + \
			C * math.cos(sigma) * (-1 + 2 * math.pow(math.cos(two_sigma_m), 2)))) 
		lambda2 = lambda1 + L 
		alpha21 = math.atan2 (Sinalpha0, (-math.sin(U1) * math.sin(sigma) + math.cos(U1) * math.cos(sigma) * math.cos(alpha12))) 
		alpha21 = alpha21 + two_pi / 2.0 # backwards azimuth
		
		if alpha21 < 0.0: 
			alpha21 = alpha21 + two_pi 
		if alpha21 > two_pi: 
			alpha21 = alpha21 - two_pi 
		
		phi2 = phi2 * 45.0 / piD4 
		lambda2 = lambda2 * 45.0 / piD4 
		alpha21 = alpha21 * 45.0 / piD4 
		
		buffer_vertices_list.append([lambda2, phi2, pair[3], pair[4]])
		
def inverse_vincenty(flattening, major_axis, phi1, lambda1, phi2, lambda2):
	global geodesic_distance_crater_area, geodesic_distance_crater_area2, direction12
	
	""" Calculation of Points 1 and 2 (phi) on auxiliary sphere and difference in latitude. """
	
	a = major_axis
	f = flattening
	piD4 = math.atan(1.0)
	two_pi = piD4 * 8.0
	phi1 = phi1 * piD4 / 45.0
	lambda1 = lambda1 * piD4 / 45.0
	phi2 = phi2 * piD4 / 45.0
	lambda2 = lambda2 * piD4 / 45.0
	
	b = a * (1.0 - f)
	TanU1 = (1-f) * math.tan(phi1)
	TanU2 = (1-f) * math.tan(phi2)
	U1 = math.atan(TanU1)
	U2 = math.atan(TanU2)
	L = lambda2 - lambda1
	lambda_new = L
	last_lembda = -9999999.9

	""" Approximation of lambda_new (difference in longitude on auxiliary sphere, sigma (distance Point1Point2 on auxiliary sphere) and alpha0 
	(azimuth of geodesic at the equator). """
	
	while abs(last_lembda - lambda_new) > 1.0e-12:
		sqr_sin_sigma = pow(math.cos(U2) * math.sin(lambda_new), 2) + pow((math.cos(U1) * math.sin(U2) - \
			math.sin(U1) *  math.cos(U2) * math.cos(lambda_new)), 2)
		Sin_sigma = math.sqrt(sqr_sin_sigma)
		Cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * math.cos(lambda_new)
		sigma = math.atan2(Sin_sigma, Cos_sigma)
		if sigma == 0:
			sigma += 0.0000000001 # to avoid division by zero with Sin_alpha0 computation
		Sin_alpha0 = math.cos(U1) * math.cos(U2) * math.sin(lambda_new) / math.sin(sigma)
		alpha0 = math.asin(Sin_alpha0)
		Cos2sigma_m = math.cos(sigma) - (2 * math.sin(U1) * math.sin(U2) / pow(math.cos(alpha0), 2))
		C = (f/16) * pow(math.cos(alpha0), 2) * (4 + f * (4 - 3 * pow(math.cos(alpha0), 2)))
		last_lembda = lambda_new
		lambda_new = L + (1-C) * f * math.sin(alpha0) * (sigma + C * math.sin(sigma) * \
			(Cos2sigma_m + C * math.cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2))))
		
	u2 = pow(math.cos(alpha0), 2) * (a * a - b * b) / (b * b)
	A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
	B = (u2/1024) * (256 + u2 * (-128+ u2 * (74 - 47 * u2)))
	delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2)) - \
		(B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * (-3 + 4 * pow(Cos2sigma_m, 2))))
	
	""" Calculation of distance and azimuth on ellipsoid. """
	
	s = b * A * (sigma - delta_sigma)
	alpha12 = math.atan2((math.cos(U2) * math.sin(lambda_new)), \
		(math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(U2) * math.cos(lambda_new)))
	alpha21 = math.atan2((math.cos(U1) * math.sin(lambda_new)), \
		(-math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(U2) * math.cos(lambda_new)))

	if alpha12 < 0.0: 
		alpha12 =  alpha12 + two_pi
	if alpha12 > two_pi: 
		alpha12 = alpha12 - two_pi
		
	alpha21 = alpha21 + two_pi / 2.0 # backwards azimuth
	if alpha21 < 0.0: 
		alpha21 = alpha21 + two_pi
	if alpha21 > two_pi: 
		alpha21 = alpha21 - two_pi

	alpha12 = alpha12 * 45.0 / piD4
	alpha21 = alpha21 * 45.0 / piD4
	geodesic_distance_crater_area = s
	geodesic_distance_crater_area2 = s
	direction12 = alpha12
	
def randomness_generate_random_points(no_of_random_points):
	global random_points
	
	random_points = []
	
	minimum_lat = -90
	maximum_lat = 90
	minimum_lon = -179.99
	maximum_lon = 179.99 
	
	""" Generate random points on a sphere (in an envelope around the reference areas). If the random point intersects the reference ares, 
	it is considered for randomness analysis. """
	
	while len(random_points) < no_of_random_points:		
		random_lon = random.uniform(minimum_lon * math.pi / 180, maximum_lon * math.pi / 180)
		random_lat = random.uniform((math.sin(math.radians(minimum_lat)) + 1)/2, (math.sin(math.radians(maximum_lat)) + 1)/2)
		
		random_lon = random_lon * 180 / math.pi
		random_lat = math.degrees(math.asin(2 * random_lat - 1))
		
		""" Add random point if it intersects the initial reference area(s) """
			
		random_points.append([random_lon, random_lat]) # random points array has the same structure as the bins array
			
def main(use_random_set_of_input_points, number_of_random_points, input_samples, input_samples_point_shapefile_output_path, voronoi_polygons_shapefile_output_path, geogr_sr_text):
	
	voronoi_polygon_areas = []
	
	""" We use +- 179.999999 rather than +- 180 deg lon to indicate a date line since -180 deg lon is sometimes automatically 
	converted to +180 deg lon by GDAL/OGR. """
	
	if use_random_set_of_input_points == True:
	
		""" Generate random Points """ 
		
		randomness_generate_random_points(number_of_random_points)
		
		input_samples = random_points
	
	if len(input_samples) < 20:
		print("Need minimum of 20 points for the construction of Voronoi polygons.")
		exit()
	
	input_craters_scipy = numpy.array([])
	
	# Parse WKT to extract ellipsoid parameters
	voronoi_out_polygons_crs = CRS.from_wkt(geogr_sr_text)
	major_axis = voronoi_out_polygons_crs.ellipsoid.semi_major_metre
	minor_axis = voronoi_out_polygons_crs.ellipsoid.semi_minor_metre
	
	# Prepare GeoJSON structures
	points_features = []
	polygon_features = []
	
	""" Project geographic coordinates to coordinates in 3D space [-1,1] """
	
	for input_crater in input_samples:
		input_crater_X = input_crater[0]
		input_crater_Y = input_crater[1]	
	
		# Add point to GeoJSON features
		point_geom = Point(input_crater_X, input_crater_Y)
		points_features.append({
			"type": "Feature",
			"geometry": mapping(point_geom),
			"properties": {}
		})
	
		input_crater_X = math.pi*input_crater_X/180 # Convert to radians
		input_crater_Y = math.pi*input_crater_Y/180
		
		input_crater_Y -= 1.570795765134 # subtract 90 degrees (in radians)
		
		input_crater_X_3D = math.sin(input_crater_Y) * math.cos(input_crater_X)
		input_crater_Z_3D = math.cos(input_crater_Y)
		input_crater_Y_3D = math.sin(input_crater_Y) * math.sin(input_crater_X)
		
		input_craters_scipy = numpy.append(input_craters_scipy, [[input_crater_X_3D, input_crater_Y_3D, input_crater_Z_3D]])
	
	input_craters_scipy.resize((len(input_samples), 3)) # reshape numpy array
	
	""" Calculate Voronoi diagrams on a sphere in 3D space. This only generates vertices on the sphere. Polygon edges cut the sphere's interior. """
		
	sphere_center_coordinates = numpy.array([0, 0, 0])
	sphere_radius = 1
	
	spherical_voronoi = SphericalVoronoi(input_craters_scipy, sphere_radius, sphere_center_coordinates)
	
	spherical_voronoi.sort_vertices_of_regions()
	
	""" Main part of geodesic Voronoi polygon generation: Calculate spherical coordinates of each Voronoi Polygon and save geometries. """
	
	voronoi_count = 0
	
	for spherical_voronoi_region_3D in spherical_voronoi.regions:
		
		voronoi_vertex_count = 1
		voronoi_polygon_area = 0
		
		# List to collect polygon parts (for multi-polygon support)
		voronoi_polygon_parts = []
		
		""" Calculate geographic coordinates of polygon vertices """
		
		lons = []
		lats = []
		
		for spherical_voronoi_vertex_3D in spherical_voronoi.vertices[spherical_voronoi_region_3D]:
			spherical_voronoi_vertex_3D_X = spherical_voronoi_vertex_3D[0]
			spherical_voronoi_vertex_3D_Y = spherical_voronoi_vertex_3D[1] 
			spherical_voronoi_vertex_3D_Z = spherical_voronoi_vertex_3D[2]
			
			spherical_voronoi_vertex_lon = math.degrees(math.atan(spherical_voronoi_vertex_3D_Y/spherical_voronoi_vertex_3D_X))
			
			""" Account for atan ambiguity. Otherwise, polygon vertices would all be on the nearside. """
			
			if spherical_voronoi_vertex_3D_X > 0:
				spherical_voronoi_vertex_lon -= 180
				if spherical_voronoi_vertex_lon < -180:
					spherical_voronoi_vertex_lon += 360
					
			spherical_voronoi_vertex_lat = math.degrees(math.asin(spherical_voronoi_vertex_3D_Z))
			
			""" Add vertex to lons[] and lats[] """
	
			lons.append(spherical_voronoi_vertex_lon)
			lats.append(spherical_voronoi_vertex_lat)
			
			if voronoi_vertex_count == 1:
				spherical_voronoi_end_vertex_lon = spherical_voronoi_vertex_lon
				spherical_voronoi_end_vertex_lat = spherical_voronoi_vertex_lat
			
			voronoi_vertex_count += 1
		
		""" Close ring """
		
		lons.append(spherical_voronoi_end_vertex_lon)
		lats.append(spherical_voronoi_end_vertex_lat)
		
		""" define reprojections from input craters (center of voronoi polygon) """
		
		center_voronoi_polygon_X = input_samples[voronoi_count][0]
		center_voronoi_polygon_Y = input_samples[voronoi_count][1]
		
		# Create Lambert Azimuthal Equal Area projection centered on the polygon center
		voronoi_LAEA_proj = Proj(proj='laea', lat_0=center_voronoi_polygon_Y, lon_0=center_voronoi_polygon_X, 
		                          x_0=0, y_0=0, a=major_axis, b=minor_axis, units='m')
		
		# Create transformer for coordinate conversions
		geographic_proj = Proj(proj='latlong', a=major_axis, b=minor_axis)
		transformer_to_laea = Transformer.from_proj(geographic_proj, voronoi_LAEA_proj, always_xy=True)
		transformer_from_laea = Transformer.from_proj(voronoi_LAEA_proj, geographic_proj, always_xy=True)
		
		""" define polygons that eventually become a voronoi polygon """
		
		# List to hold ring coordinates
		voronoi_ring_coords = []
		voronoi_ring_dateline_intersection_coords = []
		
		flattening = (major_axis-minor_axis) / major_axis
		
		""" Use vertices of voronoi polygons to calculate geodesic polygon edges. """
		
		voronoi_vertex_lat_lon_count = 0
		voronoi_date_line_intersection_count = 0
		
		for voronoi_vertex_lon in lons:
			voronoi_vertex_lat = lats[voronoi_vertex_lat_lon_count]
			
			if voronoi_vertex_lat_lon_count+1 < len(lons):
				voronoi_next_vertex_lat = lats[voronoi_vertex_lat_lon_count+1]
				voronoi_next_vertex_lon = lons[voronoi_vertex_lat_lon_count+1]
				
			if voronoi_vertex_lat_lon_count+1 >= len(lons):
				
				if voronoi_date_line_intersection_count % 2 == 0: 	
					voronoi_ring_coords.append((lons[0], lats[0]))
					
				if voronoi_date_line_intersection_count % 2 == 1: 	
					voronoi_ring_dateline_intersection_coords.append((lons[0], lats[0]))
				continue
			
			if voronoi_date_line_intersection_count % 2 == 0: 
				voronoi_ring_coords.append((voronoi_vertex_lon, voronoi_vertex_lat))
			
			if voronoi_date_line_intersection_count % 2 == 1: 	
				voronoi_ring_dateline_intersection_coords.append((voronoi_vertex_lon, voronoi_vertex_lat))
			
			""" Get distance and azimuth between voronoi vertices. """
			
			inverse_vincenty(flattening, major_axis, voronoi_vertex_lat, voronoi_vertex_lon, voronoi_next_vertex_lat, voronoi_next_vertex_lon)
			geodesic_distance_between_voronoi_vertices = geodesic_distance_crater_area
			azimuth_between_voronoi_vertices = direction12
			
			distance_geodesic_polygon_vertices = 0
			geodesic_voronoi_vertices_count = 0
			
			""" Mark scipy voronoi vertex as 'previous vertex' in geodesic edge calculation. This is done to consider Date Line intersections 
			that occur between the voronoi vertex and the first geodesic vertex (less than 15 km away).   """
			
			if geodesic_voronoi_vertices_count == 0:
				previous_voronoi_geodesic_vertex_X = voronoi_vertex_lon
				previous_voronoi_geodesic_vertex_Y = voronoi_vertex_lat
			
			""" Divide distance between vertices into 15 km segments and get coordinates of vertices. Vertices from geodesic distances are added to voronoi_ring. """
			
			geodesic_distance_voronoi_vertices = 15000 
			
			""" When vertices of a polygon are closer than the segment length of the geodesic vertices, the distance between geodesic vertices is reduced to 
			one fifth of the distance between the vertices (to get additional vertex in between). When start and end 
			vertex of a polygon are close to the date line on either hemisphere, this would cause the program to believe that there is a polar intersection 
			(because there is only one date line intersection present due to the missing additional vertices over the second date line intersection) 
			resulting in wrong voronoi polygons. """
			
			if geodesic_distance_crater_area < geodesic_distance_voronoi_vertices:
				geodesic_distance_voronoi_vertices = geodesic_distance_crater_area/5
			
			""" When the current or the next vertex is very close to the date line, make sure that the distance between the geodesic vertices in between 
			is lower than that. If not considered, a date line intersection would be recognized as a polar intersection because the second date line overlap 
			would not occur as planned. If either vertex is <-170 deg lon or >170 deg lon (10 deg lon can be very short close to the poles), check distance to date line and if either distance is lower 
			than the geodesic distance to calculate the vertices in between, modify the geodesic distance. """
			
			voronoi_date_line_correction_lon_1 = 170
			voronoi_date_line_correction_lon_2 = -170
			
			""" Close to the poles, a vertex is considered very close to the date line when it is closer than +- 150 deg lon. """
			
			if voronoi_vertex_lat > 89 or voronoi_vertex_lat < -89:
				voronoi_date_line_correction_lon_1 = 150
				voronoi_date_line_correction_lon_2 = -150			
				
			if voronoi_vertex_lon < voronoi_date_line_correction_lon_2 or voronoi_next_vertex_lon < voronoi_date_line_correction_lon_2 or voronoi_vertex_lon > voronoi_date_line_correction_lon_1 or voronoi_next_vertex_lon > voronoi_date_line_correction_lon_1:
				if voronoi_vertex_lon > voronoi_date_line_correction_lon_1 or voronoi_next_vertex_lon > voronoi_date_line_correction_lon_1:
					shapely_dateline = LineString([(180, 90),(180, 0),(180, -90)])
				if voronoi_vertex_lon < voronoi_date_line_correction_lon_2 or voronoi_next_vertex_lon < voronoi_date_line_correction_lon_2:
					shapely_dateline = LineString([(-180, 90),(-180, 0),(-180, -90)])
				
				""" Find closest point on Date Line and get geodesic distance. """
				
				shapely_voronoi_vertex = Point(voronoi_vertex_lon, voronoi_vertex_lat)
				shapely_next_voronoi_vertex = Point(voronoi_next_vertex_lon, voronoi_next_vertex_lat)
				
				voronoi_vertex_closest_point_on_dateline = shapely_dateline.interpolate(shapely_dateline.project(shapely_voronoi_vertex))
				voronoi_next_vertex_closest_point_on_dateline = shapely_dateline.interpolate(shapely_dateline.project(shapely_next_voronoi_vertex))
				 
				voronoi_vertex_closest_point_on_dateline_ogr_X = voronoi_vertex_closest_point_on_dateline.x
				voronoi_vertex_closest_point_on_dateline_ogr_Y = voronoi_vertex_closest_point_on_dateline.y
				
				voronoi_next_vertex_closest_point_on_dateline_ogr_X = voronoi_next_vertex_closest_point_on_dateline.x
				voronoi_next_vertex_closest_point_on_dateline_ogr_Y = voronoi_next_vertex_closest_point_on_dateline.y
				
				inverse_vincenty(flattening, major_axis, voronoi_vertex_lat, voronoi_vertex_lon, voronoi_vertex_closest_point_on_dateline_ogr_Y, voronoi_vertex_closest_point_on_dateline_ogr_X)		
				distance_voronoi_vertex_closest_point_on_dateline = geodesic_distance_crater_area
				
				inverse_vincenty(flattening, major_axis, voronoi_next_vertex_lat, voronoi_next_vertex_lon, voronoi_next_vertex_closest_point_on_dateline_ogr_Y, voronoi_next_vertex_closest_point_on_dateline_ogr_X)			
				distance_next_voronoi_vertex_closest_point_on_dateline = geodesic_distance_crater_area
				
				""" Check whether the current or the next vertex is closest to the date line. """
				
				minimum_distance_current_and_next_voronoi_vertex_to_dateline = min(distance_voronoi_vertex_closest_point_on_dateline, distance_next_voronoi_vertex_closest_point_on_dateline)
				
				""" Use closest distance to the date line to modify the distance between the geodesic vertices for such polygon sections. """
				
				if minimum_distance_current_and_next_voronoi_vertex_to_dateline < 2 * geodesic_distance_voronoi_vertices:
					geodesic_distance_voronoi_vertices = minimum_distance_current_and_next_voronoi_vertex_to_dateline/5
				
			while geodesic_distance_between_voronoi_vertices > geodesic_distance_voronoi_vertices:
				distance_geodesic_polygon_vertices += geodesic_distance_voronoi_vertices
				geodesic_distance_between_voronoi_vertices -= geodesic_distance_voronoi_vertices
				
				direct_vincenty(flattening, major_axis, [[voronoi_vertex_lon, voronoi_vertex_lat, azimuth_between_voronoi_vertices, 0, 0]], distance_geodesic_polygon_vertices, 1)
				
				voronoi_geodesic_vertex_X = buffer_vertices_list[0][0]
				voronoi_geodesic_vertex_Y = buffer_vertices_list[0][1]
				
				""" Consideration of Date Line intersections. If the voronoi polygon crosses the Date Line, the voronoi polygon is merged from 
				individual polygons voronoi_ring and voronoi_ring_dateline_intersection to avoid Date Line ambiguity. """
				
				""" Polygons that intersect the date line cross the date line twice, polygons that intersect the poles intersect the date line once """
					
				""" detect intersection with date line """
				
				voronoi_dateline_intersection = False
				
				if previous_voronoi_geodesic_vertex_X <= 179.999999 and voronoi_geodesic_vertex_X > 179.999999 or previous_voronoi_geodesic_vertex_X > 179.999999 and voronoi_geodesic_vertex_X <= 179.999999:# or previous_voronoi_start_vertex_X <= 179.999999 and voronoi_geodesic_vertex_X > 179.999999 or previous_voronoi_start_vertex_X > 179.999999 and voronoi_geodesic_vertex_X <= 179.999999:
					voronoi_dateline_intersection = True
					voronoi_dateline_hemisphere = "East"
					
					dateline_voronoi = LineString([(179.999999, 90), (179.999999, -90)])
				
				if previous_voronoi_geodesic_vertex_X >= -179.999999 and voronoi_geodesic_vertex_X < -179.999999 or previous_voronoi_geodesic_vertex_X < -179.999999 and voronoi_geodesic_vertex_X >= -179.999999:# or previous_voronoi_start_vertex_X >= -179.999999 and voronoi_geodesic_vertex_X < -179.999999 or previous_voronoi_start_vertex_X < -179.999999 and voronoi_geodesic_vertex_X >= -179.999999:
					voronoi_dateline_intersection = True
					voronoi_dateline_hemisphere = "West"
					
					dateline_voronoi = LineString([(-179.999999, 90), (-179.999999, -90)])		
				
				if voronoi_dateline_intersection ==  True:
					voronoi_date_line_intersection_count += 1
		
					previous_current_voronoi_vertex_line = LineString([(previous_voronoi_geodesic_vertex_X, previous_voronoi_geodesic_vertex_Y),
					                                                     (voronoi_geodesic_vertex_X, voronoi_geodesic_vertex_Y)])
					
					voronoi_dateline_intersection_point = previous_current_voronoi_vertex_line.intersection(dateline_voronoi)
					
					voronoi_dateline_intersection_point_X = voronoi_dateline_intersection_point.x
					voronoi_dateline_intersection_point_Y = voronoi_dateline_intersection_point.y
			
					""" add date line intersection point to voronoi_ring """
					
					if voronoi_date_line_intersection_count % 2 == 1:
						voronoi_ring_coords.append((voronoi_dateline_intersection_point_X, voronoi_dateline_intersection_point_Y))
					
					""" When crossing to the other hemisphere for the first time, the new polygon longitude must always be -179.999999 if it crosses 
					from west to east and +179.999999 when crossing from east to west. This is the first intersection point with the date line. """
					
					if voronoi_dateline_hemisphere == "East":
						voronoi_dateline_intersection_point_X = -179.999999	
						
					if voronoi_dateline_hemisphere == "West":
						voronoi_dateline_intersection_point_X = +179.999999	
	
					if voronoi_date_line_intersection_count % 2 == 0:
						voronoi_ring_coords.append((voronoi_dateline_intersection_point_X, voronoi_dateline_intersection_point_Y))
					
					""" When crossing back to the other hemisphere, the new polygon longitude must always be +179.999999 if it crosses from west to east and -179.999999 
					when crossing from east to west. This is the second intersection point with the date line. """
					
					if voronoi_date_line_intersection_count % 2 == 0:
						if voronoi_dateline_hemisphere == "East":
							voronoi_dateline_intersection_point_X = 179.999999	
							
						if voronoi_dateline_hemisphere == "West":
							voronoi_dateline_intersection_point_X = -179.999999
					
					""" add intersection point to second voronoi ring, voronoi_ring_dateline_intersection """
					
					voronoi_ring_dateline_intersection_coords.append((voronoi_dateline_intersection_point_X, voronoi_dateline_intersection_point_Y))
					
					if voronoi_date_line_intersection_count % 2 == 1:
						old_voronoi_dateline_intersection_point_X = voronoi_dateline_intersection_point_X
						old_voronoi_dateline_intersection_point_Y = voronoi_dateline_intersection_point_Y
					
					""" when entering back to the other hemisphere, add voronoi_ring_dateline_intersection to the voronoi polygon and create a new 
					voronoi_ring_dateline_intersection polygon in case there is another date line intersection """
					
					if voronoi_date_line_intersection_count % 2 == 0: 
						
						""" add first date line intersection point to close ring """
						
						""" When crossing back to the other hemisphere, the new polygon longitude must always be +179.999999 if it crosses from west to east and -179.999999 
						when crossing from east to west. This is the closing point for the polygon (the first intersection point). """
							
						if voronoi_dateline_hemisphere == "East":
							old_voronoi_dateline_intersection_point_X = 179.999999	
							
						if voronoi_dateline_hemisphere == "West":
							old_voronoi_dateline_intersection_point_X = -179.999999
						
						
						voronoi_ring_dateline_intersection_coords.append((old_voronoi_dateline_intersection_point_X, old_voronoi_dateline_intersection_point_Y))
						
						""" get area of voronoi_ring_dateline_intersection. area of voronoi polygon is summed."""
						
						if len(voronoi_ring_dateline_intersection_coords) >= 3:
							dateline_poly = Polygon(voronoi_ring_dateline_intersection_coords)
							# Transform to LAEA for area calculation
							dateline_coords_transformed = [transformer_to_laea.transform(x, y) for x, y in voronoi_ring_dateline_intersection_coords]
							dateline_poly_transformed = Polygon(dateline_coords_transformed)
							voronoi_ring_dateline_intersection_area = dateline_poly_transformed.area
							voronoi_polygon_area += voronoi_ring_dateline_intersection_area
							
							# Add to polygon parts
							voronoi_polygon_parts.append(voronoi_ring_dateline_intersection_coords[:])
						
						""" Create a new empty voronoi_ring_dateline_intersection in case there is another date line intersection coming up. """
						
						voronoi_ring_dateline_intersection_coords = []
				
				""" remember current vertex """
				
				previous_voronoi_geodesic_vertex_X = voronoi_geodesic_vertex_X
				previous_voronoi_geodesic_vertex_Y = voronoi_geodesic_vertex_Y
				
				""" Correct coordinates so that lon is between -179.999999 and +179.999999 deg, to avoid 'jumps' in polygons. Values outside this range are still used to 
				determine the presence of date line intersections and the calculation of geodesic coordinates (Vincenty) """ 
				
				if voronoi_geodesic_vertex_X < -179.999999:
					voronoi_geodesic_vertex_X += 360
				if voronoi_geodesic_vertex_X > 179.999999:
					voronoi_geodesic_vertex_X -= 360
					
				""" Vertices from geodesic distance calculations are added either to voronoi_ring (no date line intersection) or 
				voronoi_ring_dateline_intersection (date line intersection). """
				
				if voronoi_date_line_intersection_count % 2 == 0: 
					voronoi_ring_coords.append((voronoi_geodesic_vertex_X, voronoi_geodesic_vertex_Y))
				
				if voronoi_date_line_intersection_count % 2 == 1: 	
					voronoi_ring_dateline_intersection_coords.append((voronoi_geodesic_vertex_X, voronoi_geodesic_vertex_Y))
				
				geodesic_voronoi_vertices_count += 1
			
			voronoi_vertex_lat_lon_count += 1
		
		""" Here, the processing of a voronoi polygon is finished. If during processing there was only one intersection with the dateline 
		and no crossing back over the date line, the current voronoi polygon intersects the poles. Here, we draw a new polygon from the voronoi_ring 
		and voronoi_ring_dateline_intersection vertices. We check whether there is a north pole or south pole intersection and add the respective pole 
		and the longitudinally-ordered vertices to the polar intersection voronoi polygon. """
		
		if voronoi_date_line_intersection_count % 2 == 1:
			
			""" get all vertices in voronoi_ring and voronoi_ring_dateline_intersection and sort them according to their longitudes """
			
			voronoi_polar_intersection_vertices = []
			
			for voronoi_polar_vertex_lon, voronoi_polar_vertex_lat in voronoi_ring_dateline_intersection_coords:
				voronoi_polar_intersection_vertices.append([voronoi_polar_vertex_lon, voronoi_polar_vertex_lat])
			
			for voronoi_polar_vertex_lon, voronoi_polar_vertex_lat in voronoi_ring_coords:
				voronoi_polar_intersection_vertices.append([voronoi_polar_vertex_lon, voronoi_polar_vertex_lat])
				
			voronoi_polar_intersection_vertices = sorted(voronoi_polar_intersection_vertices, key=lambda coordinates: coordinates[0]) 
			
			""" find pole that is closest to the impact crater - we assume that this indicates whether this is an intersection with the north pole or the south pole """
			
			inverse_vincenty(flattening, major_axis, center_voronoi_polygon_Y, center_voronoi_polygon_X, 89.99999, 0)			
			distance_crater_north_pole = geodesic_distance_crater_area
	
			inverse_vincenty(flattening, major_axis, center_voronoi_polygon_Y, center_voronoi_polygon_X, -89.99999, 0)			
			distance_crater_south_pole = geodesic_distance_crater_area
			
			if distance_crater_north_pole <= distance_crater_south_pole:
				voronoi_polar_vertex_Y = 89.99999
				
			if distance_crater_north_pole > distance_crater_south_pole:
				voronoi_polar_vertex_Y = -89.99999
			
			""" Draw a new polygon (voronoi_ring) for polar intersections. The poles and date lines are not exactly at +- 179.999999 deg lon or +-90 deg lat. 
			Date lines are sometimes automatically converted from -179.999999 to +179.999999 degrees in OGR which is not very helpful in this case. 
			Polar intersections may result in gaps if the latitude is set to +- 90 deg. """
			
			voronoi_ring_coords = []
			
			voronoi_ring_coords.append((-179.99999, voronoi_polar_vertex_Y)) # closest pole
			
			for voronoi_polar_intersection_vertex in voronoi_polar_intersection_vertices:
				voronoi_polar_intersection_vertex_X = voronoi_polar_intersection_vertex[0]
				voronoi_polar_intersection_vertex_Y = voronoi_polar_intersection_vertex[1]
				
				if voronoi_polar_intersection_vertex_X == -179.999999:
					voronoi_polar_intersection_vertex_X = -179.99999
					
				if voronoi_polar_intersection_vertex_X == 179.999999:
					voronoi_polar_intersection_vertex_X = 179.99999
					
				voronoi_ring_coords.append((voronoi_polar_intersection_vertex_X, voronoi_polar_intersection_vertex_Y))
				
			voronoi_ring_coords.append((179.99999, voronoi_polar_vertex_Y)) # closest pole	
	
		voronoi_count += 1	
		
		print(voronoi_count)
		print("---")	
		
		""" Get area of voronoi_ring. The total area is eventually summed from voronoi_ring and voronoi_ring_dateline_intersection due to 
		inaccurate area measurements when one polygon intersects the date line intersections. """
		
		if len(voronoi_ring_coords) >= 3:
			# Transform to LAEA for area calculation
			ring_coords_transformed = [transformer_to_laea.transform(x, y) for x, y in voronoi_ring_coords]
			ring_poly_transformed = Polygon(ring_coords_transformed)
			voronoi_ring_area = ring_poly_transformed.area
			voronoi_polygon_area += voronoi_ring_area
			
			# Add to polygon parts
			voronoi_polygon_parts.append(voronoi_ring_coords[:])
		
		voronoi_polygon_areas.append(voronoi_polygon_area)
		
		""" Create GeoJSON feature from polygon parts """
		
		if len(voronoi_polygon_parts) == 1:
			# Simple polygon
			geom = Polygon(voronoi_polygon_parts[0])
		elif len(voronoi_polygon_parts) > 1:
			# MultiPolygon
			polys = [Polygon(part) for part in voronoi_polygon_parts]
			geom = MultiPolygon(polys)
		else:
			# Fallback to empty geometry
			geom = Polygon()
		
		polygon_features.append({
			"type": "Feature",
			"geometry": mapping(geom),
			"properties": {
				"Area": voronoi_polygon_area
			}
		})
	
	# Write GeoJSON files
	points_geojson = {
		"type": "FeatureCollection",
		"features": points_features
	}
	
	polygons_geojson = {
		"type": "FeatureCollection",
		"features": polygon_features
	}
	
	# Replace .shp extension with .geojson for output files
	points_output_path = input_samples_point_shapefile_output_path.replace('.shp', '.geojson')
	polygons_output_path = voronoi_polygons_shapefile_output_path.replace('.shp', '.geojson')
	
	with open(points_output_path, 'w') as f:
		json.dump(points_geojson, f, indent=2)
	
	with open(polygons_output_path, 'w') as f:
		json.dump(polygons_geojson, f, indent=2)
	
	print("Done.")
	print(f"Points saved to: {points_output_path}")
	print(f"Polygons saved to: {polygons_output_path}")

""" Change parameters here: """

# True - if randomly distributed points should be used / False - if individual set of points should be used (modify input_samples below)
use_random_set_of_input_points = True 

# if random points should be used - 20 points minimum
number_of_random_points = 999 

# List of points to generate Voronoi polygons from. Minimum of 20 points required. Must set use_random_set_of_input_points = False
input_samples = [[0.1,0.1],[180,0.1], [90,0], [-90,0], [0.1,89], [0.1,-89], [20,20], [-20,-20], [-40,-40], [40, 40], [-60,-60], [60, 60], [0, 80], [-100, 80], [-10, 10], [-30, -30], [-25, 50], [40, -40], [160, -50], [-70, 60]]

# Output paths - use raw strings (r'...') or forward slashes to avoid escape sequence issues
# Example: r'C:\path\to\point_file.shp' or 'C:/path/to/point_file.shp'
input_samples_point_shapefile_output_path = r'D:\Download\point_file.shp' 

# Example: r'C:\path\to\polygon_file.shp' or 'C:/path/to/polygon_file.shp'
voronoi_polygons_shapefile_output_path = r'D:\Download\polygon_file.shp' 

# The reference body is defined from a geographic coordinate system WKT. In this example, a sphere with the Earth's dimensions is used. 
geogr_sr_text = 'GEOGCS["WGS 84", DATUM["WGS_1984", SPHEROID["WGS 84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["Decimal_Degree",0.0174532925199433]]'

""" Start main function """

main(use_random_set_of_input_points, number_of_random_points, input_samples, input_samples_point_shapefile_output_path, voronoi_polygons_shapefile_output_path, geogr_sr_text)	
