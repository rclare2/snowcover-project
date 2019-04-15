'''
region_masker.py

Makes binary masks of various North American regions by multiplying crude
latitude/longitude outlines with a landmask file.  Resolution is 30 km and 
works with the WRF domain used for the Desai Lab Snow Cover Project data. 

These masks are used primarily to measure changes in regional precipitation
but may be used by multiplying their arrays with any 2D, unstaggered variable. 

Will need access to the basic_stats library, authored by Ryan Clare.

Ryan Clare
April, 2019
'''

from numpy import *
from basic_stats import *

dom = [209,143]

# Read WRF LANDMASK

f = open("/air/rclare/wrf_grid_data/LANDMASK.txt", "r")
content = map(float, f)
f.close()

LANDMASK = asarray(content).reshape(dom[1], dom[0])

LANDMASK[LANDMASK > 1.0] = 0.0


def ll_reader():
	# Parse latitude and longitude data from lat/lon index file
	
	f = open("/air/rclare/wrf_grid_data/wrf_ll.txt", "r")
	raw = f.readlines()
	f.close()
	
	longitudes, latitudes, grid_id = [], [], []
	
	for line in raw:
		mid = line.find("|")
		grid_nums = line[:mid-1]
		latlons = line[mid+2:]
		ll_mid = latlons.find(",")
		lons = float(latlons[:ll_mid])
		lats = float(latlons[ll_mid+2:])	
		longitudes.append(lons)
		latitudes.append(lats)
		grid_id.append(ast.literal_eval(grid_nums))
		
	return longitudes, latitudes, grid_id
	

def LL2grid(longitudes, latitudes):
	# Converts longitude-latitude list to wrf-grid coordinates
	
	lon, lat, grd = ll_reader()
	path_x, path_y = [], []
	xs, ys = [], []
	
	# Scan points to see which yours are closes to
	for i in range(len(longitudes)):
		min_dist = 1e9					# Set impossibly high minimum distance
		for j in range(len(grd)):
			r = pythagorean(latitudes[i], longitudes[i], lat[j], lon[j])
			if r < min_dist:			# If current distance is less than the minimum:
				nearest_pt = grd[j]		# Save current grid point
				min_dist = r			# Set new minimum distance as current
		
		xs.append(nearest_pt[0])
		ys.append(nearest_pt[1])
			
	return xs, ys
	

def wrfGrid2ll(xs, ys):
	# Converts wrf-grid list to longitude-latitude coordinates
	
	lon, lat, grd = ll_reader()
	path_x, path_y = [], []
	
	for i in range(len(xs)):
		for j in range(len(grd)):
			if [xs[i], ys[i]] == list(grd[j]):
				path_x.append(lon[j])
				path_y.append(lat[j])
				break
			
	return path_x, path_y


def closer(xs, ys):
	# If lat/lon specifications do not form a closed polygon, close it
	
	if xs[-1] != xs[0] and ys[-1] != ys[0]:
		
		xs.append(xs[-1])
		ys.append(dom[1]-4)
		
		xs.append(xs[0])
		ys.append(dom[1]-4)
		
		xs.append(xs[0])
		ys.append(ys[0])
		
	return xs, ys


def form_polygon(xs, ys):
	# Fill lines between polygon vectors, forming a closed polygon
	
	all_xs, all_ys = [], []
	
	for i in range(len(xs)-1):
		
		segment_xs = (xs[i], xs[i+1])
		segment_ys = (ys[i], ys[i+1])
		segment_length = segment_ys[1] - segment_ys[0]
		
		all_xs.append(segment_xs[0])
		all_ys.append(segment_ys[0])
		
		if segment_xs[0] != segment_xs[1]:
			# Calculate segment slope
			m = float((segment_ys[1] - segment_ys[0])/float((segment_xs[1] - segment_xs[0])))
			
			# Calculate x all the way up
			if abs(segment_length) > 1:
				
				for j in range(abs(segment_length)-1):
					
					if segment_length > 0:
						y = segment_ys[0] + j
					elif segment_length < 0:
						y = segment_ys[0] - j
						
					x = (y - segment_ys[0])/m + segment_xs[0]
					
					all_xs.append(int(round(x)))
					all_ys.append(y)
		
		else:
			if abs(segment_length) > 1:
				
				for j in range(abs(segment_length)-1):
					
					if segment_length > 0:
						y = segment_ys[0] + j
					elif segment_length < 0:
						y = segment_ys[0] - j
					
					all_xs.append(segment_xs[0])
					all_ys.append(y)
			
		all_xs.append(segment_xs[1])
		all_ys.append(segment_ys[1])
		
	return all_xs, all_ys


def ray_casting(xs, ys):
	# Find which points are within the polygon (not perfect...)
	# https://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
	
	region = zeros((dom[1], dom[0]))
	
	polygon_walls = []
	for i in range(len(xs)):
		polygon_walls.append([xs[i],ys[i]])
		
	x_span, y_span = range(dom[0]), range(dom[1])
	
	for i in x_span:
		
		for j in y_span:
			
			intersects = 0
			for point in x_span[i:]:
				
				POI = [point, j]
				for k in polygon_walls:
					if k == POI:
						intersects+=1
				
			if intersects%2 == 1:
				region[j,i] = 1.
	
	# Patchwork fix (erasure)
	 # Could be fixed by eliminating vertices from sides BUT
	  # then you'd also need to deal with parallel lines
	for j in y_span:
		
		if region[j,0] == 1.:
			
			for i in x_span:
				
				wall = False
				for k in polygon_walls:
					if k == [i, j]:
						wall = True
						break
						
				if wall == True:
					break
				else:
					region[j,i] = 0.
	
	return region
	
	
def write_REGIONMASK(region, name):
	# Write region mask to file (after landmask has been applied)
	
	f = open("region_masks/%s.txt" % name, "w")
	
	x_span, y_span = range(dom[0]), range(dom[1])
	
	for j in y_span:
		
		for i in x_span:
		
			f.write(str(region[j,i]) + "\n")
			
	f.close()
	
	print name


###########################################################
#                                                         #
#                     Canadian Plains                     #
#                                                         #
###########################################################

lats = [  60.003,  60.008,  53.702,  49.000,  49.000,  49.000, 49.000, 52.840, 56.915]

lons = [-123.921,-120.008,-120.008,-114.068,-110.005,-101.363,-95.154,-95.154,-88.893] 

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = closer(Xs, Ys)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "CanPlains")


###########################################################
#                                                         #
#                         Ontario                         #
#                                                         #
###########################################################

lats = [ 49.000, 52.840, 56.915, 54.934, 47.110, 45.456, 45.411, 43.648, 43.456, 42.839, 41.822, 43.594, 45.354, 48.305, 48.034, 49.000]

lons = [-95.154,-95.154,-88.893,-79.880,-79.501,-75.721,-73.988,-77.124,-79.191,-78.935,-83.183,-82.122,-82.517,-88.376,-89.600,-95.154]

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "Ontario")


###########################################################
#                                                         #
#                    Quebec and Labrador                  #
#                                                         #
###########################################################

lats = [ 58.556, 54.934, 47.110, 45.456, 45.411, 45.008, 45.022, 47.537, 48.662, 50.288, 50.114, 52.638]

lons = [-83.053,-79.880,-79.501,-75.721,-73.988,-74.669,-71.526,-69.293,-63.684,-65.379,-59.816,-54.307]

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = closer(Xs, Ys)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "Quebec")


###########################################################
#                                                         #
#                     Maritime Canada                     #
#                                                         #
###########################################################

lats = [ 47.537, 48.662, 50.288, 50.114, 52.638, 49.372, 41.049, 43.191, 45.725, 47.359, 47.537]
                                                                                        
lons = [-69.293,-63.684,-65.379,-59.816,-54.307,-51.718,-64.648,-65.991,-67.716,-67.919,-69.293]

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "Maritime")


###########################################################
#                                                         #
#                        The South                        #
#                                                         #
###########################################################

lats = [  32.013,  32.006,  37.004,  37.001,  40.003, 40.003, 38.977, 36.508, 36.481, 34.995, 35.010, 28.351, 25.609,  29.445,  32.013]
                                                                                                                              
lons = [-106.718,-103.064,-103.057,-102.050,-102.062,-95.314,-94.599,-94.615,-89.495,-90.327,-88.096,-88.553,-97.171,-104.306,-106.718] 

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "South")


###########################################################
#                                                         #
#                      The Southeast                      #
#                                                         #
###########################################################

lats = [ 35.010, 28.351, 23.332, 30.203, 37.897, 37.936, 38.423, 38.957, 39.690, 37.693, 36.598, 36.598, 35.009, 35.010]
                                                                                                                       
lons = [-88.096,-88.553,-79.077,-80.038,-70.542,-76.240,-76.981,-76.931,-78.113,-80.301,-83.696,-81.685,-84.354,-88.096]

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "Southeast")


###########################################################
#                                                         #
#                       Central U.S.                      #
#                                                         #
###########################################################

lats = [ 39.690, 37.693, 36.598, 36.598, 35.009, 34.995, 36.481, 36.508, 38.977, 40.586, 40.595, 42.518, 42.519, 41.784, 41.822, 42.128, 39.721, 39.690]
                                                                                                                                                       
lons = [-78.113,-80.301,-83.696,-81.685,-84.354,-90.327,-89.495,-94.615,-94.599,-95.817,-91.494,-90.590,-87.665,-87.047,-83.183,-80.520,-80.508,-78.113] 

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "CentralUS")


###########################################################
#                                                         #
#                      Upper Midwest                      #
#                                                         #
###########################################################

lats = [ 40.586, 40.595, 42.518, 42.519, 41.784, 41.822, 43.594, 45.354, 48.305, 48.034, 49.000, 49.000, 42.682, 40.586]
                                                                                                                       
lons = [-95.817,-91.494,-90.590,-87.665,-87.047,-83.183,-82.122,-82.517,-88.376,-89.600,-95.154,-97.233,-96.590,-95.817]

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "Midwest")


###########################################################
#                                                         #
#                      Northeast U.S.                     #
#                                                         #
###########################################################

lats = [ 37.897, 37.936, 38.423, 38.957, 39.690, 39.721, 42.128, 42.839, 43.456, 43.648, 45.411, 45.008, 45.022, 47.537, 47.359, 45.725, 43.191, 37.897]
                                                                                                                                                       
lons = [-70.542,-76.240,-76.981,-76.931,-78.113,-80.508,-80.520,-78.935,-79.191,-77.124,-73.988,-74.669,-71.526,-69.293,-67.919,-67.716,-65.991,-70.542]

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "Northeast")


###########################################################
#                                                         #
#                       West Central                      #
#                                                         #
###########################################################

lats = [  49.000,  37.013,  36.997,  40.003, 40.018, 42.682, 49.000,  49.000,  49.000,  49.000]

lons = [-117.034,-109.038,-102.043,-102.052,-95.384,-96.590,-97.233,-101.363,-114.068,-117.034]

Xs, Ys = LL2grid(lons, lats)

Xs, Ys = form_polygon(Xs, Ys)

REGION = ray_casting(Xs, Ys) * LANDMASK

write_REGIONMASK(REGION, "MtnWest")



'''
# Tool for plotting outlines on map

from matplotlib.pyplot import *
from mpl_toolkits.basemap import Basemap

figure()

m = Basemap(llcrnrlon=-126.5696,llcrnrlat=19.3052,urcrnrlon=-46.00888,urcrnrlat=53.53284,\
				rsphere=(6378137.00,6356752.3142),resolution='l',projection='lcc',\
				lat_0=43.5,lon_0=-98.,lat_1=30.,lat_2=60.)
				
m.drawcoastlines()
m.drawcountries()
m.drawstates()

path_x, path_y = m(lons,lats)

plot(path_x, path_y, 'r-')

show()
'''
