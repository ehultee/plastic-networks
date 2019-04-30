## Nov 2018: Plotting map of MEaSUREs outlets over satellite imagery
## Apr 2019 edit: Visualise hindcast validation, e.g. total retreat and misfit with observation
## EHU

from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np
import matplotlib.pyplot as plt
import shapefile
 
##--------------------------
## SET UP MAP
##--------------------------

#selection of which map appearance you want, projection settings
services = ['World_Physical_Map', 'World_Shaded_Relief', 'World_Topo_Map', 'NatGeo_World_Map', 'ESRI_Imagery_World_2D', 'World_Street_Map', 'World_Imagery', 'Ocean_Basemap']
wgs84 = pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
psn_gl = pyproj.Proj("+init=epsg:3413") # Polar Stereographic North used by BedMachine 



def Greenland_map(service='ESRI_Imagery_World_2D', epsg=3413, xpixels=5000):
    """Function using Basemap to plot map for all of Greenland
    """
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=300, llcrnrlat=57, urcrnrlon=20, urcrnrlat=80, resolution='h')
    
    plt.figure()
    m.arcgisimage(service=service, xpixels=xpixels)
    plt.show()
    return m

##Convert coords into lat/lon so that Basemap can convert them back (don't know why this is necessary, but it works)
def flowline_latlon(coords, fromproj=pyproj.Proj("+init=epsg:3413"), toproj=pyproj.Proj("+init=EPSG:4326")):
    """Convert coords into lat/lon so that Basemap can convert them back for plotting (don't know why this is necessary, but it works)
    Defaults:
        fromproj = NSIDC Polar Stereographic North
        toproj = WGS84 lat-lon
    """
    xs = coords[:,0]
    ys = coords[:,1]
    x_lon, y_lat = pyproj.transform(fromproj, toproj, xs, ys)
    latlon_coords = np.asarray(zip(x_lon, y_lat))
    return latlon_coords
    
##--------------------------
## SET UP OUTLET MARKERS
##--------------------------

def read_termini(filename, year):
    """Make a dictionary of terminus positions, indexed by MEaSUREs ID. Outputs dictionary"""
    print 'Reading in MEaSUREs terminus positions for year ' + str(year)
    sf = shapefile.Reader(filename)
    fields = sf.fields[1:] #excluding the mute "DeletionFlag"
    field_names = [field[0] for field in fields]
    term_recs = sf.shapeRecords()
    termpts_dict = {}
    for r in term_recs:
        atr = dict(zip(field_names, r.record)) #dictionary of shapefile fields, so we can access GlacierID by name rather than index.  Index changes in later years.
        key = atr['GlacierID'] #MEaSUREs ID number for the glacier, found by name rather than index
        termpts_dict[key] = np.asarray(r.shape.points) #save points spanning terminus to dictionary
    return termpts_dict

gl_termpos_fldr = 'Documents/GitHub/Data_unsynced/MEaSUREs-termini'
terminus_basefile = '/termini_1415_v01_2'
init_year = 2014
fn = gl_termpos_fldr + terminus_basefile #filename to read in for termini that will be traced
termini_init = read_termini(fn, init_year)

all_coords_latlon = {}
for gid in termini_init.keys():
    latlon_coords = flowline_latlon(termini_init[gid])
    all_coords_latlon[gid] = np.asarray(latlon_coords)

##--------------------------
## MAKE PLOT
##--------------------------

gld_backdrop = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=3413, llcrnrlon=300, llcrnrlat=57, urcrnrlon=20, urcrnrlat=80, resolution='h')

#plt.figure()
#gld_backdrop.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=5000)
##for k in termini_init.keys():
##    pt = all_coords_latlon[k][0]
##    gld_backdrop.scatter(pt[0], pt[1], s=60, marker='*', color='LightCyan', edgecolors='b', latlon=True)
##    term_marker = gld_backdrop(pt[0], pt[1])
##    #offset = 100 * np.mod(k,2)
##    plt.annotate(s=str(k), xy=term_marker, fontsize='small', fontweight='demi', color='MediumBlue')
#plt.show()

plt.figure()
gld_backdrop.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=5000)
for k in termini_init.keys():
    pt = all_coords_latlon[k][0]
    size = k #testing variable size
    gld_backdrop.scatter(pt[0], pt[1], s=size, marker='*', color='b', edgecolors='DarkViolet', latlon=True)
    #    term_marker = gld_backdrop(pt[0], pt[1])
    #    #offset = 100 * np.mod(k,2)
    #    plt.annotate(s=str(k), xy=term_marker, fontsize='small', fontweight='demi', color='MediumBlue')
plt.show()