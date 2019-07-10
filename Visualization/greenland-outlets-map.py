## Nov 2018: Plotting map of MEaSUREs outlets over satellite imagery
## Apr 2019 edit: Visualise hindcast validation, e.g. total retreat and misfit with observation
## EHU

from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np
import matplotlib.pyplot as plt
import shapefile
from GL_model_tools import Greenland_map, flowline_latlon, read_termini
 

##--------------------------
## SET UP OUTLET MARKERS
##--------------------------

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

plt.figure()
gld_backdrop.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=5000)
for k in termini_init.keys():
    pt = all_coords_latlon[k][0]
    gld_backdrop.scatter(pt[0], pt[1], s=40, marker='o', color='Blue', edgecolors='DarkViolet', latlon=True)
    #term_marker = gld_backdrop(pt[0], pt[1])
    #offset = 100 * np.mod(k,2)
    #plt.annotate(s=str(k), xy=term_marker, fontsize='small', fontweight='demi', color='MediumBlue')
# Now plot the case study glaciers on top
for k in (3, 109, 137):
    pt = all_coords_latlon[k][0]
    gld_backdrop.scatter(pt[0], pt[1], s=180, marker='*', color='Yellow', edgecolors='Gold', latlon=True)
plt.show()