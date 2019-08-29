## Finding flowline networks for all ~200 MEaSUREs outlet glaciers
## 24 Sept 2018  EHU
## 10 Jan 19 - edit to use fuller MEaSUREs velocity composite, catching last glaciers stored empty
import numpy as np
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter #use Savitzky-Golay filter to smooth catchments in Trace_wWidth
from osgeo import gdal
import sys #allowing GDAL to throw Python exceptions
import pandas as pd #want to treat velocity maps as Pandas dataframes
import shapefile
import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
#from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from network_selection import *


##Reading in velocities -- function lifted from Greenland-vel-compositing.py
##Function to read MEaSUREs velocity GeoTIFFs
def read_velocities(filename, return_grid=True, return_proj=False):
    """Extract x, y, v from a MEaSUREs GeoTIFF"""
    ds = gdal.Open(filename)
    #Get dimensions
    nc = ds.RasterXSize
    nr = ds.RasterYSize
    
    geotransform = ds.GetGeoTransform()
    xOrigin = geotransform[0]
    xPix = geotransform[1] #pixel width in x-direction
    yOrigin = geotransform[3]
    yPix = geotransform[5] #pixel height in y-direction
    
    lons = xOrigin + np.arange(0, nc)*xPix
    lats = yOrigin + np.arange(0, nr)*yPix
    
    x, y = np.meshgrid(lons, lats)
    
    vband = ds.GetRasterBand(1)
    varr = vband.ReadAsArray()
    
    #if return_grid and return_proj:
    #    return x, y, varr, ds.GetProjection()
    #elif return_grid:
    if return_grid:
        return x, y, varr
    else: 
        return varr



## Read in MEaSUREs velocity composite
print 'Reading MEaSUREs velocities'
x_comp, y_comp, v_comp_raw = read_velocities('Documents/GitHub/Data_unsynced/gld-velocity-composite-10Jan19.tif')
vx_comp_raw = read_velocities('Documents/GitHub/Data_unsynced/gld-x_velocity-composite-10Jan19.tif', return_grid=False)
vy_comp_raw = read_velocities('Documents/GitHub/Data_unsynced/gld-y_velocity-composite-10Jan19.tif', return_grid=False)
v_comp = np.ma.masked_invalid(v_comp_raw)
vx_comp = np.ma.masked_invalid(vx_comp_raw)
vy_comp = np.ma.masked_invalid(vy_comp_raw)
v_excludemasked = np.ma.filled(v_comp, fill_value=0)
vx_excludemasked = np.ma.filled(vx_comp, fill_value=0)
vy_excludemasked = np.ma.filled(vy_comp, fill_value=0)
days_per_annum = 365.242 #need to convert units of MEaSUREs velocity to align with what we used from Sentinel before
v_a2d = np.array(v_excludemasked) / days_per_annum
vx_a2d = np.array(vx_excludemasked) / days_per_annum
vy_a2d = np.array(vy_excludemasked) / days_per_annum


##Make 2D-interpolated function of velocity field for tracing
print 'Interpolating MEaSUREs velocity composites for tracing'
x_flat = x_comp[0,:]
y_flat = y_comp[:,0]
func_vxcomp = interpolate.interp2d(x_flat, y_flat[::-1], vx_a2d) #these need to be flipped along y-axis
func_vycomp = interpolate.interp2d(x_flat, y_flat[::-1], vy_a2d)
func_vcomp = interpolate.interp2d(x_flat, y_flat[::-1], v_a2d)

#xtest = np.linspace(min(x_flat), max(x_flat), 1000)
#ytest = np.linspace(min(y_flat), max(y_flat), 500)
#plt.figure()
#plt.contour(xtest, ytest, func_vcomp(xtest, ytest))
#plt.show()

## Using read_termini from MEaSUREs-validation.py
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
terminus_basefile = '/termini_0607_v01_2'
init_year = 2006
fn = gl_termpos_fldr + terminus_basefile #filename to read in for termini that will be traced
termini_init = read_termini(fn, init_year)

##iterate over keys in termini_init to make dictionary of lines for each GlacierID
#ids_to_trace = termini_init.keys() #trace all points of all glaciers
#ids_to_trace = (3, 153, 175) # IDs for only Jakobshavn, Kangerlussuaq, Helheim
#ids_to_trace = (143, 144, 145, 146, 147, 148, 149, 151, 159, 177, 195) #branches that previously returned a KeyError when filtered (4 Oct)
ids_stored_empty = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177) #glaciers whose termini lie off previous velocity composite (10 Jan 19)
#ids_to_trace = range(201, 230)

all_lines = {}
for gid in ids_stored_empty:
    lines = {}
    termcoords = termini_init[gid] #points spanning terminus for this glacier
    sampling = int(np.floor(len(termcoords)/10))+1 #allowing up to 10 branches traced, and with +1 ensuring that we don't step by 0
    for j in range(len(termcoords))[::sampling]: #need to down-sample, since some glaciers have 50 points across terminus (and we expect most will have fairly simple networks)
        print 'Tracing terminus point {} of {} in Glacier ID {}'.format(j, len(termcoords), gid)
        p = termcoords[j]
        line_coords, width = Trace_wWidth(p[0], p[1], trace_up=True, xarr=x_comp, yarr=y_comp, Vx = func_vxcomp, Vy = func_vycomp, V = func_vcomp) #Uses Trace_wWidth and FilterMainTributaries from network_selection.py
        if np.any(np.isnan(line_coords)):
            pass #Don't save the line if it contains nan points
        else:
            xyw = [(line_coords[n][0], line_coords[n][1], width[n]) for n in range(len(line_coords))]
            lines[j] = (xyw)
    filtered_tribs = FilterMainTributaries(lines, Vx = func_vxcomp, Vy = func_vycomp) 
    all_lines[gid] = filtered_tribs #saving this to the dictionary of all lines we're working with
    
    #writing to CSV (based on Write_Network)
    out_dir = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/'
    output_name = out_dir + 'Gld-autonetwork-GID{}-date_{}.csv'.format(gid, datetime.date.today())
    with open(output_name, 'wb') as csvfile:
        linewriter = csv.writer(csvfile, delimiter=',')
        linewriter.writerow(['Line-number', 'x-coord', 'y-coord', 'width'])
        for n in range(len(filtered_tribs)):
            for m in range(len(filtered_tribs[n])):
                linewriter.writerow([str(n), filtered_tribs[n][m][0], filtered_tribs[n][m][1], filtered_tribs[n][m][2]])


#plt.figure()
#for k in filtered_tribs.keys():
#    trib = np.asarray(filtered_tribs[k])
#    #print shape(trib), type(trib)
#    plt.plot(trib[:,0], trib[:,1], label='Line {}'.format(k))
#plt.legend()
#plt.show()

#plt.figure()
#for j in range(len(termcoords))[::10]:
#    p=termcoords[j]
#    plt.plot(p[0], p[1], label='Point {}'.format(j), marker='*', ms=20)
#plt.legend()
#plt.show()