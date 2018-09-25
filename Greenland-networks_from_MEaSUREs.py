## Finding flowline networks for all ~200 MEaSUREs outlet glaciers
## 24 Sept 2018  EHU
import numpy as np
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter #use Savitzky-Golay filter to smooth catchments in Trace_wWidth
from osgeo import gdal
import sys #allowing GDAL to throw Python exceptions
import pandas as pd #want to treat velocity maps as Pandas dataframes
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
#from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter


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
x_comp, y_comp, v_comp_raw = read_velocities('Documents/GitHub/gld-velocity-composite.tif')
vx_comp_raw = read_velocities('Documents/GitHub/gld-x_velocity-composite.tif', return_grid=False)
vy_comp_raw = read_velocities('Documents/GitHub/gld-y_velocity-composite.tif', return_grid=False)
v_comp = np.ma.masked_invalid(v_comp_raw)
vx_comp = np.ma.masked_invalid(vx_comp_raw)
vy_comp = np.ma.masked_invalid(vy_comp_raw)
v_excludemasked = np.ma.filled(v_comp, fill_value=0)
vx_excludemasked = np.ma.filled(vx_comp, fill_value=0)
vy_excludemasked = np.ma.filled(vy_comp, fill_value=0)


##Make 2D-interpolated function of velocity field for tracing
print 'Interpolating MEaSUREs velocity composites for tracing'
#xm = x_comp[~v_comp.mask]
#ym = y_comp[~v_comp.mask]
x_flat = x_comp[0,:]
y_flat = y_comp[:,0]
func_vxcomp = interpolate.interp2d(x_flat, y_flat[::-1], vx_excludemasked) #check whether these need to be flipped along y-axis
func_vycomp = interpolate.interp2d(x_flat, y_flat[::-1], vy_excludemasked)
func_vcomp = interpolate.interp2d(x_flat, y_flat[::-1], v_excludemasked)

xtest = np.linspace(min(x_flat), max(x_flat), 1000)
ytest = np.linspace(min(y_flat), max(y_flat), 500)
plt.figure()
plt.contour(xtest, ytest, func_vcomp(xtest, ytest))
plt.show()

### Using read_termini from MEaSUREs-validation.py
#gl_termpos_fldr = 'Documents/GitHub/plastic-networks/Data/MEaSUREs-termini'
#terminus_basefile = '/termini_0607_v01_2'
#init_year = 2006
#fn = gl_termpos_fldr + terminus_basefile #filename to read in for termini that will be traced
#termini_init = read_termini(fn, init_year)
#
###iterate over keys in termini_init to make dictionary of lines for each GlacierID
##ids_to_trace = termini_init.keys() #trace all points of all glaciers
#ids_to_trace = (3, 153, 175) # IDs for only Jakobshavn, Kangerlussuaq, Helheim
#
#all_lines = {}
#for gid in ids_to_trace:
#    lines = {}
#    termcoords = termini_init[gid] #points spanning terminus for this glacier
#    for j in range(len(termcoords)):
#        p = termcoords[j]
#        line_coords, width = Trace_wWidth(p[0], p[1], trace_up=True, xarr=x_comp, yarr=y_comp, Vx = func_vxcomp, Vy = func_vycomp, V = func_vcomp) #Uses Trace_wWidth and FilterMainTributaries from network_selection.py
#        xyw = [(line_coords[n][0], line_coords[n][1], width[n]) for n in range(len(line_coords))]
#        lines[j] = (xyw)
#    filtered_tribs = FilterMainTributaries(lines, Vx = func_vxcomp, Vy = func_vycomp)
#    all_lines[gid] = lines

        