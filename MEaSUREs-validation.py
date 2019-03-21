## Loading in MEaSUREs terminus position data for Greenland to assess utility for validating hindcasts
## 12 Sept 2018  EHU
## 20 Mar 2019 edit: use these functions to compare 2006-2014 termini with hindcasted

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
import shapefile
#from matplotlib.colors import LogNorm
from matplotlib import cm
#from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import shapely.geometry as geom
from plastic_utilities_v2 import *
from GL_model_tools import *
from flowline_class_hierarchy import *

##-------------------
### READING IN BED, VELOCITY, SURFACE CHANGE ETC.
### COMMENT OUT IF DATA IS ALREADY READ IN TO YOUR SESSION
##-------------------

#print 'Reading in surface topography'
#gl_bed_path ='Documents/1. Research/2. Flowline networks/Model/Data/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
#fh = Dataset(gl_bed_path, mode='r')
#xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
#yy = fh.variables['y'][:].copy() #y-coord
#s_raw = fh.variables['surface'][:].copy() #surface elevation
#h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
#b_raw = fh.variables['bed'][:].copy() # bed topo
#thick_mask = fh.variables['mask'][:].copy()
#ss = np.ma.masked_where(thick_mask !=2, s_raw)#mask values: 0=ocean, 1=ice-free land, 2=grounded ice, 3=floating ice, 4=non-Greenland land
#hh = np.ma.masked_where(thick_mask !=2, h_raw) 
#bb = np.ma.masked_where(thick_mask !=2, b_raw)
### Down-sampling
#X = xx[::2]
#Y = yy[::2]
#S = ss[::2, ::2]
#H = hh[::2, ::2] 
#B = bb[::2, ::2]
### Not down-sampling
##X = xx
##Y = yy
##S = ss
#fh.close()
#
##Smoothing bed to check effect on dLdt
#unsmoothB = B
#smoothB = gaussian_filter(B, 2)
##B_processed = np.ma.masked_where(thick_mask !=2, smoothB)
#
#S_interp = interpolate.RectBivariateSpline(X, Y[::-1], S.T[::, ::-1])
#H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
#B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])

### Reading in SENTINEL velocity map
#print 'Now reading in (vector) velocity map'
#v_path = 'Documents/1. Research/2. Flowline networks/Model/Data/ESA-Greenland/greenland_iv_500m_s1_20161223_20170227_v1_0.nc'
#fh2 = Dataset(v_path, mode='r')
#xv = fh2.variables['x'][:].copy()
#yv = fh2.variables['y'][:].copy()
##yv = yv_flipped[::-1]
#v_raw = fh2.variables['land_ice_surface_velocity_magnitude'][:].copy() #this is v(y, x)
#vx_raw = fh2.variables['land_ice_surface_easting_velocity'][:].copy()
#vy_raw =fh2.variables['land_ice_surface_northing_velocity'][:].copy()
#v_upper = np.ma.masked_greater(v_raw, 10000)
#vx_upper = np.ma.masked_greater(vx_raw, 10000)
#vy_upper = np.ma.masked_greater(vy_raw, 10000)
#fh2.close()
#
### Interpolate SENTINEL and sample at BedMachine points
#print 'Now interpolating to same grid'
#vf_x = interpolate.interp2d(xv, yv[::-1], vx_upper[::-1,::])
#vf_y = interpolate.interp2d(xv, yv[::-1], vy_upper[::-1,::])
#vf = interpolate.interp2d(xv, yv[::-1], v_upper[::-1, ::])
#

## Read in MEaSUREs velocity composite
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


print 'Reading in MEaSUREs reference file' 
gl_gid_fldr = 'Documents/GitHub/Data_unsynced/MEaSUREs-GlacierIDs'
sf_ref = shapefile.Reader(gl_gid_fldr+'/GlacierIDs_v01_2') #Specify the base filename of the group of files that makes up a shapefile


## Reading terminus positions consistently
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
basefiles = ['/termini_0001_v01_2', '/termini_0506_v01_2', '/termini_0607_v01_2', '/termini_0708_v01_2', '/termini_0809_v01_2', '/termini_1213_v01_2', '/termini_1415_v01_2', '/termini_1516_v01_2', '/termini_1617_v01_2']
years = [2000, 2005, 2006, 2007, 2008, 2012, 2014, 2015, 2016]

termini = {}
for i,b in enumerate(basefiles):
    yr = years[i]
    fn = gl_termpos_fldr+b
    termini[yr] = read_termini(fn, yr) #creating dictionary for each year
    print len(termini[yr])
    
### Test earliest year of appearance for each glacier
#master_initial_termini = {}
#keys_05 = []
#keys_06 = []
#keys_07 = []
#keys_08 = []
#keys_12 = []
#
#for k in termini[2014].keys():
#    if k in termini[2000].keys():
#        master_initial_termini[k] = termini[2000][k]
#    elif k in termini[2005].keys():
#        print 'Glacier ID ' + str(k) + ' taken from year 2005'
#        master_initial_termini[k] = termini[2005][k]
#        keys_05.append(k)
#    elif k in termini[2006].keys():
#        print 'Glacier ID ' + str(k) + ' taken from year 2006'
#        master_initial_termini[k] = termini[2006][k]
#        keys_06.append(k)
#    elif k in termini[2007].keys():
#        print 'Glacier ID ' + str(k) + ' taken from year 2007'
#        master_initial_termini[k] = termini[2007][k]
#        keys_07.append(k)
#    elif k in termini[2008].keys():
#        print 'Glacier ID ' + str(k) + 'taken from year 2008'
#        master_initial_termini[k] = termini[2008][k]
#        keys_08.append(k)
#    elif k in termini[2012].keys():
#        print 'Glacier ID ' + str(k) + ' taken from year 2012'
#        master_initial_termini[k] = termini[2012][k]
#        keys_12.append(k)
#    else:
#        print 'Glacier ID ' + str(k) + ' not found before 2014'
 
       
#gl_termpos_fldr = 'Documents/GitHub/Data_unsynced/MEaSUREs-termini'
##print 'Reading in MEaSUREs terminus positions for year 2000'
##sf_termpos = shapefile.Reader(gl_termpos_fldr+'/termini_0001_v01_2') #Specify the base filename of the group of files that makes up a shapefile
#
#print 'Reading in MEaSUREs terminus positions for year 2014'
#sf_termpos_1415 = shapefile.Reader(gl_termpos_fldr+'/termini_1415_v01_2') #Specify the base filename of the group of files that makes up a shapefile
#
#print 'Reading in MEaSUREs terminus positions for year 2015'
#sf_termpos_1516 = shapefile.Reader(gl_termpos_fldr+'/termini_1516_v01_2') #Specify the base filename of the group of files that makes up a shapefile


## Find centroid rates of retreat for observed period 2006-2014
def Centroid_dLdt(termpts1, termpts2, time_interval, vx_func, vy_func):
    """Given two arrays of terminus-spanning coordinates, finds the distance between their centroids and converts to a retreat rate in the given time_interval.
    Returns retreat rate in distance per annum - check that termpts are in units of km (UTM) and time_interval in anni
    Inputs:
        termpts1 - points describing terminus at time t1
        termpts2 - points describing terminus at time t2
        time_interval - interval in anna, t2-t1
        vx_func - a 2D interpolated function for the x-component of ice velocity, used to determine sign of dL/dt
        vy_func - a 2D interpolated function for the y-component of ice velocity
    """
    term1 = geom.LineString(termpts1)
    term1_centr = term1.centroid
    centr1_coords = np.squeeze(list(term1_centr.coords))
    term2 = geom.LineString(termpts2)
    term2_centr = term2.centroid
    centr2_coords = np.squeeze(list(term2_centr.coords))
    
    disp_vector = centr1_coords - centr2_coords
    termchange = np.linalg.norm(disp_vector)
    dLdt_abs = termchange / time_interval
    
    v_vector = (vx_func(centr1_coords[0], centr1_coords[1]), vy_func(centr1_coords[0], centr1_coords[1]))
    dotprod = np.vdot(disp_vector, v_vector)
    sgn_dL = sign(dotprod) #sign of dot product indicates whether displacement of termini is parallel to velocity (advance) or antiparallel (retreat)
    
    return dLdt_abs*sgn_dL
    
retreat_rates = {yr:{} for yr in years[1::]}
for i in range(1, len(years)):
    current_termini = termini[years[i]]
    previous_termini = termini[years[i-1]]
    for gid in current_termini.keys():
        try:
            term1 = current_termini[gid]
            term2 = previous_termini[gid]
            retreat_rates[years[i]][gid] = Centroid_dLdt(term1, term2, time_interval=years[i]-years[i-1], vx_func=func_vxcomp, vy_func=func_vycomp) #calculate retreat rate at each glacier for each year 
        except KeyError:
            print 'No terminus found in year {} for glacier {}'.format(years[i], gid) # happens when glacier terminus is recorded in one year but not the previous

