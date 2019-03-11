## Investigating problems in empty network tracing from MEaSUREs termini
## 6 Nov 2018  EHU
## Edited 2 Jan 2019 - investigating networks that claim very small or negative terminus elevation
import numpy as np
from scipy import interpolate, ndimage, misc
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
from plastic_utilities_v2 import *
from GL_model_tools import *
from flowline_class_hierarchy import *

##-------------------------------
##   READING IN & PROCESSING DATA
##-------------------------------

###Reading in velocities -- function lifted from Greenland-vel-compositing.py
###Function to read MEaSUREs velocity GeoTIFFs
#def read_velocities(filename, return_grid=True, return_proj=False):
#    """Extract x, y, v from a MEaSUREs GeoTIFF"""
#    ds = gdal.Open(filename)
#    #Get dimensions
#    nc = ds.RasterXSize
#    nr = ds.RasterYSize
#    
#    geotransform = ds.GetGeoTransform()
#    xOrigin = geotransform[0]
#    xPix = geotransform[1] #pixel width in x-direction
#    yOrigin = geotransform[3]
#    yPix = geotransform[5] #pixel height in y-direction
#    
#    lons = xOrigin + np.arange(0, nc)*xPix
#    lats = yOrigin + np.arange(0, nr)*yPix
#    
#    x, y = np.meshgrid(lons, lats)
#    
#    vband = ds.GetRasterBand(1)
#    varr = vband.ReadAsArray()
#    
#    #if return_grid and return_proj:
#    #    return x, y, varr, ds.GetProjection()
#    #elif return_grid:
#    if return_grid:
#        return x, y, varr
#    else: 
#        return varr
#
### Read in MEaSUREs velocity composite
#print 'Reading MEaSUREs velocities'
#x_comp, y_comp, v_comp_raw = read_velocities('Documents/GitHub/Data_unsynced/gld-velocity-composite.tif')
#vx_comp_raw = read_velocities('Documents/GitHub/Data_unsynced/gld-x_velocity-composite.tif', return_grid=False)
#vy_comp_raw = read_velocities('Documents/GitHub/Data_unsynced/gld-y_velocity-composite.tif', return_grid=False)
#v_comp = np.ma.masked_invalid(v_comp_raw)
#vx_comp = np.ma.masked_invalid(vx_comp_raw)
#vy_comp = np.ma.masked_invalid(vy_comp_raw)
#v_excludemasked = np.ma.filled(v_comp, fill_value=0)
#vx_excludemasked = np.ma.filled(vx_comp, fill_value=0)
#vy_excludemasked = np.ma.filled(vy_comp, fill_value=0)
#days_per_annum = 365.242 #need to convert units of MEaSUREs velocity to align with what we used from Sentinel before
#v_a2d = np.array(v_excludemasked) / days_per_annum
#vx_a2d = np.array(vx_excludemasked) / days_per_annum
#vy_a2d = np.array(vy_excludemasked) / days_per_annum

## Read in BedMachine topography
print 'Reading in surface topography'
gl_bed_path ='Documents/1. Research/2. Flowline networks/Model/Data/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
b_raw = fh.variables['bed'][:].copy() # bed topo
thick_mask = fh.variables['mask'][:].copy()
ss = np.ma.masked_where(thick_mask !=2, s_raw)#mask values: 0=ocean, 1=ice-free land, 2=grounded ice, 3=floating ice, 4=non-Greenland land
hh = np.ma.masked_where(thick_mask !=2, h_raw) 
#bb = np.ma.masked_where(thick_mask !=2, b_raw)
bb = b_raw #don't mask, to allow bed sampling from modern bathymetry (was subglacial in ~2006)
## Down-sampling
X = xx[::2]
Y = yy[::2]
S = ss[::2, ::2]
H = hh[::2, ::2] 
B = bb[::2, ::2]
M = thick_mask[::2,::2]
## Not down-sampling
#X = xx
#Y = yy
#S = ss
fh.close()

unsmoothB = B
smoothing_sigma = 2 #chooses how much smoothing to apply
smoothB = gaussian_filter(B, smoothing_sigma)

S_new = np.add(B, H) #alternative surface elevation: bed elevation + ice thickness

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], S.T[::, ::-1]) #interpolating surface elevation provided
#S_interp = interpolate.RectBivariateSpline(X, Y[::-1], S_new.T[::, ::-1]) #interpolating bed + thickness
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])



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


## New function to read in CSV results of optimal-yield-strength analysis
def read_optimization_analysis(filename, read_yieldtype=False):
    """Read a CSV file listing optimal values of yield strength for auto-selected Greenland glaciers
    Input: 
        filename
    Default arg: 
        read_yieldtype=False: determines whether we want to read and save the yield type (constant vs. Coulomb variable)
    Output: 
        Dictionary of lists including
        -Glacier ID (referenced to MEaSUREs Greenland outlets)
        -Optimal tau_y
        -Terminal bed
        -Terminal surface elevation
    """
    
    f = open(filename, 'r')
    
    header = f.readline() #header line
    #hdr = header.strip('\r\n')
    #keys = hdr.split(',') #get names of columns
    #data = {k: [] for k in keys}
    data = {'Glacier_IDs': [], #shorter keys than the names stored in CSV header
    'Optimal_taus': [],
    'Yieldtype': [],
    'Terminal_bed': [],
    'Terminal_SE': [],
    'Terminal_H': []} #adding field for ice thickness
    
    lines = f.readlines()
    f.close
    
    for i, l in enumerate(lines):
        linstrip = l.strip('\r\n')
        parts = linstrip.split(',')
        
        bed_el = float(parts[3])
        surf_el = float(parts[4])
        
        data['Glacier_IDs'].append(int(parts[0]))
        data['Optimal_taus'].append(float(parts[1]))
        if read_yieldtype: #generally won't need this
            data['Yieldtype'].append(parts[2])
        else:
            pass
        data['Terminal_bed'].append(bed_el)
        data['Terminal_SE'].append(surf_el)
        data['Terminal_H'].append(surf_el-bed_el) #calculating ice thickness from reported bed and surface elevations
    
    return data

analysis_fn = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/Optimization_analysis/bestfit_taus-B_S_smoothing-fromdate_2019-01-17.csv'
opt_data = read_optimization_analysis(analysis_fn)
#IDs_neg_termini = [g for i,g in enumerate(opt_data['Glacier_IDs']) if opt_data['Terminal_SE'][i]<0]
#IDs_thin_termini = [g for i,g in enumerate(opt_data['Glacier_IDs']) if opt_data['Terminal_H'][i]<10]
#taus_neg_termini = [t for i, t in enumerate(opt_data['Optimal_taus']) if opt_data['Terminal_SE'][i]<0]   
#SEs_neg_termini = [se for i, se in enumerate(opt_data['Terminal_SE']) if opt_data['Terminal_SE'][i]<0]

IDs_foranalysis = (9,)
taus_foranalysis = [t for i, t in enumerate(opt_data['Optimal_taus']) if opt_data['Glacier_IDs'][i] in IDs_foranalysis]

## Read in the networks of glaciers with problematic termini--based on Greenland-network_processing_routine.py
base_fpath = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/Gld-autonetwork-GID'
flowlines_foranalysis = []
load_all = True #decide whether to analyse all flowlines of a network or only the central branch
for j, gid in enumerate(IDs_foranalysis):
    print 'Reading in glacier ID: '+str(gid)
    if gid<160:
        filename = base_fpath+str(gid)+'-date_2018-10-03.csv'
    else:
        filename = base_fpath+str(gid)+'-date_2018-10-04.csv' #workaround because I ran these in batches and saved them with the date
    
    coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
    nlines = len(coords_list)
    branch_0 = Branch(coords=coords_list[0], index=0, order=0) #saving central branch as main
    branch_list = [branch_0]
    if load_all:
        for k in range(1, nlines):
            branch = Branch(coords=coords_list[k], index=k)
            branch_list.append(branch)
    else:
        pass
    nw = PlasticNetwork(name='GID'+str(gid), init_type='Branch', branches=branch_list, main_terminus=branch_0.coords[0])
    #fl = Flowline(coords=coords_list[0], name='GID'+str(gid), index=0) #saving central branch as main
    print 'Now processing glacier ID: '+str(gid)
    nw.make_full_lines()
    nw.process_full_lines(B_interp, S_interp, H_interp)
    nw.remove_floating()
    nw.make_full_lines()
    nw.process_full_lines(B_interp, S_interp, H_interp)
    nw.optimal_tau = taus_foranalysis[j]
    flowlines_foranalysis.append(nw)
    

##-------------------------------
##   EXPLORATORY CALCULATIONS
##-------------------------------

## How far away are the MEaSUREs termini of networks stored empty from a valid MEaSUREs velocity composite value?
#IDs_stored_empty = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)
#def find_nearest_edge(masked_field, xgrid, ygrid, point, searchbuffer = 3):
#    '''Edge detection for masked quantities, e.g. surface, velocity
#    Default arg:
#        searchbuffer: buffer of surrounding pixels to search.  Default is 3--so will search a 6x6 pixel box
#    '''
#    sobel_x = ndimage.sobel(masked_field, axis=0, mode='constant')
#    sobel_y = ndimage.sobel(masked_field, axis=1, mode='constant')
#    sobel_im = np.hypot(sobel_x, sobel_y)
#    
#    nearest_gridx = np.argmin(xgrid - point[0])
#    nearest_gridy = np.argmin(ygrid - point[1]) #find indices of the point in the grid closest to the given coords
#    
#    search_indices = ((nearest_gridx + i, nearest_gridy + j) for i in range(-1*searchbuffer, searchbuffer) for j in range(-1*searchbuffer, searchbuffer))
#    search_pts = xgrid    
#    
    
    
  
##-------------------------------
##   EXPLORATORY PLOTTING
##-------------------------------
### Do networks stored empty have termini that lie 'off the map' of velocity for tracing?
#plt.figure()
#plt.contourf(x_comp, y_comp[::], v_comp[::-1, ::])
#for gid in IDs_stored_empty:
#    plt.scatter(termini_init[gid][:,0], termini_init[gid][:,1], marker='*', color='r')
#plt.show()
### Yes, they do.


## Do the surface profiles of glaciers with small/negative terminus elevation look wonky in other ways?
#IDs_neg_termini = (3, 9, 13, 19, 20, 22, 29, 30, ) #Glacier IDs with negative terminus elevation according to bestfit_taus-fromdate_2018-11-02.csv (extracted surface elevation from interpolation of B+H)
#for nw in flowlines_foranalysis: 
#    fl = nw.flowlines[0]
#    xarr = linspace(0, fl.length, num=1000)   
#    plt.figure()
#    plt.title('Glacier ID: {}.  Yield strength: {} kPa. Terminal SE: {} m'.format(nw.name, 0.001*nw.optimal_tau, fl.surface_function(0)))
#    plt.plot(10*xarr, fl.bed_function(xarr), color='Chocolate')
#    plt.plot(10*xarr, fl.surface_function(xarr), color='Gainsboro')
#    plt.axes().set_xlim(left=10*xarr[-1], right=0)
#    plt.show()
## Yes, many do look wonky.  Many also report a thickness (from fl.thickness_function) different from their surface-bed. More structured analysis to come.


## Do glaciers with surprisingly large sea level contributions have bed features that affect their stability?
for nw in flowlines_foranalysis:
    for j in range(len(nw.flowlines)): 
        fl = nw.flowlines[j]
        xarr = linspace(0, fl.length, num=1000)   
        plt.figure('GID{}, flowline {}'.format(nw.name, j))
        plt.title('Glacier ID: {}.  Terminal SE: {} m.  Bed smoothing sigma = {}'.format(nw.name, fl.surface_function(0), smoothing_sigma))
        plt.plot(10*xarr, fl.bed_function(xarr), color='Chocolate')
        plt.plot(10*xarr, fl.surface_function(xarr), color='Gainsboro')
        plt.fill_between(10*xarr, y1=fl.surface_function(xarr), y2 = fl.bed_function(xarr), color='Gainsboro', alpha=0.7)
        plt.fill_between(10*xarr, y1=fl.bed_function(xarr), y2=plt.axes().get_ylim()[0], color='Chocolate', alpha=0.7)
        plt.axes().set_xlim(left=10*xarr[-1], right=0)
        plt.axes().set_aspect(0.01)
        plt.axes().set_xlabel('Along-flowline distance [km]', fontsize=18)
        plt.axes().set_ylabel('Elevation [m a.s.l.]', fontsize=18)
        plt.tick_params(axis='both', labelsize=16)
        plt.show()
## On GID9, some prominent lumps that are much steeper on downstream side; retrograde bed far upstream
## GID10 shows retrograde bed as well, very thin nose sitting in water