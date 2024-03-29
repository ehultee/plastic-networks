## Script to track net advance on seaward-projected flowlines

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj
import csv
import collections
import glob
import datetime
#from matplotlib.colors import LogNorm
from matplotlib import cm
#from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
## Special import for SERMeQ modules
import sys
sys.path.insert(0, 'Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *


##-------------------
### READING IN BED
### COMMENT OUT IF DATA IS ALREADY READ IN TO YOUR SESSION
##-------------------

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
bb = np.ma.masked_where(thick_mask !=2, b_raw)
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

#Smoothing bed
unsmoothB = B
smoothB = gaussian_filter(B, 2)
smoothS = gaussian_filter(S, 2) #17 Jan 19 - smoothing S as well for consistency with auto-selected networks
#B_processed = np.ma.masked_where(thick_mask !=2, smoothB)

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothS.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])
def NearestMaskVal(x0,y0):
    x_idx = np.abs(X-x0).argmin()
    y_idx = np.abs(Y-y0).argmin()
    return M[y_idx, x_idx]


print 'Reading in surface mass balance from 1981-2010 climatology'
gl_smb_path ='Documents/GitHub/Data_unsynced/HIRHAM5-SMB/DMI-HIRHAM5_GL2_ERAI_1980_2014_SMB_YM.nc'
fh2 = Dataset(gl_smb_path, mode='r')
x_lon = fh2.variables['lon'][:].copy() #x-coord (latlon)
y_lat = fh2.variables['lat'][:].copy() #y-coord (latlon)
#zs = fh2.variables['height'][:].copy() #height in m - is this surface elevation or SMB?
ts = fh2.variables['time'][:].copy()
smb_raw = fh2.variables['smb'][:].copy()
fh2.close()


print 'Now transforming coordinate system of SMB'
wgs84 = pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
psn_gl = pyproj.Proj("+init=epsg:3413") # Polar Stereographic North used by BedMachine (as stated in NetCDF header)
xs, ys = pyproj.transform(wgs84, psn_gl, x_lon, y_lat)
Xmat, Ymat = np.meshgrid(X, Y)

## Hindcasting SMB: year-specific 2006-2014
SMB_dict = {} #set up a dictionary of surface mass balance fields indexed by year
for year in range(2006, 2015):
    index = year - 2015 #so that 2014 will be smb_raw[-1], etc.
    smb_year = smb_raw[index][0]
    regridded_smb_year = interpolate.griddata((xs.ravel(), ys.ravel()), smb_year.ravel(), (Xmat, Ymat), method='nearest')
    SMB_dict[year] = interpolate.interp2d(X, Y, regridded_smb_year, kind='linear')

##-------------------
### LOADING SAVED GLACIERS
##-------------------
print 'Reading in optimal yield strength dictionary'
optimal_taus_fpath = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/Optimization_analysis/bestfit_taus-B_S_smoothing-fromdate_2019-01-17.csv'
f_ot = open(optimal_taus_fpath, 'r')
header = f_ot.readline()
hdr = header.strip('\r\n')
optimal_taus = {}
lines = f_ot.readlines()
for i, l in enumerate(lines):
    linstrip = l.strip('\r\n')
    parts = linstrip.split(',')
    gid = int(parts[0]) #MEaSUREs glacier ID is first
    tau_y = float(parts[1]) #yield strength in Pa
    yieldtype = parts[2] #'constant' or 'variable' string
    optimal_taus[gid] = (tau_y, yieldtype)

##-------------------
### SIMULATING GLACIER NETWORKS
##-------------------

print 'Loading and simulating glacier networks'
gids_totest = (61, 64, 82, 83, 99, 130, 132, 139, 140, 141, 156, 157, 158, 161, 167, 170, 178, 179, 180, 184) 

base_fpath = 'Documents/GitHub/Data_unsynced/Auto_selected-networks/Gld-autonetwork-GID'
seaward_coords_fpath = 'Documents/GitHub/Data_unsynced/Auto_selected-networks/Seaward_coords/Gld-advnetwork-GID' 

## Simulation settings
testyears = arange(0, 9, 0.25) #test only 2006-2015, for comparison
start_year=2006 #determined by which MEaSUREs termini we used to initialize a given set
branch_sep_buffer = 10000/L0 #buffer between tributary intersections
db = True
test_A, icetemp = 3.5E-25, 'min10C' # -10 C, good guess for Greenland
scenario = 'persistence'

output_heavy = False #pref for output--should be False for efficient running, True for saving full surface profiles
network_output = []
bad_gids = []

for gid in gids_totest:
    print 'Reading in glacier ID: '+str(gid)
    filename = glob.glob(base_fpath+'{}-date_*.csv'.format(gid))[0] #using glob * to select files of multiple run dates
    coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
    nlines = len(coords_list)
    seaward_fn = seaward_coords_fpath+'{}-fwd_2000_m.csv'.format(gid)
    seaward_coords = Flowline_CSV(seaward_fn, has_width=True, flip_order=True)[0]
    
    branch_0 = Branch(coords=np.concatenate((seaward_coords, coords_list[0])), index=0, order=0) #saving central branch as main
    branch_list = [branch_0]

    print 'Now processing glacier ID: '+str(gid)
    nw = PlasticNetwork(name='GID'+str(gid), init_type='Branch', branches=branch_list, main_terminus=coords_list[0][0])
    nw.make_full_lines()
    nw.process_full_lines(B_interp, S_interp, H_interp)    
    nw.network_tau = optimal_taus[gid][0]
    nw.network_yield_type = optimal_taus[gid][1]
    nw.network_ref_profiles()
    
    ## Simulations forced by SMB
    #gl_smb_init = [0.001*365.26*SMB_i(nw.flowlines[0].coords[i,0], nw.flowlines[0].coords[i,1]) for i in range(len(nw.flowlines[0].coords))]
    time_varying_smb = [[0.001*(1000/nw.rho_ice)*SMB_dict[int(floor(yr))](nw.flowlines[0].coords[i,0], nw.flowlines[0].coords[i,1]) for i in range(len(nw.flowlines[0].coords))] for yr in start_year+testyears]
    catchment_smb_vals = [np.mean(time_varying_smb[i]) for i in range(len(time_varying_smb))] #because grid spacing is very fine, mean of all points should approximate integral/length
    print 'shape(catchment_smb_vals) = {}'.format(shape(catchment_smb_vals))
    nw.smb_alphadot = catchment_smb_vals[0] #initial SMB
    print 'a_dot from SMB: {}'.format(nw.smb_alphadot)
    nw.terminus_adot = time_varying_smb[0][0]
    print 'Terminus a_dot: {}'.format(nw.terminus_adot)
    try:
        nw.terminus_time_evolve(testyears=testyears, initial_term=max(ArcArray(seaward_coords)), alpha_dot_variable=catchment_smb_vals, dL=1/L0, separation_buffer=10000/L0, has_smb=True, terminus_balance=nw.terminus_adot, submarine_melt = 0, debug_mode=db, rate_factor=test_A, output_heavy=output_heavy)
    except:
        bad_gids.append(gid)   
        continue 
    
    print 'Saving output for {}'.format(nw.name)
    fn = str(nw.name)
    fn1 = fn.replace(" ", "")
    fn2 = fn1.replace("[", "-")
    fn3 = fn2.replace("/", "_")
    fn4 = fn3.replace("]", "")
    fn5 = 'Documents/GitHub/Data_unsynced/Hindcasted_networks/advance_test/'+fn4+'-{}-{}-{}ice-{}a_dt025a.pickle'.format(datetime.date.today(), scenario, icetemp, int(max(testyears)))
    nw.save_network(filename=fn5)

    #network_output.append(nw.model_output)

##Output errors to csv file
error_fn = 'Documents/GitHub/Data_unsynced/Hindcasted_networks/error_gids-{}.csv'.format(datetime.date.today())
np.savetxt(error_fn, np.asarray(bad_gids), delimiter=',')