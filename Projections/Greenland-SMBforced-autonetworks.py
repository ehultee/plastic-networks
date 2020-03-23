# Forward projections with calving flux on Greenland glaciers, forced by surface mass balance
# Sept 2018  EHU
# Mar 2020 usage: large-scale 21st century projections

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj
import datetime
import glob
from matplotlib import cm
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
gl_bed_path ='Documents/GitHub/Data_unsynced/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
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

#print 'Reading in RCP 4.5 SMB'
#gl_smb_2081_path = 'Documents/GitHub/Data_unsynced/DMI-HIRHAM5_G6s2_ECEARTH_RCP45_2081_2100_gld_YM.nc'
#fh3 = Dataset(gl_smb_2081_path, mode='r')
#x_lon_81 = fh3.variables['lon'][:].copy() #x-coord (latlon)
#y_lat_81 = fh3.variables['lat'][:].copy() #y-coord (latlon)
##zs = fh2.variables['height'][:].copy() #height in m - is this surface elevation or SMB?
#ts_81 = fh3.variables['time'][:].copy()
#smb_2081_raw = fh3.variables['gld'][:].copy() #acc SMB in mm/day weq...need to convert
#fh3.close()


print 'Now transforming coordinate system of SMB'
wgs84 = pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
psn_gl = pyproj.Proj("+init=epsg:3413") # Polar Stereographic North used by BedMachine (as stated in NetDCF header)
xs, ys = pyproj.transform(wgs84, psn_gl, x_lon, y_lat)
#xs_81, ys_81 = pyproj.transform(wgs84, psn_gl, x_lon_81, y_lat_81)

#Xs = xs[0:,] #flattening; note that x-dimension is 402 according to file header
#Ys = ys[:,0] #flattening; note that y-dimension is 602 according to file header

smb_1980 = smb_raw[0][0]
smb_2014 = smb_raw[-1][0]
#smb_init_interpolated = interpolate.interp2d(ys, xs, smb_init, kind='linear')
Xmat, Ymat = np.meshgrid(X, Y)
regridded_smb_1980 = interpolate.griddata((xs.ravel(), ys.ravel()), smb_1980.ravel(), (Xmat, Ymat), method='nearest')
regridded_smb_2014 = interpolate.griddata((xs.ravel(), ys.ravel()), smb_2014.ravel(), (Xmat, Ymat), method='nearest')
SMB_1980 = interpolate.interp2d(X, Y, regridded_smb_1980, kind='linear')
SMB_2014 = interpolate.interp2d(X, Y, regridded_smb_2014, kind='linear')

### Hindcasting SMB: year-specific 2006-2014
#SMB_dict = {} #set up a dictionary of surface mass balance fields indexed by year
#for year in range(2006, 2015):
#    index = year - 2015 #so that 2014 will be smb_raw[-1], etc.
#    smb_year = smb_raw[index][0]
#    regridded_smb_year = interpolate.griddata((xs.ravel(), ys.ravel()), smb_year.ravel(), (Xmat, Ymat), method='nearest')
#    SMB_dict[year] = interpolate.interp2d(X, Y, regridded_smb_year, kind='linear')
    
#smb_2081_rcp4pt5 = smb_2081_raw[0]
#smb_2100_rcp4pt5 = smb_2081_raw[-1] #2100
#regridded_smb_2081 = interpolate.griddata((xs_81.ravel(), ys_81.ravel()), smb_2081_rcp4pt5.ravel(), (Xmat, Ymat), method='nearest')
#regridded_smb_2100 = interpolate.griddata((xs_81.ravel(), ys_81.ravel()), smb_2100_rcp4pt5.ravel(), (Xmat, Ymat), method='nearest')
#SMB_2081_RCP4pt5 = interpolate.interp2d(X, Y, regridded_smb_2081, kind='linear')
#SMB_2100_RCP4pt5 = interpolate.interp2d(X, Y, regridded_smb_2100, kind='linear')

## Add RCP 8.5 when available

##-------------------
### LOADING SAVED GLACIERS
##-------------------
print 'Reading in optimal yield strength dictionary'
optimal_taus_fpath = 'Documents/GitHub/Data_unsynced/Auto_selected-networks/Optimization_analysis/bestfit_taus-B_S_smoothing-fromdate_2019-01-17.csv'
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
glacier_ids = range(1,195) #MEaSUREs glacier IDs to process.
not_present = (93, 94, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 169) #glacier IDs missing from set
added_jan19 = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)
errors = (5, 17, 18, 19, 29, 51, 71, 92, 95, 97, 100, 101, 102, 106, 107, 108, 109, 110, 113, 115, 117, 120, 121, 134, 168, 171) #glacier IDs that crashed in hindcasting 12 Mar 2019 *or* showed network problems 21 May 2019
rmv = np.concatenate((not_present, errors))
for n in rmv:
    try:
        glacier_ids.remove(n)
    except ValueError:
        pass

## Simulation settings
testyears = arange(0, 100, 0.25)
start_year=2006 #determined by which MEaSUREs termini we used to initialize a given set
branch_sep_buffer = 10000/L0 #buffer between tributary intersections
db = True
#test_A, icetemp = 1.7E-24, 'min2C' # -2 C, warm ice
test_A, icetemp = 3.5E-25, 'min10C' # -10 C, good guess for Greenland
#test_A, icetemp = 3.7E-26, 'min30C' #-30 C, cold ice that should show slower response

scenario, SMB_i, SMB_l = 'persistence', SMB_2014, SMB_2014 #choose climate scenario - persistence of 1981-2014 climatology
#scenario, SMB_i, SMB_l = 'RCP4pt5', SMB_2014, SMB_2100_RCP4pt5 #or RCP 4.5
#scenario, SMB_i, SMB_l = 'RCP8pt5', SMB_2014, SMB_2100_RCP8pt5 #or RCP 8.5

#gids_totest = glacier_ids #test all
gids_totest = (12, 13, 14, 17) #test a selected subset
bad_gids = [] #store and write out poorly behaved glaciers

for gid in gids_totest:
    print 'Reading in glacier ID: '+str(gid)
    for gid in glaciers_simulated:
        filename = glob.glob('Documents/GitHub/Data_unsynced/Auto_selected-networks/Gld-autonetwork-GID{}-date_*'.format(gid))[0] #using glob * to select from files of multiple save dates

    
    coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
    nlines = len(coords_list)
    branch_0 = Branch(coords=coords_list[0], index=0, order=0) #saving central branch as main
    branch_list = [branch_0]
    if nlines>0:
        for l in range(1, nlines):
            branch_l = Branch(coords = coords_list[l], index=l, order=1, flows_to=0)
            branch_list.append(branch_l)
    nw = PlasticNetwork(name='GID'+str(gid), init_type='Branch', branches=branch_list, main_terminus=branch_0.coords[0])
    nw.make_full_lines()

    print 'Now processing glacier ID: '+str(gid)
    nw.process_full_lines(B_interp, S_interp, H_interp)
    nw.remove_floating()
    nw.make_full_lines()
    nw.process_full_lines(B_interp, S_interp, H_interp)

    print 'Now finding BedMachine terminus'
    idx, term_bm = next((i,c) for i,c in enumerate(nw.flowlines[0].coords) if NearestMaskVal(c[0], c[1])==2) #need to do error handling in case this is nowhere satisfied
    print 'BedMachine terminus is at index {}, coords {}'.format(idx, term_bm)
    #idx, term_bm = next((i, c) for i, c in enumerate(nw.flowlines[0].coords) if S_interp(c[0], c[1])>2.0) #testing with surface cutoff (2m) for now instead of retrieving mask
    term_arcval = ArcArray(nw.flowlines[0].coords)[idx]
    term_bed = nw.flowlines[0].bed_function(term_arcval)
    term_surface = nw.flowlines[0].surface_function(term_arcval)
    
    nw.network_tau = optimal_taus[gid][0]
    nw.network_yield_type = optimal_taus[gid][1]
    nw.network_ref_profiles()
    
    ## Simulations forced by SMB
    gl_smb_init = [0.001*SMB_i(nw.flowlines[0].coords[i,0], nw.flowlines[0].coords[i,1]) for i in range(len(nw.flowlines[0].coords))]
    nw.smb_alphadot = np.mean(gl_smb_init) #work on summing over all branches later
    print 'a_dot from SMB: {}'.format(nw.smb_alphadot)
    nw.terminus_adot = gl_smb_init[0]
    print 'Terminus a_dot: {}'.format(nw.terminus_adot)
    gl_smb_2100 = [0.001*365.26*SMB_l(nw.flowlines[0].coords[i,0], nw.flowlines[0].coords[i,1]) for i in range(len(nw.flowlines[0].coords))]
    nw.smb_2100_alphadot = np.mean(gl_smb_2100)
    steps_til_2100 = (2100-start_year) / mean(diff(testyears)) #length of linspace array that will determine forcing up to 2100 (HIRHAM-5 end 21st Century time)
    forcing_til_2100 = linspace(start=nw.smb_alphadot, stop=nw.smb_2100_alphadot, num=steps_til_2100) 
    if len(testyears)>steps_til_2100:
        steps_after_2100 = len(testyears) - steps_til_2100 #length of array that determines forcing after 2100
        forcing_after_2100 = np.full(shape=steps_after_2100, fill_value=nw.smb_2100_alphadot) #persist with 2100 SMB for timesteps after 2100
        variable_forcing = np.concatenate((forcing_til_2100, forcing_after_2100))
    else:
        variable_forcing = forcing_til_2100[:len(testyears)]  
    try:
        nw.terminus_time_evolve(testyears=testyears, alpha_dot_variable=variable_forcing, dL=1/L0, separation_buffer=10000/L0, has_smb=True, terminus_balance=nw.terminus_adot, submarine_melt = 0, debug_mode=db, rate_factor=test_A, output_heavy=False)
    except: ## move on to next glacier if this one won't behave
        bad_gids.append(gid)
        continue

    print 'Saving output for {}'.format(nw.name)
    fpath = 'Documents/GitHub/Data_unsynced/SERMeQ_output/'
    fn = str(nw.name)
    fn1 = fn.replace(" ", "")
    fn2 = fn1.replace("[", "-")
    fn3 = fn2.replace("/", "_")
    fn4 = fn3.replace("]", "")
    fn5 = fn4+'-{}-{}-{}ice-{}a_dt025a.pickle'.format(datetime.date.today(), scenario, icetemp, int(max(testyears)))
    nw.save_network(filename=fpath+fn5)
    
    
##Output errors to csv file
error_fn = 'Documents/GitHub/Data_unsynced/SERMeQ_output/error_gids-{}.csv'.format(datetime.date.today())
np.savetxt(error_fn, np.asarray(bad_gids), delimiter=',')