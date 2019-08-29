# Reading in HIRHAM SMB, reprojecting to match BedMachine/SENTINEL
# 24 May 2018  EHU

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj
import csv
import collections
#from matplotlib.colors import LogNorm
from matplotlib import cm
#from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from plastic_utilities_v2 import *
from GL_model_tools import *
from flowline_class_hierarchy import *

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
## Not down-sampling
#X = xx
#Y = yy
#S = ss
fh.close()

#Smoothing bed to check effect on dLdt
unsmoothB = B
smoothB = gaussian_filter(B, 2)
#B_processed = np.ma.masked_where(thick_mask !=2, smoothB)

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], S.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])

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

print 'Reading in surface mass balance'
gl_smb_path ='Documents/GitHub/plastic-networks/Data/DMI-HIRHAM5_GL2_ERAI_1980_2014_SMB_YM.nc'
fh2 = Dataset(gl_smb_path, mode='r')
x_lon = fh2.variables['lon'][:].copy() #x-coord (latlon)
y_lat = fh2.variables['lat'][:].copy() #y-coord (latlon)
#zs = fh2.variables['height'][:].copy() #height in m - is this surface elevation or SMB?
ts = fh2.variables['time'][:].copy()
smb_raw = fh2.variables['smb'][:].copy()
fh2.close()

print 'Now transforming coordinate system of SMB'
wgs84 = pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
psn_gl = pyproj.Proj("+init=epsg:3413") # Polar Stereographic North used by BedMachine (as stated in NetDCF header)
xs, ys = pyproj.transform(wgs84, psn_gl, x_lon, y_lat)

#Xs = xs[0:,] #flattening; note that x-dimension is 402 according to file header
#Ys = ys[:,0] #flattening; note that y-dimension is 602 according to file header

smb_init = smb_raw[0][0]
smb_latest = smb_raw[-1][0]
#smb_init_interpolated = interpolate.interp2d(ys, xs, smb_init, kind='linear')
Xmat, Ymat = np.meshgrid(X, Y)
regridded_smb_init = interpolate.griddata((xs.ravel(), ys.ravel()), smb_init.ravel(), (Xmat, Ymat), method='nearest')
regridded_smb_latest = interpolate.griddata((xs.ravel(), ys.ravel()), smb_latest.ravel(), (Xmat, Ymat), method='nearest')
SMB_i = interpolate.interp2d(X, Y, regridded_smb_init, kind='linear')
SMB_l = interpolate.interp2d(X, Y, regridded_smb_latest, kind='linear')


print 'Reading in surface elevation change'
gl_sec_path ='Documents/GitHub/plastic-networks/Data/CS2-SEC_2yr.nc'
fh3 = Dataset(gl_sec_path, mode='r')
x_sec = fh3.variables['x'][:].copy() #x-coord (polar stereo)
y_sec = fh3.variables['y'][:].copy() #y-coord (polar stereo)
t_sec = fh3.variables['t'][:].copy() #average time of slice (days since 1 JAN 2000)
sec_raw = fh3.variables['SEC'][:].copy()
fh3.close()

sec_i_masked = np.ma.masked_greater(sec_raw[:,:,0], 9000)
sec_i_excludemasked = np.ma.filled(sec_i_masked, fill_value=np.mean(sec_i_masked))
#sec_i_regrid = interpolate.griddata((x_sec.ravel(), y_sec.ravel()), sec_i_masked.ravel(), (Xmat, Ymat), method='nearest')
SEC_i = interpolate.RectBivariateSpline(x_sec, y_sec, sec_i_excludemasked)



jakcoords_main = Flowline_CSV('Documents/GitHub/plastic-networks/jakobshavn-mainline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_0 = Flowline(coords=jakcoords_main, index=0, name='Jak mainline', has_width=True)
Jakobshavn_main = PlasticNetwork(name='Jakobshavn Isbrae [main/south]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_main[0])
Jakobshavn_main.load_network(filename='JakobshavnIsbrae-main_south.pickle')

Jakobshavn_main.process_full_lines(B_interp, S_interp, H_interp)
for fln in Jakobshavn_main.flowlines:
    fln.yield_type = Jakobshavn_main.network_yield_type
    fln.optimal_tau = Jakobshavn_main.network_tau
Jakobshavn_main.network_ref_profiles()

#Jakobshavn_smb_l = [0.001*SMB_l(Jakobshavn_main.flowlines[0].coords[i,0], Jakobshavn_main.flowlines[0].coords[i,1]) for i in range(len(Jakobshavn_main.flowlines[0].coords))] #multiplied by 0.001 to convert from mm to m
#Jak_smb_alphadot = np.mean(Jakobshavn_smb_l)
#Jak_terminus_adot = Jakobshavn_smb_l[0]
#
Jak_sec_mainline = np.asarray([SEC_i(Jakobshavn_main.flowlines[0].coords[i,0], Jakobshavn_main.flowlines[0].coords[i,1]) for i in range(len(Jakobshavn_main.flowlines[0].coords))])
away_from_edge = np.argmin(Jak_sec_mainline)
Jak_sec_alphadot = np.mean(Jak_sec_mainline[away_from_edge::])
Jak_terminus_sec = float(min(np.asarray(Jak_sec_mainline).flatten())) #using min because close to edge, values get disrupted by mask interpolation

Jakobshavn_main.terminus_time_evolve(testyears=arange(10), alpha_dot=Jak_sec_alphadot/Jakobshavn_main.H0, has_smb=True, terminus_balance=Jak_terminus_sec/Jakobshavn_main.H0, debug_mode=True)