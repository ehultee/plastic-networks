#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Demonstrate temperature effect on Helheim
Created on Tue Oct  6 09:57:23 2020

@author: lizz
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import glob
from matplotlib import cm
from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
## Special import for SERMeQ modules
import sys
sys.path.insert(0, 'Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *

### Topography needed to remove floating points from saved coords
###
print('Reading in surface topography')
gl_bed_path ='/Users/lizz/Documents/GitHub/Data_unsynced/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
b_raw = fh.variables['bed'][:].copy() # bed topo
e_raw = fh.variables['errbed'][:].copy() # error in bed topo or ice thickness
thick_mask = fh.variables['mask'][:].copy()
ss = np.ma.masked_where(thick_mask !=2, s_raw)#mask values: 0=ocean, 1=ice-free land, 2=grounded ice, 3=floating ice, 4=non-Greenland land
hh = np.ma.masked_where(thick_mask !=2, h_raw) 
bb = b_raw #don't mask, to allow bed sampling from modern bathymetry (was subglacial in ~2006)
## Down-sampling
X = xx[::2]
Y = yy[::2]
S = ss[::2, ::2]
H = hh[::2, ::2] 
B = bb[::2, ::2]
E = e_raw[::2, ::2]
M = thick_mask[::2,::2]
fh.close()

#Smoothing bed and surface
smoothB = gaussian_filter(B, 2)
smoothS = gaussian_filter(S, 2)

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothS.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])


## Read in output from 2 Helheim simulations 
testyears = np.arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
min10C_output = {'Testyears': testyears}
min30C_output = {'Testyears': testyears}
fn_10 = '/Users/lizz/Documents/GitHub/Data_unsynced/Hindcasted_networks/GID175-2019-03-13-persistence-min10Cice-8a_dt025a.pickle'
fn_30 = '/Users/lizz/Dropbox/MBP-MBAir shared files/cold_ice/GID175-2020-07-06-persistence-min30Cice-8a_dt025a.pickle'
lightload(fn_10, glacier_name = 'GID175', output_dictionary = min10C_output)
lightload(fn_30, glacier_name = 'GID175', output_dictionary = min30C_output)

gl_termpos_fldr = '/Users/lizz/Documents/GitHub/Data_unsynced/MEaSUREs-termini'
basefiles = ['/termini_0607_v01_2', '/termini_0708_v01_2', '/termini_0809_v01_2', '/termini_1213_v01_2', '/termini_1415_v01_2', '/termini_1516_v01_2']
obs_years = [2006, 2007, 2008, 2012, 2014, 2015] #compare with term of hindcast, 2006-2014

termini = {}
for i,b in enumerate(basefiles):
    yr = obs_years[i]
    fn = gl_termpos_fldr+b
    termini[yr] = read_termini(fn, yr) #creating dictionary for each year
    print len(termini[yr])

projected_termini = []
print 'Reading in glacier ID 175'
filename = glob.glob('/Users/lizz/Documents/GitHub/Data_unsynced/Auto_selected-networks/Gld-autonetwork-GID175-date_*.csv')[0] #using glob * to select files of different run dates
coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
branch_0 = Branch(coords=coords_list[0], index=0, order=0) #saving central branch as main
branch_list = [branch_0]
nw = PlasticNetwork(name='GID175', init_type='Branch', branches=branch_list, main_terminus=branch_0.coords[0])
nw.make_full_lines()
nw.process_full_lines(B_interp, S_interp, H_interp)
nw.remove_floating()
mainline = LineString(nw.flowlines[0].coords)

for yr in obs_years:
    try:
        termpts = termini[yr][175] #get terminus points for each year
        t = projected_term_obs(termpts, mainline) #project onto main flowline
        r = retterm(termpts, mainline) #find most retreated point
        a = advterm(termpts, mainline) #find most advanced point
        print 'GID 175, t={}, r={}, a={}'.format(t, r, a)
        projected_termini.append(np.asarray((a, t, r))) #add these to dictionary of projected termini per glacier
    except KeyError:
        print 'No terminus found in year {} for GID 175.'.format(yr)
        projected_termini.append((0, np.nan, 0))    
    

###--------------------------------------
#### PLOTTING
###--------------------------------------   
plot_years = 2006+np.array(testyears)
sim_termini = min10C_output['GID175'][0]['Termini'][1::]
cold_termini = min30C_output['GID175'][0]['Termini'][1::]
obs_termini = np.asarray(projected_termini) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
obs_term_centr = obs_termini[:,1]
e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)

plt.figure('Main line terminus change, GID175')
plt.plot(plot_years, -0.001*np.array(sim_termini), linewidth=2, color='k', linestyle='-', label='Modelled')
plt.plot(plot_years, -0.001*np.array(cold_termini), linewidth=2, color='DarkGrey', linestyle='-', label='-30 C ice')
plt.errorbar(obs_years, -1*obs_term_centr, yerr = e, color='b', fmt='D')
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=18)
plt.axes().set(xlim=(2006, 2014.5), ylim=(-20, 2), 
               xticks=[2006, 2008, 2010, 2012, 2014], yticks=[-20, -15, -10, -5, 0],
               xticklabels=['2006', '2008', '2010', '2012', '2014'],
               yticklabels=['-20', '', '-10', '', '0'])
plt.axes().set_xlabel('Year', size=20)
plt.axes().set_ylabel('Terminus change [km]', size=20)
plt.axes().set_aspect(0.3)
plt.tight_layout()
plt.show()
