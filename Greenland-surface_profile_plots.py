## Plotting sequential surface profiles for a simulated glacier
## 21 June 2019  EHU

import numpy as np
import matplotlib.pyplot as plt
import csv
import shapefile
#import collections
import glob
#from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib.patches as mpatches
from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from plastic_utilities_v2 import *
from GL_model_tools import *
from flowline_class_hierarchy import *

###---------------------------------------
#### READ IN BED TOPO
###---------------------------------------
print 'Reading in bed topography'
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
fh.close()

#Smoothing bed
unsmoothB = B
smoothB = gaussian_filter(B, 2)
smoothS = gaussian_filter(S, 2) #17 Jan 19 - smoothing S as well for consistency with auto-selected networks

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothS.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])

###---------------------------------------
#### DEFINE WHERE THINGS LIVE, HOW TO READ
###---------------------------------------
flowlines_fpath = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/'
model_output_fpath = 'Documents/GitHub/Data_unsynced/Hindcasted_networks/'
yield_strength_fn = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/Optimization_analysis/bestfit_taus-B_S_smoothing-fromdate_2019-01-17.csv'

#function modified from Greenland-network-troubleshooting to read CSV of yield strengths
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
        
        data['Glacier_IDs'].append(int(parts[0]))
        data['Optimal_taus'].append(float(parts[1]))
        if read_yieldtype: #generally won't need this
            data['Yieldtype'].append(parts[2])
        else:
            pass
    
    return data

taudata = read_optimization_analysis(yield_strength_fn, read_yieldtype=True)


def ReadPlasticProfiles(gid, load_all=False):
    """Reads in data and model output related to a given glacier, and supplies a dictionary of surface profiles for plotting
    Arguments:
        gid: MEaSUREs glacier ID that identifies all related data
    Returns a fully-loaded plastic network from storage
    """
    flowlines_fn = glob.glob(flowlines_fpath+'Gld-autonetwork-GID{}-*.csv'.format(gid))[0] #using glob * to select files of any save date--there should be only one CSV of flowlines per GID
    output_fn = glob.glob(model_output_fpath+'GID{}-*.pickle'.format(gid))[0] #note that this will load anything of this GID - only one in hindcasted, but revisit in forward scenario projection
    tau_idx = (np.abs(np.asarray(taudata['Glacier_IDs']) - gid)).argmin()
    
    if load_all: #load all the flowlines
        coords = Flowline_CSV(flowlines_fn, has_width=True, flip_order=False)
        lines = []
        for j in range(len(coords)):
            line_j = Flowlines(coords=coords[j], index=j, name='GID {} line {}'.format(gid, j), has_width=True)
            lines.append(line_j)
        nw = PlasticNetwork(name='GID{}'.format(gid), init_type='Flowline', branches=lines, main_terminus=coords[0][0])
        nw.load_network(filename=output_fn, load_mainline_output=True, load_tributary_output=True)
    else: #load only one
        coords = Flowline_CSV(flowlines_fn, has_width=True, flip_order=False)[0]
        line = Flowline(coords=coords, index=0, name='GID {} line 0', has_width=True)
        nw = PlasticNetwork(name='GID{}'.format(gid), init_type='Flowline', branches=(line), main_terminus=coords[0])
        nw.load_network(filename=output_fn, load_mainline_output=True, load_tributary_output=False)
    
    nw.network_tau = taudata['Optimal_taus'][tau_idx]
    nw.network_yield_type = taudata['Yieldtype'][tau_idx]
    for fl in nw.flowlines:
        fl.optimal_tau = nw.network_tau
        fl.yield_type = nw.network_yield_type
    return nw
    
###---------------------------------------
#### DEFINE PLOTTING DEFAULTS
###---------------------------------------

testyears = arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
plotsize = (12,7)

def PlotSnapshots(network, years, plot_all=False, stored_profiles=False):
    """create snapshot for each year requested, only for 'main' flowline for now"""
    for year in years:
        if stored_profiles: #if model was run with output_heavy=True, profiles are already stored and don't need to be reconstructed
            profile_dict = network.model_output[0][year]
            xarr = profile_dict[0]
            SE_arr = profile_dict[1]
            bed_arr = profile_dict[2]
        else:
            #network.make_full_lines()
            #network.process_full_lines(B_interp, S_interp, H_interp)
            #network.remove_floating()
            #network.make_full_lines()
            network.process_full_lines(B_interp, S_interp, H_interp)
            output_dict = network.model_output[0] #output of line 0, the 'main' flowline
            idx = (np.abs(testyears - year)).argmin() # identify index of year requested
            terminus_position = output_dict['Termini'][idx]
            terminal_bed = network.flowlines[0].bed_function(terminus_position)
            Bingham_num = network.flowlines[0].Bingham_num(elev=0, thick=0) #ignore elevation/thickness dependence of Bingham number for this reconstruction
            profile_array = network.flowlines[0].plastic_profile(endpoint=terminus_position, hinit=BalanceThick(terminal_bed, Bingham_num))        
            xarr = profile_array[0]
            SE_arr = profile_array[1]
            bed_arr = profile_array[2]

            
        plt.figure(year, figsize=plotsize)
        plt.plot()
        plt.title('Glacier ID: {}, year {}'.format(network.name, year))
        plt.plot(10*xarr, bed_arr, color='Chocolate')
        plt.plot(10*xarr, SE_arr, color='Gainsboro')
        plt.fill_between(10*xarr, y1=SE_arr, y2 = bed_arr, color='Gainsboro', alpha=0.7)
        plt.fill_between(10*xarr, y1=bed_arr, y2=plt.axes().get_ylim()[0], color='Chocolate', alpha=0.7)
        plt.axes().set_xlim(left=10*xarr[-1], right=0)
        plt.axes().set_aspect(0.01)
        plt.axes().set_xlabel('Along-flowline distance [km]', fontsize=18)
        plt.axes().set_ylabel('Elevation [m a.s.l.]', fontsize=18)
        plt.tick_params(axis='both', labelsize=16)
        plt.show()


###---------------------------------------
#### GENERATE PLOTS
###---------------------------------------   

gid185 = ReadPlasticProfiles(185)
PlotSnapshots(gid185, (0, 2, 4, 6))