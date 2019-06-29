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
#### DEFINE WHERE THINGS LIVE, HOW TO READ
###---------------------------------------
flowlines_fpath = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/'
model_output_fpath = 'Documents/GitHub/Data_unsynced/Hindcasted_networks/'

def ReadPlasticProfiles(gid, load_all=False):
    """Reads in data and model output related to a given glacier, and supplies a dictionary of surface profiles for plotting
    Arguments:
        gid: MEaSUREs glacier ID that identifies all related data
    Returns a fully-loaded plastic network from storage
    """
    flowlines_fn = glob.glob(flowlines_fpath+'Gld_autonetwork-GID{}-*.csv'.format(gid))[0] #using glob * to select files of any save date--there should be only one CSV of flowlines per GID
    output_fn = glob.glob(model_output_fpath+'GID{}-*.pickle'.format(gid))[0] #note that this will load anything of this GID - only one in hindcasted, but revisit in forward scenario projection
    
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
        lines = (line,)
        nw = PlasticNetwork(name='GID{}'.format(gid), init_type='Flowline', branches=lines, main_terminus=coords[0])
        nw.load_network(filename=output_fn, load_mainline_output=True, load_tributary_output=False)
    
    return nw
    
    

