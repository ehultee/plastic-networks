## Extrapolate some flowlines forward and save as line -1
## 18 Nov 2019  EHU

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
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


def compute_seaward_coords(coords_for_network, distance=2000, spacing=50):
    """
    distance: distance in m to project forward. Default 2000
    spacing: distance in m between saved points
    """
    mainbranch_coords = coords_for_network[0]
    diffvec = np.asarray(mainbranch_coords[0])-np.asarray(mainbranch_coords[1]) #directional vector
    diffnorm = np.linalg.norm(diffvec)
    projector = (1./diffnorm)*np.array(diffvec)
    
    new_coords = [mainbranch_coords[0]+(d*projector) for d in arange(0, stop=distance+1, step=spacing)]
    return new_coords


def write_new_coords(new_coords_array, output_name, default_width=0):
    """
    new_coords_list: ndarray of seaward coordinates identified with compute_seaward_coords
    output_name: CSV filename to write coords to
    default_width: width in m to assign to this flowline.  Assigning a default width because we cannot use the Trace_wWidth method without surface velocity field
    """
    
    with open(output_name, 'wb') as csvfile: 
        linewriter = csv.writer(csvfile, delimiter=',')
        linewriter.writerow(['Line-number', 'x-coord', 'y-coord', 'width']) #use same format as was used to write flowlines with network_selection.py
        for n in range(len(new_coords_array)):
            linewriter.writerow([-1, new_coords_array[n][0], new_coords_array[n][1], default_width])
    

input_fpath = 'Documents/GitHub/Data_unsynced/Auto_selected-networks/Gld-autonetwork-GID'
output_fpath = 'Documents/GitHub/Data_unsynced/Auto_selected-networks/Seaward_coords/Gld-advnetwork-GID'

gids_in = (61, 64, 82, 83, 99, 130, 132, 139, 140, 141, 156, 157, 158, 161, 167, 170, 178, 179, 180, 184)
added_jan19 = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)


for gid in gids_in:
    print 'Reading in glacier ID: '+str(gid)
    if gid in added_jan19:
        filename = input_fpath+str(gid)+'-date_2019-01-10.csv'
    elif gid<160:
        filename = input_fpath+str(gid)+'-date_2018-10-03.csv'
    else:
        filename = input_fpath+str(gid)+'-date_2018-10-04.csv' #workaround because I ran these in batches and saved them with the date
    
    coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
    
    fwd_dist = 2000 #m, distance forward to project
    new_coords = compute_seaward_coords(coords_list, distance=fwd_dist)
    out_fn = output_fpath+'{}-fwd_{}_m.csv'.format(gid, fwd_dist)
    write_new_coords(new_coords, output_name=out_fn)