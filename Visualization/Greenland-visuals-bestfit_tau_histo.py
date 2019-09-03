## Making histogram of best-fit yield strengths for Greenland networks
## 31 Oct 2018 - EHU

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
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


#--------------------------------------
# OPENING CSV OF STORED OPTIMAL YIELD STRENGTHS
#--------------------------------------

filename = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/Optimization_analysis/29Oct18.csv'

bins = range(50, 500, 5)
histo_dict = {k: 0 for k in bins}
error_dict = {}

tau_list = []

with open(filename, 'r') as f:
    header = f.readline()
    hdr = header.strip('\r\n')
    keys = hdr.split(',') #get names of variables

    lines = f.readlines()
    
    temp = []
    j =0
    for i, l in enumerate(lines):
        linstrip = l.strip('\r\n')
        parts = linstrip.split(',')
        
        GID = parts[0] #ID of glacier on this line
        bestfit_tau = int(float(parts[1]))/1000
        tau_type = parts[2]
        bed_elevation = float(parts[3])
        surface_elevation = float(parts[4])
        
        if surface_elevation - bed_elevation <0:
            error_dict[GID] = 'Negative ice thickness'
        else:
            tau_list.append(bestfit_tau)
            taubin = next(b for b in bins if b==bestfit_tau)
            try:
                histo_dict[taubin] +=1
            except: #blanket exception-catching because unsure what kind of error missing taubin will create
                pass
            
values = []
for b in bins:
    values.append(histo_dict[b])
#--------------------------------------
# HISTOGRAM BASED ON COUNTS
#--------------------------------------
plt.figure()
plt.bar(bins, values, align='center', alpha=0.5, facecolor='Indigo')
plt.xlabel('Optimal yield strength [kPa]', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.axes().set_xticklabels([0, 100, 200, 300, 400, 500], fontsize=18)
plt.axes().set_yticks([0, 10, 20, 30, 40])
plt.axes().set_yticklabels([0, 10, 20, 30, 40], fontsize=18)
plt.title('Greenland outlet glacier yield strengths found', fontsize=22)
plt.show()

plotting_bins = (50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
plt.figure()
plt.hist(tau_list, bins=plotting_bins, facecolor='Indigo', alpha=0.5)
plt.xlabel('Optimal yield strength [kPa]', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.axes().set_xticklabels(plotting_bins, fontsize=18)
plt.axes().set_yticks([0, 25, 50, 75])
plt.axes().set_yticklabels([0, 25, 50, 75], fontsize=18)
plt.title('Histogram of Greenland outlet glacier yield strengths', fontsize=22)
plt.show()