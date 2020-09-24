#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Extract bed topography error along Greenland flowlines
Created on Thu Sep 24 10:41:48 2020

@author: lizz
"""

from netCDF4 import Dataset
import numpy as np
import glob
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import csv as csv


### Topography needed to remove floating points from saved coords
print('Reading in surface topography')
gl_bed_path ='/Users/lizz/Documents/GitHub/Data_unsynced/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
b_raw = fh.variables['bed'][:].copy() # bed topo
e_raw = fh.variables['errbed'][:].copy()
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
E = e_raw[::2, ::2]
M = thick_mask[::2,::2]
## Not down-sampling
#X = xx
#Y = yy
#S = ss
fh.close()

#Smoothing bed for sampling
smoothB = gaussian_filter(B, 2)
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])
E_interp = interpolate.RectBivariateSpline(X, Y[::-1], E.T[::,::-1])

## Read in flowline coordinates and write out error
glacier_ids = range(1,195) #MEaSUREs glacier IDs to process.
not_present = (93, 94, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 169) #glacier IDs missing from set
errors = (5, 17, 18, 19, 29, 51, 71, 92, 95, 97, 100, 101, 102, 106, 107, 108, 109, 110, 113, 115, 117, 120, 121, 134, 168, 171) #glacier IDs that crashed in hindcasting 12 Mar 2019 *or* showed network problems 21 May 2019
rmv = np.concatenate((not_present, errors))
for n in rmv:
    try:
        glacier_ids.remove(n)
    except ValueError:
        pass

in_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Auto_selected-networks/Gld-autonetwork-GID'
out_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Autonetwork_data_QC/Gld-lines_errors-GID'
for gid in glacier_ids: 
    input_fn = glob.glob(in_fpath+'{}-date_*.csv'.format(gid))[0] #using glob * to select files of different run dates
    output_fn = out_fpath+'{}.csv'.format(gid)
    with open(input_fn, 'r') as read_obj, \
        open(output_fn, 'w') as write_obj:
            reader = csv.DictReader(read_obj)
            headers = reader.fieldnames
            headers.extend(['bed-topo', 'bed-error'])
            writer = csv.DictWriter(write_obj, fieldnames=headers)
            for row in reader:
                ln = row['Line-number']
                x = row['x-coord']
                y = row['y-coord']
                w = row['width'] 
                b = B_interp(x,y) 
                e = E_interp(x,y)
                row['bed-topo'] = float(b)
                row['bed-error'] = float(e)
                writer.writerow(row)
    