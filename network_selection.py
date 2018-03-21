# Formalizing flowline selection algorithm for Greenland glaciers
# Adding functionality to find widths
# 16 March 2018  EHU

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
from scipy.ndimage import gaussian_filter

## Reading in BedMachine map
print 'Reading in surface topography'
gl_bed_path ='Documents/1. Research/2. Flowline networks/Model/Data/BedMachine-Greenland/MCdataset-2015-04-27.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
#h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
thick_mask = fh.variables['mask'][:].copy()
ss = np.ma.masked_where(thick_mask !=2, s_raw)
#hh = np.ma.masked_where(thick_mask !=2, h_raw) #mask values: 0=ocean, 1=ice-free land, 2=grounded ice, 3=floating ice, 4=non-Greenland land
## Down-sampling
#X = xx[::2]
#Y = yy[::2]
#S = ss[::2, ::2]
##H = hh[::2, ::2] 
## Not down-sampling
X = xx
Y = yy
#S = ss
fh.close()

## Smoothing surface elevation
unsmoothS = ss
smoothS = gaussian_filter(ss, 2)
S = np.ma.masked_where(thick_mask !=2, smoothS)

## Reading in SENTINEL velocity map
print 'Now reading in (vector) velocity map'
v_path = 'Documents/1. Research/2. Flowline networks/Model/Data/ESA-Greenland/greenland_iv_500m_s1_20161223_20170227_v1_0.nc'
fh2 = Dataset(v_path, mode='r')
xv = fh2.variables['x'][:].copy()
yv = fh2.variables['y'][:].copy()
#yv = yv_flipped[::-1]
v_raw = fh2.variables['land_ice_surface_velocity_magnitude'][:].copy() #this is v(y, x)
vx_raw = fh2.variables['land_ice_surface_easting_velocity'][:].copy()
vy_raw =fh2.variables['land_ice_surface_northing_velocity'][:].copy()
v_upper = np.ma.masked_greater(v_raw, 10000)
vx_upper = np.ma.masked_greater(vx_raw, 10000)
vy_upper = np.ma.masked_greater(vy_raw, 10000)
fh2.close()

## Interpolate SENTINEL and sample at BedMachine points
print 'Now interpolating to same grid'
vf_x = interpolate.interp2d(xv, yv[::-1], vx_upper[::-1,::])
vf_y = interpolate.interp2d(xv, yv[::-1], vy_upper[::-1,::])
vf = interpolate.interp2d(xv, yv[::-1], v_upper[::-1, ::])



## Functions to search for flowlines
def check_sides(x, y, Vx = vf_x, Vy = vf_y, V=vf, side_cutoff = 0.25, dw=20):
    """Reports coordinates of glacier 'sides' by searching velocity of points normal to direction of flow and comparing to a cutoff value. 
    x: x-coordinate of current point along flowline
    y: y-coordinate of current point along flowline
    Vx: 2D interpolated function of the x-velocity field
    Vy: 2D interpolated function of the y-velocity field
    V: 2d interpolated function of magnitude of surface velocity
    side_cutoff_v: the minimum velocity, in units of V field, for ice to be included in branch.  Default is 1 m/d.
    dw: spatial step size, in meters, for checking across flow.  Default is 20 m.
    """
    vx = Vx(x, y)
    vy = Vy(x, y)
    vm = V(x, y)
    
    nr = 1/vm*np.array([-vy, vx]) #right normal
    nl = 1/vm*np.array([vy, -vx]) #left normal
    
    
    Cr = 1 # counter for width increments on right
    Cl = 1 # counter for width increments on left
    v_rightside = vm
    v_leftside = vm
    
    while v_rightside > side_cutoff:
        x_r = float(x + Cr*dw*nr[0])
        y_r = float(y + Cr*dw*nr[1])
        v_rightside = V(x_r, y_r)
        Cr += 1
    while v_leftside > side_cutoff:
        x_l = float(x + Cl*dw*nl[0])
        y_l = float(y + Cl*dw*nl[1])
        v_leftside = V(x_l, y_l)
        Cl += 1
    
    print Cr, Cl
    print v_leftside, v_rightside
    
    #width = (Cr+Cl-2)*dw
    rightcoords = (x_r, y_r)
    leftcoords = (x_l, y_l)
    
    return rightcoords, leftcoords


def Trace_wWidth(startcoord_x, startcoord_y, trace_up=False, xarr=X, yarr=Y, Vx = vf_x, Vy = vf_y, V = vf, dx=150, side_cutoff = 0.25, dw=20, output_sidecoords=False):
    """Traces flowlines down from a tributary head or up from terminus.  Stores evenly spaced coordinates and width at each point.
    startcoord_x: x-coordinate of the starting point
    startcoord_y: y-coordinate of the starting point
    trace_up: indicates whether tracing upstream (True) or downstream (False)
    xarr: x-coordinates of the velocity field
    yarr: y-cordinates of the velocity field
    Vx: 2D interpolated function of the x-velocity field
    Vy: 2D interpolated function of the y-velocity field
    V: 2d interpolated function of magnitude of surface velocity
    dx: spatial step size, in meters, for stepping algorithm.  Default is 150 m.
    side_cutoff_v: the minimum velocity, in units of V field, for ice to be included in branch.  Default is 1m/d.
    dw: spatial step size, in meters, for checking width.  Default is 20 m.
    output_sidecoords: 
    """
    outarr = []
    rightside = []
    leftside = []
    
    side_tolerance = 3*dx
    
    
    currentpt = (startcoord_x, startcoord_y)
    current_r = np.asarray(check_sides(startcoord_x, startcoord_y, Vx, Vy, V, side_cutoff=side_cutoff)[0])
    current_l = np.asarray(check_sides(startcoord_x, startcoord_y, Vx, Vy, V, side_cutoff=side_cutoff)[1])
    nstep = 0
    
    while nstep<1000:
        vmax = 50 # m/day, the maximum velocity we believe
        
        vx = Vx(currentpt[0], currentpt[1])
        vy = Vy(currentpt[0], currentpt[1])
        vm = V(currentpt[0], currentpt[1]) #speed - need to have V-field for width check
        
        if vm > vmax:
            print 'Speed exceeds maximum recognised.  Exiting step routine.'
            break
        else:
            r, l = check_sides(currentpt[0], currentpt[1], Vx, Vy, V, side_cutoff = side_cutoff)
            r = np.asarray(r)
            l = np.asarray(l)
            if trace_up: #if we are going upstream from terminus
                x_n = float(currentpt[0] - (vx/vm)*dx)
                y_n = float(currentpt[1] - (vy/vm)*dx)
                if np.linalg.norm(r-current_r)>side_tolerance: #check if the width is going to be something wacky, and replace with parallel projection if so
                    #DO SOMETHING ABOUT REPLACEMENT VALUE
                    x_r = float(current_r[0] - (vx/vm)*dx)
                    y_r = float(current_r[1] - (vy/vm)*dx)
                    r = np.array((x_r, y_r))
                if np.linalg.norm(l-current_l)>side_tolerance:
                    x_l = float(current_l[0] - (vx/vm)*dx)
                    y_l = float(current_l[1]-(vy/vm)*dx)
                    l = np.array((x_l, y_l))
            else:
                x_n = float(currentpt[0] + (vx/vm)*dx)
                y_n = float(currentpt[1] + (vy/vm)*dx)
                if np.linalg.norm(r-current_r)>side_tolerance:#check if the width is going to be something wacky, and replace with parallel projection if so
                    x_r = float(current_r[0] + (vx/vm)*dx)
                    y_r = float(current_r[1] + (vy/vm)*dx)
                    r = np.array((x_r, y_r))
                if np.linalg.norm(l-current_l)>side_tolerance:
                    x_l = float(current_l[0] + (vx/vm)*dx)
                    y_l = float(current_l[1] + (vy/vm)*dx)
                    l = np.array((x_l, y_l))
            nextpt = (x_n, y_n)
            w = np.linalg.norm(r-l)
                     
            currentpt = nextpt
            current_r, current_l = r, l
            outarr.append((currentpt[0], currentpt[1], w)) 
            rightside.append((current_r))
            leftside.append((current_l)) 
            nstep += 1
        
    outarr = np.asarray(outarr)
    rightside = np.asarray(rightside)
    leftside = np.asarray(leftside)
    if output_sidecoords:
        return outarr, rightside, leftside
    else:
        return outarr


def WriteNetwork(startcoords, trace_up=False, output_name='glacier-network-lines.csv'):
    """Given a list of starting coordinates, applies Trace_wWidth sequentially and writes output to a CSV file.
    """
    outdict = {}
    for j in range(len(startcoords)):
        k = startcoords[j]
        line_coords = Trace_wWidth(k[0], k[1], trace_up)
        outdict[j] = (line_coords)
        
    with open(output_name, 'wb') as csvfile:
        linewriter = csv.writer(csvfile, delimiter=',')
        linewriter.writerow(['Line-number', 'x-coord', 'y-coord', 'width'])
        for n in range(len(outdict)):
            for m in range(len(outdict[n])):
                linewriter.writerow([str(n), outdict[n][m][0], outdict[n][m][1], outdict[n][m][2]])
        