# Formalizing flowline selection algorithm for Greenland glaciers
# Adding functionality to find widths
# 16 March 2018  EHU
## Edits 20 Sept 2018: adding automated identification of centerline

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
import shapely.geometry as geom
from shapely.ops import nearest_points
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter

##------------------------------
## Data used for these functions
##------------------------------

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
def check_sides(x, y, Vx = vf_x, Vy = vf_y, V=vf, side_cutoff = 1.0, dw=20):
    """Reports coordinates of glacier 'sides' by searching velocity of points normal to direction of flow and comparing to a cutoff value. 
    x: x-coordinate of current point along flowline
    y: y-coordinate of current point along flowline
    Vx: 2D interpolated function of the x-velocity field
    Vy: 2D interpolated function of the y-velocity field
    V: 2d interpolated function of magnitude of surface velocity
    side_cutoff: the minimum velocity, in units of V field, for ice to be included in branch.  Default is 1 m/d.
    dw: spatial step size, in meters, for checking across flow.  Default is 20 m.
    """
    vx = Vx(x, y)
    vy = Vy(x, y)
    vm = V(x, y)
    vc = (vx, vy) #centerline velocity vector for testing direction of flow
    
    nr = 1/vm*np.array([-vy, vx]) #right normal
    nl = 1/vm*np.array([vy, -vx]) #left normal
    
    
    Cr = 1 # counter for width increments on right
    Cl = 1 # counter for width increments on left
    speed_rightside = vm #initial value of side speed
    speed_leftside = vm #initial value of side speed
    vr = (Vx(x,y), Vy(x,y)) #initial vr = vc
    vl = (Vx(x,y), Vy(x,y)) #initial vl = vc
        

    while np.vdot(vr, vc)/(speed_rightside*vm) > 0.5 and speed_rightside > side_cutoff:
        x_r = float(x + Cr*dw*nr[0])
        y_r = float(y + Cr*dw*nr[1])
        vr = (Vx(x_r, y_r), Vy(x_r, y_r))
        #speed_rightside = V(x_r, y_r)
        speed_rightside = np.linalg.norm(vr)
        #print Cr
        Cr += 1
        
    while np.vdot(vl, vc)/(speed_leftside*vm) > 0.5 and speed_leftside > side_cutoff:
        x_l = float(x + Cl*dw*nl[0])
        y_l = float(y + Cl*dw*nl[1])
        vl = (Vx(x_l, y_l), Vy(x_l, y_l))
        #speed_leftside = V(x_l, y_l)
        speed_leftside = np.linalg.norm(vl)
        #print Cl
        Cl += 1
        
    #print 'Right side: step coefficient {}, speed {}, parallel velocity {}'.format(Cr, speed_rightside, np.vdot(vr, vc)/np.linalg.norm(vc))
    #print 'Left side: step coefficient {}, speed {}, parallel velocity {}'.format(Cl, speed_leftside, np.vdot(vl, vc)/np.linalg.norm(vc))
    try:
        rightcoords = (x_r, y_r)
    except UnboundLocalError: #will raise this error if first step to right side hit cutoff and x_r, y_r not assigned
        rightcoords = (x, y)
   
    try:
        leftcoords = (x_l, y_l)
    except UnboundLocalError: #will raise this error if first step to left side hit cutoff and x_l, y_l not assigned
        leftcoords = (x, y)
    
    return rightcoords, leftcoords


def Trace_wWidth(startcoord_x, startcoord_y, trace_up=True, xarr=X, yarr=Y, Vx = vf_x, Vy = vf_y, V = vf, dx=150, side_cutoff = 0.5, dw=20, output_sidecoords=False, smooth_width=True):
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
    output_sidecoords: determines whether to return raw coordinates of left and right side of glacier
    """
    outarr = []
    rightside = []
    leftside = []
    raw_width = []
    
    #side_tolerance = 3*dx
    
    
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
                #if np.linalg.norm(r-current_r)>side_tolerance: #check if the width is going to be something wacky, and ignore it if so
                #    #DO SOMETHING ABOUT REPLACEMENT VALUE
                #    #x_r = float(current_r[0] - (Vx(r[0], r[1])/V(r[0], r[1]))*side_tolerance)
                #    #y_r = float(current_r[1] - (Vy(r[0], r[1])/V(r[0], r[1]))*side_tolerance)
                #    #r = np.array((x_r, y_r))
                #    r = np.array((np.nan, np.nan))
                #if np.linalg.norm(l-current_l)>side_tolerance:
                #    #x_l = float(current_l[0] - (Vx(l[0], l[1])/V(l[0], l[1]))*side_tolerance)
                #    #y_l = float(current_l[1]-(Vy(l[0], l[1])/V(l[0], l[1]))*side_tolerance)
                #    #l = np.array((x_l, y_l))
                #    l = np.array((np.nan, np.nan))
            else:
                x_n = float(currentpt[0] + (vx/vm)*dx)
                y_n = float(currentpt[1] + (vy/vm)*dx)
                #if np.linalg.norm(r-current_r)>side_tolerance:#check if the width is going to be something wacky, and replace with parallel projection if so
                #    x_r = float(current_r[0] + (Vx(r[0], r[1])/V(r[0], r[1]))*side_tolerance)
                #    y_r = float(current_r[1] + (Vy(r[0], r[1])/V(r[0], r[1]))*side_tolerance)
                #    r = np.array((x_r, y_r))
                #if np.linalg.norm(l-current_l)>side_tolerance:
                #    x_l = float(current_l[0] + (Vx(l[0], l[1])/V(l[0], l[1]))*side_tolerance)
                #    y_l = float(current_l[1] + (Vy(l[0], l[1])/V(l[0], l[1]))*side_tolerance)
                #    l = np.array((x_l, y_l))
            nextpt = (x_n, y_n)
            w = np.linalg.norm(r-l)
                     
            currentpt = nextpt
            current_r, current_l = r, l
            outarr.append((currentpt[0], currentpt[1])) 
            rightside.append((current_r))
            leftside.append((current_l))
            raw_width.append((w)) 
            nstep += 1
        
    outarr = np.asarray(outarr)
    rightside = np.asarray(rightside)
    leftside = np.asarray(leftside)
    raw_width = np.asarray(raw_width)
    
    if smooth_width:
        width = savgol_filter(raw_width, 21, 4, mode='mirror')
    else:
        width = raw_width
        
    if output_sidecoords:
        return outarr, width, rightside, leftside
    else:
        return outarr, width

def IDCentralBranch(lines):
    """Given a dictionary of flowline-like objects, oriented from terminus to head, finds which one is 'centerline' (terminal point closest to centroid of terminus) and returns its dictionary key
    """
    if len(lines.keys())>1:
        end_points = [lines[k][0][0:2] for k in lines.keys()] #(x,y) coordinates at downstream end of each line
        end_ls = geom.LineString(end_points) #need LineString to identify centroid
        end_mp = geom.MultiPoint(end_points) #need MultiPoint to find closest point to centroid
        end_centr = end_ls.centroid
        central_point = nearest_points(end_mp, end_centr)[0] #(x,y) coords of point closest to centroid
        ## Now check which line has that terminus point
        for i,p in enumerate(end_mp):
            if p.distance(central_point) > 0.01: #allowing tolerance for floating-point errors in Shapely
                pass
            elif p.distance(central_point) < 0.01:
                print 'Central branch has key {}'.format(lines.keys()[i])
                central_key = lines.keys()[i] # be careful that lines.keys() may be ordered differently from what we expect--but should match end_mp by construction
            else:
                print 'Error finding central key'
    elif len(lines.keys())==1:
        central_key = lines.keys()[0]
    else:
        print 'Error finding central key: no valid lines in input'
    
    return central_key
    
        
        

def FilterMainTributaries(lines, Vx = vf_x, Vy = vf_y):
    """Given a dictionary of flowline-like objects, oriented from terminus to head, finds where they should be joined into mainline-tributary network.
    Expects input in the form of {'linenumber': ((x1, y1, w1), (x2, y2, w2))} - coordinates with width
    """
    branches = {}
    central_branch_key = IDCentralBranch(lines)
    branches[central_branch_key] = lines[central_branch_key] #Keep full length.  Note central line will need to be re-indexed to have index 0 for other functions in PlasticNetwork
    
    for ln in lines.keys():
        #print ln
        if ln==central_branch_key: #central line is reference, doesn't need filtering
            pass
        else:
            mainline = branches[central_branch_key] #should be an array of (x,y, w) coords describing flowline
            full_line = lines[ln]
            comparison_length = min(len(mainline), len(full_line))
            
            #Check directionality of velocity
            Vmain = [(Vx(mainline[j][0], mainline[j][1]), Vy(mainline[j][0], mainline[j][1])) for j in range(comparison_length)]
            Vbranch = [(Vx(full_line[j][0], full_line[j][1]), Vy(full_line[j][0], full_line[j][1])) for j in range(comparison_length)]
            
            for j in range(comparison_length):
                dotprod = np.vdot(Vmain[j], Vbranch[j])
                #print dotprod
                mag = np.linalg.norm(Vmain[j])*np.linalg.norm(Vbranch[j])
                #print mag
                cos_between = dotprod/mag
                if cos_between < 0.95: #if angle between velocity vectors exceeds ~18 degrees
                    trunc_index = j
                    print 'Truncation index of line {} is {}.  Exiting.'.format(ln, trunc_index)
                    break
                else:
                    continue
            try:
                truncated_line = full_line[trunc_index::]
            except NameError: #no truncation index found
                print 'Line {} may parallel main line for full length. Removing from set.'.format(ln)
                truncated_line = None
            if trunc_index==0:
                print 'Line {} may be dynamically independent from centerline at terminus.  Removing from set.'.format(ln)
                truncated_line = None
            branches[ln] = truncated_line
                
            #diff = [np.linalg.norm(np.asarray(full_line[j][0:2])-np.asarray(mainline[j][0:2])) for j in range(comparison_length)] #array of pointwise distances.  Needs [0:2] so that width is not included in norm
            #width = [mainline[j][2] for j in range(comparison_length)]
            #comparison_width = 0.5*np.array(width)
            #if np.all(diff < comparison_width): #check whether lines always within main branch width of each other - if so probably same branch
            #    print 'Line {} may parallel main line for full length.'.format(ln)
            #    branches[ln] = lines[ln]
            #else:
            #    trunc_index = np.argwhere(diff > comparison_width)[0]
            #    print trunc_index
            #    truncated_line = full_line[trunc_index::]
            #    branches[ln] = truncated_line
      
    branches_NoNone = {key: val for key,val in branches.iteritems() if val is not None} # dictionary with None entries removed  
    remaining_keys = branches_NoNone.keys()
    remaining_keys.remove(central_branch_key)
    branches_cleaned = {int(k+1): branches_NoNone[remaining_keys[k]] for k in range(len(remaining_keys))} #making output dictionary have sequential keys, starting from 1 so central branch can be 0 index
    branches_cleaned[0] = branches_NoNone[central_branch_key]
    
    #return branches_NoNone
    return branches_cleaned
            
    

def WriteNetwork(startcoords, trace_up=False, output_name='glacier-network-lines.csv'):
    """Given a list of starting coordinates, applies Trace_wWidth sequentially and writes output to a CSV file.
    """
    lines = {}
    for j in range(len(startcoords)):
        k = startcoords[j]
        line_coords, width = Trace_wWidth(k[0], k[1], trace_up=trace_up)
        if sum(np.isnan(line_coords))>0:
            pass #Don't save the line if it contains nan points
        else:
            xyw = [(line_coords[n][0], line_coords[n][1], width[n]) for n in range(len(line_coords))]
            lines[j] = (xyw)

   
    if len(lines) > 1:    
        outdict = FilterMainTributaries(lines)
    else:
        outdict = lines
        
    with open(output_name, 'wb') as csvfile:
        linewriter = csv.writer(csvfile, delimiter=',')
        linewriter.writerow(['Line-number', 'x-coord', 'y-coord', 'width'])
        for n in range(len(outdict)):
            for m in range(len(outdict[n])):
                linewriter.writerow([str(n), outdict[n][m][0], outdict[n][m][1], outdict[n][m][2]])
        