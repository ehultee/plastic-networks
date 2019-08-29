#Read-in utilities to help process Huss and Kienholz data
from numpy import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import collections
from types import *

# Define Dimensional and Dimensionless parameters
H0=1e3 #characteristic height for nondimensionalisation 
L0=10e3 #characteristic length (10km)
#tau_yield = 200e3 #actual yield strength, 100kPa?
g = 9.8 
rho_ice = 920.0 #ice density kg/m^3
rho_sea=1020.0 #seawater density kg/m^3

#B = tau_yield/(rho_ice*g*H0**2/L0) #nondimensional Bingham stress

def beddata(DEM, thicknessmap, bedsmoothing=300):
    """Function to read in datasets provided by Huss and solve for bed topography
    
    Inputs: 'DEM' and 'thicknessmap' should be filenames, 
    'bedsmoothing' should be the radius (in m) of smoothing to apply to bed.  Default is 300m
    
    Outputs: xx, yy, zz arrays; 
    bedinterp, thickinterp, seinterp 2d functions;
    xbounds, ybounds values
    """
    
    SE = []
    H = []
    Z = []
    
    hdrs_s = []
    # Read and interpolate DEM
    with open(DEM, 'r') as f:
        #headarray = [ncols, nrows, xll, yll, cellsize, nodata_val]
        print 'Currently reading DEM'
        for num in xrange(0, 6, 1):
            headern = f.readline().split() #read info in header
            headerval = float(headern[1])
            hdrs_s.append(headerval)
        ncols_s = int(hdrs_s[0])
        nrows_s = int(hdrs_s[1])
        xll_s = hdrs_s[2]
        yll_s = hdrs_s[3]
        delta_s = hdrs_s[4]
        
        lines = f.read().split('\r\n') #lines in Huss files end with \r\n rather than \n alone
        for n in xrange(0, (len(lines)-1), 1): #leconte_dem file reads with empty line at end, need to ignore
            line = lines[n]
            splitup = line.split() #no splitter defined--will automatically see whitespace
            for ent in xrange(0, len(splitup), 1): #every entry of DEM is surface elevation, first entry is tabbed
                SEval = float(splitup[ent])
                SE.append(SEval)
    SE = np.array(SE)
    se = np.reshape(SE, (-1, ncols_s))
    
    xs = []
    ys = []
    for i in xrange(0, ncols_s, 1): #x values change by column
        xs.append(xll_s + i*delta_s)
    for j in xrange(0, nrows_s, 1): #y values change by row
        ys.append(yll_s + j*delta_s) 
    #ys = ys[::-1]
    seinterp = interpolate.interp2d(xs,ys,se, kind='linear', copy=True, bounds_error=False) #interpolated DEM
    
    
    # Read and interpolate thickness map
    hdrs_h = []
    with open(thicknessmap, 'r') as f:
        print 'Currently reading ice thickness map'
        for num in xrange(0, 6, 1):
            headern = f.readline().split() #read info in header
            headerval = float(headern[1])
            hdrs_h.append(headerval)
        ncols_h = int(hdrs_h[0])
        nrows_h = int(hdrs_h[1])
        xll_h = hdrs_h[2]
        if xll_h != xll_s:
            print 'Warning: DEM and thickness map grids start in different locations'
        yll_h = hdrs_h[3]
        delta_h = hdrs_h[4]
    
        lines = f.read().split('\r\n') #lines in Huss files end with \r\n rather than \n alone
        for n in xrange(0, (len(lines)-1), 1): #leconte_thick file reads with empty line at end, need to ignore
            line = lines[n]
            splitup = line.split() #no splitter defined--will automatically see whitespace
            for ent in xrange(0, len(splitup), 1): 
                Hval = float(splitup[ent])
                H.append(Hval)
    H = np.array(H)
    hh = np.reshape(H, (-1, ncols_h))
    
    xh= []
    yh = []
    for i in xrange(0, ncols_h, 1): #x values change by column
        xh.append(xll_h + i*delta_h)
    for j in xrange(0, nrows_h, 1): #y values change by row
        yh.append(yll_h + j*delta_h) 
    #xh = xh[::-1]
    #yh = yh[::-1]
    thickinterp = interpolate.interp2d(xh,yh,hh, kind='linear', copy=True, bounds_error=False)
    
    
    #Find appropriate grid for interpolation; get bed
    xx = []
    yy = []
    ncols = max(ncols_s, ncols_h)
    nrows = max(nrows_s, nrows_h)
    if ncols == ncols_s:
        if nrows == nrows_s:
            print 'Finer grid is DEM'
            delta = delta_s
            xx = np.array(xs)
            yy = np.array(ys)
            yy = yy[::-1]
            hnew = thickinterp(xs, ys)
            zz = np.array(se) - np.array(hnew)
        else:
            print 'Error 1: grid sizes too different'
    if ncols == ncols_h:
        if nrows == nrows_h:
            print 'Finer grid is ice thickness'
            delta = delta_h
            xx = np.array(xh)
            yy = np.array(yh)
            #yy = yy[::-1]
            senew = seinterp(xh, yh)
            zz = np.array(senew) - np.array(hh)
        else:
            print 'Error 2: grid sizes too different'
    
    #print 'Reshaping y.  REMEMBER to check output graphically!'
    #yy = yy[::-1] #not sure why this should be...surely for positive yll, y should be increasing?  This makes the picture of Chenega correct...may need to check again
    
    #Smoothing bed
    if bedsmoothing == None:    
        unsmoothZ = zz #just return zz unsmoothed
    else:
        dsmooth = ceil(bedsmoothing/delta) #smoothing
        zz = gaussian_filter(zz, dsmooth) #smoothing the bed
    bedinterp = interpolate.interp2d(xx,yy,zz, kind='linear', copy=True, bounds_error=True) #might need to do transpose of Z as in Columbia?    
    
    #finding bounds
    xbounds = (min(xx), max(xx))
    ybounds = (min(yy), max(yy))
    
    return(xx, yy, zz, bedinterp, thickinterp, seinterp, xbounds, ybounds)
#------------------------------

def kienholzlines(filename, nlines, xbounds, ybounds):
    """Function to read in flowlines as sent by C. Kienholz and processed to add vertex points in ArcMap.
    Input: filename; nlines=number of flowlines; xbounds, ybounds of relevant DEM
    Output: list of flowlines
    """
    
    f = open(filename,'r')
    
    header = f.readline() #header line
    keys = header.split(',') #get names of variables
    
    data = {k : [] for k in keys[:-1]} #end of line has hidden characters, so 'point_m' does not get read
    data['POINT_M'] = []
    data['Length_ID'] = collections.OrderedDict() #new dictionary that counts how many points (i.e. lines of file) are in each flowline.  Must be ordered for later iteration!
    data['Lineslist'] = [[] for k in range(nlines)] 
    
    lines = f.readlines()
    f.close()
    
    temp = []
    j = 0
    for i,l in enumerate(lines):
        parts = l.split(',')
        data['FID'].append(parts[0])
        data['OBJECTID'].append(parts[1])
        data['GLIMSID'].append(parts[2])
        data['LENGTH_BR'].append(parts[3])
        data['MAIN'].append(parts[4])
        data['ORIG_FID'].append(parts[5])
        data['POINT_X'].append(parts[6])
        data['POINT_Y'].append(parts[7])
        data['POINT_M'].append(parts[8])
        
        x_coord = float(parts[6])
        y_coord = float(parts[7])
    
        if parts[1] not in data['Length_ID'].keys(): #finding out where lines separate 
            temp = []
            data['Length_ID'][parts[1]] = 1
            j+=1
        else:
            data['Length_ID'][parts[1]] += 1
        if xbounds[0]<x_coord<xbounds[1]: #taking out values outside of map area
            if ybounds[0]<y_coord<ybounds[1]:            
                temp.append((x_coord, y_coord))
        data['Lineslist'][j-1] = np.array(temp) 
       
    centrelines_list = [np.array(data['Lineslist'][j])[::-1] for j in range(nlines)] #making arrays, reversed to start at terminus rather than peak
    
    return centrelines_list
      