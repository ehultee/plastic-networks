# Rewriting plastic model to add SMB options
# Model functions will now be capitalised for more legible code
#28 Nov 16 EHU
from numpy import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from types import *

# Define Dimensional and Dimensionless parameters
H0=1e3 #characteristic height for nondimensionalisation 
L0=10e3 #characteristic length (10km)
#tau_yield = 200e3 #actual yield strength, 100kPa?
g = 9.8 
rho_ice = 920.0 #ice density kg/m^3
rho_sea=1020.0 #seawater density kg/m^3

#Bingham stress functions-----------------
def B_var(elev, thick, pos=None, time=None): #variation by altitude and ice thickness (effective pressure at base)...pos, time arguments required by plasticmodel
    if elev<0:
        D = -elev #Water depth D the nondim bed topography value when Z<0
    else:
        D = 0
    N = rho_ice*g*thick*H0 - rho_sea*g*D*H0
    mu = 0.01 #Coefficient between 0 and 1
    tau_y = tau_val + mu*N
    return tau_y/(rho_ice*g*H0**2/L0) 
    
def B_const(elev, thick, pos=None, time=None):  #functional form of B if using a constant yield strength
    return tau_val/(rho_ice*g*H0**2/L0)

#------------------------------
def BalanceThick(bed,B):
    """Water balance ice thickness.
    Arguments:
        bed is nondim bed value at point of interest
        B is nondim yield strength
    Returns nondim ice thickness for water balance.
    """
    if bed<0:
        D = -1*bed
    else:
        D = 0
    return (2*B*H0/L0) + math.sqrt((rho_sea*(D**2)/rho_ice)+(H0*B/L0)**2)

#------------------------------
def FlotationThick(bed):
    """Minimum ice thickness before flotation.
    Argument:
        bed is nondim bed value at point of interest
    """
    if bed<0:
        D = -bed
    else:
        D = 0
    return (rho_sea/rho_ice)*D

#------------------------------
def FlowlineArr(filename):
    """Reads in a list of points on a flowline from a .txt file and makes an array of them for use in modelling.
    Argument: filename a string indicating where the list of points is stored.
    Returns an array of shape (n, 2) from a list of n points.
    """
    linraw = filename
    linarr = []
    lines = open(linraw, 'r').read().split('\n')
    for n in xrange(0, len(lines), 1):
        linstrip = lines[n].strip('( ),')
        splitup = linstrip.split(',')
        xval = float(splitup[0])
        yval = float(splitup[1])
        linarr.append((xval, yval))
    linarr = np.array(linarr)
    return linarr

#------------------------------
def ds(num, flarr):
    """Arclength utility finding path length between adjacent points of an array. To be used in for-loop.
    Arguments: num is index of array element (calculates arclength between element num and element num-1)
    flarr is array of interest
    Returns float of the distance."""
    if num==0:
        return 0
    else:
        return math.sqrt((flarr[num,0]-flarr[num-1,0])**2 + (flarr[num,1]-flarr[num-1,1])**2)/L0

#------------------------------
def ArcArray(flowline):
    """Takes a flowline of n points (xi, yi) and returns array of n corresponding arc length values (floats)"""
    arcarr = []
    xs = xrange(0,len(flowline),1) # Just a counter for the for-loop
    for n in xs:
        if n==0:
            arcarr.append(0)
        else:
            arcarr.append(arcarr[-1]+ds(n, flowline))
    return arcarr

#------------------------------
def Evenspace(points, spacing):
    """Evenspace gives an evenly spaced flowline with the specified spacing, given a set of points along the line.
    Arguments:
        points should be a list of points [(x1, y1), ..., (xn, yn)] along the line.
        spacing should be spacing (in m) desired for the sampling.
    
    Returns an array of evenly spaced points (xi, yi).
    """
    pointsarray = np.array(points)
    arc = ArcArray(pointsarray)
    
    #Now making functions to express flowline position w.r.t. arclength, and come back out with even spacing
    spacex = interpolate.interp1d(arc, pointsarray[:,0], 'linear', copy=True)
    spacey = interpolate.interp1d(arc, pointsarray[:,1], 'linear', copy=True)
    nstep = floor((arc[-1]*L0)/spacing)
    sample = linspace(0, arc[-1], nstep)
    evenarr=[]
    for j in xrange(0, len(sample), 1):
        evenx=float(spacex(sample[j]))
        eveny=float(spacey(sample[j]))
        evenarr.append((evenx, eveny))
    evenarr=np.array(evenarr)
    return evenarr

#------------------------------
def FlowProcess(flowline, interpolatedquantity):
    """Flowbed samples an observed quantity along a flowline and makes a function for use in plasticmodel.  Useful for bed, ice thickness, surface elevation, etc.
    Arguments:
        flowline should be an array of points (xi, yi) defining the flowline.  Denser samplings of points are better.
        interpolatedbed should be the 2D interpolation of the desired values from data, on the same grid as flowline.
    
    Returns a function to be evaluated in terms of arclength along the flowline. (so the return takes an argument of arclength)
    """
    vals = []
    arc = ArcArray(flowline)
    xs = xrange(0,len(flowline),1) # Just a counter for the for-loop
    for n in xs:
        vals.append(interpolatedquantity(flowline[n,0], flowline[n,1]))
    
    sqzd = np.squeeze(vals) #need to alter shape to get into 1d interpolation
    funct=interpolate.interp1d(arc, sqzd, 'linear', copy=True)
    return funct

#------------------------------
def PlasticProfile(bedfunction, Bfunction, startpoint, hinit, endpoint, Npoints, obsheightfunction):
    """Arguments used:
        bedfunction should be function of arclength returning bed elevation of the glacier.
        Bfunction is nondim yield strength.  Should be function with arguments elevation, ice thickness, position, and time (can just not use last two if no variation)
        Startpoint is where (in arclength space) integration should start.
        hinit is nondim initial height.  Could be given by water balance, obs, or some thinning from reference height.
        Endpoint is where integration should stop.
        Npoints is how many model points to use (suggested 25000+)
        #Resolution (in m) is how closely we want to space the model sample points (CURRENTLY USING NPOINTS INSTEAD OF RESOLUTION)
        Obsheightfunction is the observations for comparison.  May need to process from a thickness measurement. (working on functionality to make this argument optional)
     
     PlasticModel function returns arrays: 
        'horiz' which gives the arclength values of each modelled point, 
        'SEarr' which gives modelled surface elevation at each point, 
        'basearr' which gives bed values at each point,
        'obsarr' which gives values of observed surface elevation at each point.
        Note if running several times over same bed or modelling same retreat with only 1 set of obs, probably will not need to use basearr and obsarr
     """
    #N = ceil(abs((endpoint-startpoint)*L0/resolution))
    N = Npoints
    horiz = linspace(startpoint, endpoint, N)
    dx = mean(diff(horiz))

    if dx<0:
        print 'Detected: running from upglacier down to terminus.'
    elif dx>0:
        print 'Detected: running from terminus upstream.'
    
    SEarr = []
    thickarr = []
    basearr = []
    obsarr = []
    
    SEarr.append(hinit)
    thickarr.append(hinit-(bedfunction(startpoint)/H0))
    basearr.append(bedfunction(startpoint)/H0)
    obsarr.append((obsheightfunction(startpoint))/H0)
    #if abs(hinit-obsarr[-1])>10/H0:
        #print 'Initial value does not match obs at starting point.  If you are modelling with e.g. prescribed thinning from a reference point, RMS and CV(RMSE) values may be less relevant.'
    # Need to find a better way to do this warning
    for x in horiz[1::]:
        bed = bedfunction(x)/H0  # value of interpolated bed function
        obsheight = (obsheightfunction(x))/H0
        modelthick = thickarr[-1]
        B = Bfunction(tau_val, bed, modelthick, None, None)
        #Break statements for thinning below yield, water balance, or flotation
        if dx<0:
            if modelthick<balancethick(bed,B):
                print 'Thinned below water balance at x=' + str(10*x)+'km'
                break
        if modelthick<flotationthick(bed):
            print 'Thinned below flotation at x=' + str(10*x)+'km'
            break
        if modelthick<4*B*H0/L0:
            print 'Thinned below yield at x=' +str(10*x)+'km'
            break
        else:
            basearr.append(bed)
            SEarr.append(SEarr[-1]+(B/modelthick)*dx) 
            thickarr.append(SEarr[-1]-basearr[-1])
            obsarr.append(obsheight)
    
    error = np.sqrt(((np.array(SEarr)-np.array(obsarr))**2).mean())
    CVrms = error/mean(SEarr)
    print 'RMS error: '+ str(error) +',   CV(RMSE): ' + str(CVrms)
    return (horiz[0:len(SEarr)], SEarr, basearr, obsarr)

##------------------------------
#def PlasticEvol(glaciername, bedfunction, Bfunction, startpoint, hinit, endpoint, Npoints, obsheightfunction):
#        """Arguments used:
#        bedfunction should be function of arclength returning bed elevation of the glacier.
#        Bfunction is nondim yield strength.  Should be function with arguments elevation, ice thickness, position, and time (can just not use last two if no variation)
#        Startpoint is where (in arclength space) integration should start.
#        hinit is nondim initial height.  Could be given by water balance, obs, or some thinning from reference height.
#        Endpoint is where integration should stop.
#        Npoints is how many model points to use (suggested 25000+)
#        #Resolution (in m) is how closely we want to space the model sample points (CURRENTLY USING NPOINTS INSTEAD OF RESOLUTION)
#        Obsheightfunction is the observations for comparison.  May need to process from a thickness measurement. (working on functionality to make this argument optional)
#     
#     PlasticEvol runs PlasticModel several times for  
#        'horiz' which gives the arclength values of each modelled point, 
#        'SEarr' which gives modelled surface elevation at each point, 
#        'basearr' which gives bed values at each point,
#        'obsarr' which gives values of observed surface elevation at each point.
#        Note if running several times over same bed or modelling same retreat with only 1 set of obs, probably will not need to use basearr and obsarr
#     """