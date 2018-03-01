## Utilities file for Greenland modelling functions
# 30 Nov 2017  EHU

from numpy import *
#from netCDF4 import Dataset
import numpy as np
#from scipy import interpolate
from scipy import spatial
#from scipy.ndimage import gaussian_filter
from plastic_utilities_v2 import*
import collections

##-------------------------------
##  FINDING AND WRITING LINES
##-------------------------------

def Continuous_DownVStep(startcoord_x, startcoord_y, surface, xarr, yarr, Vx, Vy):
    """Traces flowlines down from an upstream point.
    surface: gridded surface
    xarr: gridded x-values for area of interest
    yarr: gridded y-vals
    Vx: 2d interpolated function for x-component of velocity
    Vy: ditto y-component
    """
    outarr = []
    
    currentpt = (startcoord_x, startcoord_y)
    nstep = 0
    
    while nstep<1000:
        #dt = 10 #days, we integrate velocity presented in m/day to follow path up day by day
        dx = 150 #m, normalizing our step size
        vmax = 50 # m/day, the maximum velocity we believe
        
        vx = Vx(currentpt[0], currentpt[1])
        vy = Vy(currentpt[0], currentpt[1])
        vm = np.linalg.norm((vx, vy))
        
        if vx > vmax:
            print 'X-velocity exceeds maximum recognised.  Exiting step routine.'
            break
        elif vy > vmax:
            print 'Y-velocity exceeds maximum recognised.  Exiting step routine.'
            break
        else:
            x_n = float(currentpt[0] + (vx/vm)*dx)
            y_n = float(currentpt[1] + (vy/vm)*dx)
            nextpt = (x_n, y_n)
            print nextpt
            currentpt = nextpt
            outarr.append((currentpt[0], currentpt[1]))  
            nstep += 1
        
    outarr = np.asarray(outarr)

    return outarr

#def WriteNetworkDict()

##-------------------------------
##  READING IN SAVED LINES
##-------------------------------

## Read in line from CSV
def Flowline_CSV(filename, nlines):
    """Function to read in flowlines in CSV format, similar to those sent by C. Kienholz for Alaska.
    Input: filename; nlines=number of flowlines
    Output: list of flowlines
    """
    
    f = open(filename,'r')
    
    header = f.readline() #header line
    hdr = header.strip('\r\n')
    keys = hdr.split(',') #get names of variables
    #keys[-1] = keys[-1].strip('\r\n')
    
    data = {k : [] for k in keys} #end of line has hidden characters, so 'point_m' does not get read
    #data['Line number'] = []
    data['Length_ID'] = collections.OrderedDict() #new dictionary that counts how many points (i.e. lines of file) are in each flowline.  Must be ordered for later iteration!
    data['Lineslist'] = [[] for k in range(nlines)] 
    
    lines = f.readlines()
    f.close()
    
    temp = []
    j = 0
    for i,l in enumerate(lines):
        linstrip = l.strip('\r\n')
        parts = linstrip.split(',')
        
        #data['Line-number'].append(parts[0])
        #data['x-coord'].append(parts[1])
        #data['y-coord'].append(parts[2])
        
        x_coord = float(parts[1])
        y_coord = float(parts[2])
    
        if parts[0] not in data['Length_ID'].keys(): #finding out where lines separate 
            temp = []
            data['Length_ID'][parts[0]] = 1
            j+=1 
        else:
            data['Length_ID'][parts[0]] += 1
        #if xbounds[0]<x_coord<xbounds[1]: #taking out values outside of map area
        #    if ybounds[0]<y_coord<ybounds[1]:   
        
        temp.append((x_coord, y_coord))
        data['Lineslist'][j-1] = np.array(temp) 

        #data['Lineslist'][j] = np.array(temp) 

              
    centrelines_list = [np.array(data['Lineslist'][j])[::-1] for j in range(nlines)] #making arrays, reversed to start at terminus rather than peak
    
    return centrelines_list


def TopToTerm(branchlist):
    mainline = branchlist[0]
    maintree = spatial.KDTree(mainline)
    
    full_lines = {}
    j = 0
    
    while j<len(branchlist):
        branch = branchlist[j]
        pt = branch[0]
        dist, idx = maintree.query(pt, distance_upper_bound=5000)  #find distances and indices of nearest neighbours along main line
        #print dist, idx
        if idx==len(mainline): #if branch does not intersect with the main line
            print 'Branch {} does not intersect main line.  Searching nearest trib.'.format(j)
            tribtree = spatial.KDTree(full_lines[j-1]) #line of nearest trib
            dist_t, idx_t = tribtree.query(pt, distance_upper_bound=1000)
            if idx==len(full_lines[j-1]):
                print 'Branch {} also does not intersect tributary {}.  Appending raw line.  Use with caution.'.format(j, j-1)
                full_lines[j] = branch
            else:
                tribfrag = branchlist[j-1][:idx_t]
                fullbranch = np.concatenate((tribfrag, branch))
                full_lines[j] = fullbranch
            j+=1
        else:
            print mainline[idx]
            mainfrag = mainline[:idx]
            fullbranch = np.concatenate((mainfrag, branch))
            full_lines[j] = fullbranch
            j+=1
        
    return full_lines

##-------------------------------
##  YIELD STRENGTH OPTIMISATION
##-------------------------------

# Dimensional and Dimensionless parameters
H0=1e3 #characteristic height for nondimensionalisation 
L0=10e3 #characteristic length (10km)
#tau_yield = 100e3 #initial guess to initialize
#tau_0 = 100e3
g = 9.8 
rho_ice = 920.0 #ice density kg/m^3
rho_sea=1020.0 #seawater density kg/m^3
#Bingham stress function
def B_var(tau_0, elev, thick, pos=None, time=None): #variation by altitude and ice thickness (effective pressure at base)...pos, time arguments required by plasticmodel
    if elev<0:
        D = -elev #Water depth D the nondim bed topography value when Z<0
    else:
        D = 0
    N = rho_ice*g*thick*H0 - rho_sea*g*D*H0
    mu = 0.01 #Coefficient between 0 and 1
    tau_y = tau_0 + mu*N
    return tau_y/(rho_ice*g*H0**2/L0) 
    
def B_const(tau_yield, elev, thick, pos=None, time=None):  #functional form of B if using a constant yield strength
    return tau_yield/(rho_ice*g*H0**2/L0)


def plasticmodel_error(bedfunction, tau_val, Bfunction, startpoint, hinit, endpoint, Npoints, obsheightfunction):
    """Arguments used:
        bedfunction should be function of arclength returning bed elevation of the glacier.
        Bfunction is nondim yield strength.  Should be function with arguments elevation, ice thickness, position, and time (can just not use last two if no variation)
        Startpoint is where (in arclength space) integration should start.
        hinit is nondim initial height.  Could be given by water balance, obs, or some thinning from reference height.
        Endpoint is where integration should stop.
        Npoints is how many model points to use (suggested 25000+)
        #Resolution (in m) is how closely we want to space the model sample points (CURRENTLY USING NPOINTS INSTEAD OF RESOLUTION)
        Obsheightfunction is the observations for comparison.  May need to process from a thickness measurement. (working on functionality to make this argument optional)
     
     plasticmodel_error function returns values: 
        RMS error
        CV(RMSE)
     """
    #N = ceil(abs((endpoint-startpoint)*L0/resolution))
    N = Npoints
    horiz = linspace(startpoint, endpoint, N)
    dx = mean(diff(horiz))

    #if dx<0:
    #    print 'Detected: running from upglacier down to terminus.'
    #elif dx>0:
    #    print 'Detected: running from terminus upstream.'
    
    SEarr = []
    thickarr = []
    basearr = []
    obsarr = []
    
    SEarr.append(hinit)
    thickarr.append(hinit-(bedfunction(startpoint)/H0))
    basearr.append(bedfunction(startpoint)/H0)
    obsarr.append((obsheightfunction(startpoint))/H0)
    for x in horiz[1::]:
        bed = bedfunction(x)/H0  # value of interpolated bed function
        obsheight = (obsheightfunction(x))/H0
        modelthick = thickarr[-1]
        B = Bfunction(tau_val, bed, modelthick, None, None)
        #Break statements for thinning below yield, water balance, or flotation
        if dx<0:
            if modelthick<BalanceThick(bed,B):
                print 'Thinned below water balance at x=' + str(10*x)+'km'
                break
        if modelthick<FlotationThick(bed):
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
    return (error, CVrms)

def CV_Optimise(linename, lineref, testrange):
    #
    #ty_arr = arange(100e3, 451e3, 5e3)
    #t0_arr = arange(80e3, 431e3, 5e3)
    
    CV_const_arr = []
    CV_var_arr = []
    
    bedf = lineref['bed']
    thickf = lineref['thickness']
    sef = lineref['surface']
    
    arcmax = ArcArray(lineref['line'])[-1]
    
    for tau in testrange:
        #tau_yield = ty_arr[j]
        #tau_0 = t0_arr[j]
        #tau_yield = tau
        #tau_0 = tau
            
        print str(linename) +', no smoothing, Ty = {} Pa'.format(tau)
        tau_yield = tau
        #model_const = plasticmodel_error(tau, bedf, B_const, 0, BalanceThick((bedf(0)/H0), B_const(tau, bedf(0)/H0, thickf(0)/H0, 0, 0))+(bedf(0)/H0), arcmax, 25000, sef)
        model_const = plasticmodel_error(bedf, tau, B_const, 0, sef(0)/H0, arcmax, 25000, sef) #prescribed terminus thickness
        #RMS_const = model_const[0]
        CV_const = model_const[1]
        print str(linename) + ', no smoothing, variable with T0 = {} Pa'.format(tau)
        tau_0 = tau
        #model_var = plasticmodel_error(tau, bedf, B_var, 0, BalanceThick((bedf(0)/H0), B_var(tau, bedf(0)/H0, thickf(0)/H0, 0, 0))+(bedf(0)/H0), arcmax, 25000, sef)
        model_var = plasticmodel_error(bedf, tau, B_var, 0, sef(0)/H0, arcmax, 25000, sef) #prescribed terminus thickness
        #RMS_var = model_var[0]
        CV_var = model_var[1]
        
        CV_const_arr.append(CV_const)
        CV_var_arr.append(CV_var)
    
    constopt_index = np.argmin(CV_const_arr)
    varopt_index = np.argmin(CV_var_arr)
    
    constopt = testrange[constopt_index]
    varopt = testrange[varopt_index]
    
    print 'OPTIMAL VALUE FOR CONSTANT TY: '+str(0.001*constopt)+' kPa'
    print 'OPTIMAL VALUE FOR VARIABLE T0: '+str(0.001*varopt)+' kPa'
    
    return constopt, varopt

def Network_CV_Optimise(networklist, taurange, glaciername='Glacier'):
    """networklist = list of dictionaries, where each dictionary is a flowline
    taurange = range of yield strength values to try 50-500 kpa is a good choice
    
    Modifies the dictionaries in networklist to add entries for optimal tau values
    Returns list of best-fit tau_y, tau_0 values for all branches tested
    """
    bestfitarr = []
    for j, d in enumerate(networklist):    
        optimal_ty, optimal_t0 = CV_Optimise(glaciername+str(j), d, taurange)
        d['Best tau_y'] = optimal_ty
        d['Best tau_0'] = optimal_t0
        bestfitarr.append((optimal_ty, optimal_t0))
    
    return bestfitarr

##-------------------------------
##  TIME EVOLUTION MODELLING
##-------------------------------

def ProcessDicts(linedicts, keys, fields, bestfit_tau):
    """Processing list of flowline dicts to be ready for PlasticEvol"""
    for d in linedicts:
        for j,k in enumerate(keys):
            d[k] = FlowProcess(d['line'], fields[j])
    
    for n,d in enumerate(linedicts):
        tau_0 = bestfit_tau[n][0]
        tau_y = bestfit_tau[n][1]
        arcmax = ArcArray(d['line'])[-1]
        modelprof = PlasticProfile(d['bed'], tau_0, B_var, 0, d['surface'](0)/H0, arcmax, 10000, d['surface'])
        modelint = interpolate.interp1d(modelprof[0], modelprof[1], kind='linear', copy=True)
        d['Modelled'] = modelprof
        d['Ref-profile-func'] = modelint
        d['Best tau_y'] = tau_y
        d['Best tau_0'] = tau_0
        
    return linedicts


def PlasticEvol(linedicts, testyears, upgl_ref=15000/L0, thinrate=10/H0, thinvalues=None):
    """linedicts: a list of flowline dictionaries.  These should be already optimised and include reference profiles from a ref model run 
    testyears: a range of years to test
    upgl_ref: where to apply upglacier thinning.  Default is 15km upstream, or top of glacier if flowline <15km
    thinrate: thinning rate (constant) to apply at reference point
    thinfunc: the option to define thinning as a function fit to obs (e.g. sinusoid) or as extreme climate scenario (e.g. exponential increase in thinning)
    
    returns list of dictionaries with model output
    """
    if thinvalues is None:  
        thinvals = np.full(len(testyears), thinrate)
    else:
        thinvals = thinvalues
    
    modeldicts = [{} for j in range(len(linedicts))]
    for j,d in enumerate(linedicts):
        print 'Currently running line {}'.format(j)
        sarr = d['Modelled'][0] #calling list of master glacier dicts for initialization before getting into modeldicts...
        amax = sarr[-1] #can change if want to model shorter profile
        refpt = min(amax, upgl_ref)
        refht = d['Ref-profile-func'](refpt)
        
        bedf = d['bed']
        sef = d['surface']
        
        tau_j = d['Best tau_0']
        
        dmodel = modeldicts[j]
        dmodel['Termini'] = [L0*min(sarr)]
        dmodel['Termrates'] = []
        
        for j, yr in enumerate(testyears):
            #thinning = yr*thinrate
            thinning = np.sum(thinvals[:j])
            fwdmodel = PlasticProfile(bedf, tau_j, B_var, refpt, refht-thinning, 0, 25000, sef)
            bkmodel = PlasticProfile(bedf, tau_j, B_var, refpt, refht-thinning, amax, 25000, sef)
            modelterm = L0*min(fwdmodel[0]) #in m
            dL = modelterm - dmodel['Termini'][-1]
            #dmodel[yr] = fwdmodel #showing profile downstream of refpt
            dmodel['Termini'].append(modelterm)
            dmodel['Termrates'].append(dL) #dt = 1 by definition
    
    return modeldicts


#def PlasticFluxEvol(linedicts, testyears, upgl_ref=15000/L0, thinrate=10/H0, thinfunc=None):
#    """linedicts: a list of flowline dictionaries.  These should be already optimised and include reference profiles from a ref model run 
#    testyears: a range of years to test
#    upgl_ref: where to apply upglacier thinning.  Default is 15km upstream, or top of glacier if flowline <15km
#    thinrate: thinning rate (constant) to apply at reference point
#    thinfunc: the option to define thinning as a function fit to obs (e.g. sinusoid) or as extreme climate scenario (e.g. exponential increase in thinning)
#    
#    returns list of dictionaries with model output
#    """
#    lines = [ld['line'][:] for ld in linedicts]
#    catchment = spatial.ConvexHull(lines)
#    
#    dt = mean(diff(testyears))
#    bdot = #accumulation rate 
#    