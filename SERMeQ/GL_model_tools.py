## Utilities file for Greenland modelling functions
## 30 Nov 2017  EHU
## 10 Jul 2019  Adding visualization tools

from numpy import *
#from netCDF4 import Dataset
import numpy as np
#from scipy import interpolate
from scipy import spatial
#from scipy.ndimage import gaussian_filter
from shapely.geometry import *
from mpl_toolkits.basemap import Basemap
# import mpl_toolkits.basemap.pyproj as pyproj
import pyproj
import matplotlib.colors as colors
import shapefile
from osgeo import gdal
import cPickle as pickle
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
def Flowline_CSV(filename, nlines=None, has_width=False, flip_order=True):
    """Function to read in flowlines in CSV format, similar to those sent by C. Kienholz for Alaska.
    Input: filename; nlines=number of flowlines.  Edited to make this argument unnecessary, but leaving it in for backward compatibility
        has_width: default False for older files that have only (x,y) rather than (x,y,width) saved
        flip_order: default False for lines that already run from terminus to peak
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
    #if nlines is not None:
    #    data['Lineslist'] = [[] for k in range(nlines)] 
    data['Lineslist'] = [] #initialize as empty list
    
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
            data['Lineslist'].append(temp) #initialize new empty array that can be modified in-place later
            data['Length_ID'][parts[0]] = 1
            j+=1 
        else:
            data['Length_ID'][parts[0]] += 1
        #if xbounds[0]<x_coord<xbounds[1]: #taking out values outside of map area
        #    if ybounds[0]<y_coord<ybounds[1]:   
        
        if has_width:
            width = float(parts[3])
            temp.append((x_coord, y_coord, width))
        else:
            temp.append((x_coord, y_coord))
            
        data['Lineslist'][j-1] = np.array(temp) #need to modify an existing array rather than append to keep correct indexing

        #data['Lineslist'][j] = np.array(temp) 
    
    if nlines is None:
        nlines = len(data['Length_ID'].keys())
    
    if flip_order:          
        centrelines_list = [np.array(data['Lineslist'][j])[::-1] for j in range(nlines)] #making arrays, reversed to start at terminus rather than peak
    else:
        centrelines_list = [np.array(data['Lineslist'][j]) for j in range(nlines)] # arrays already start at terminus

    
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


def plasticmodel_error(bedfunction, tau_val, Bfunction, startpoint, hinit, endpoint, Npoints, obsheightfunction, allow_upstream_breakage=True):
    """Arguments used:
        bedfunction should be function of arclength returning bed elevation of the glacier.
        Bfunction is nondim yield strength.  Should be function with arguments elevation, ice thickness, position, and time (can just not use last two if no variation)
        Startpoint is where (in arclength space) integration should start.
        hinit is nondim initial height.  Could be given by water balance, obs, or some thinning from reference height.
        Endpoint is where integration should stop.
        Npoints is how many model points to use (suggested 25000+)
        #Resolution (in m) is how closely we want to space the model sample points (CURRENTLY USING NPOINTS INSTEAD OF RESOLUTION)
        Obsheightfunction is the observations for comparison.  May need to process from a thickness measurement. (working on functionality to make this argument optional)
        allow_upstream_breakage: determines whether profiles should be allowed to break (due to yield) when stepping upstream--default is True, but False may allow optimization with more sparse data
        
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
            if modelthick<FlotationThick(bed) and allow_upstream_breakage: #new control on whether breaking happens
                print 'Thinned below flotation at x=' + str(10*x)+'km'
                break
            if modelthick<4*B*H0/L0 and allow_upstream_breakage:
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
##  TIME EVOLUTION MODELLING - superseded by flowline_class_hierarchy code
##-------------------------------
#
#def ProcessDicts(linedicts, keys, fields, bestfit_tau):
#    """Processing list of flowline dicts to be ready for PlasticEvol"""
#    for d in linedicts:
#        for j,k in enumerate(keys):
#            d[k] = FlowProcess(d['line'], fields[j])
#    
#    for n,d in enumerate(linedicts):
#        tau_0 = bestfit_tau[n][0]
#        tau_y = bestfit_tau[n][1]
#        arcmax = ArcArray(d['line'])[-1]
#        modelprof = PlasticProfile(d['bed'], tau_0, B_var, 0, d['surface'](0)/H0, arcmax, 10000, d['surface'])
#        modelint = interpolate.interp1d(modelprof[0], modelprof[1], kind='linear', copy=True)
#        d['Modelled'] = modelprof
#        d['Ref-profile-func'] = modelint
#        d['Best tau_y'] = tau_y
#        d['Best tau_0'] = tau_0
#        
#    return linedicts
#
#
#def PlasticEvol(linedicts, testyears, upgl_ref=15000/L0, thinrate=10/H0, thinvalues=None):
#    """linedicts: a list of flowline dictionaries.  These should be already optimised and include reference profiles from a ref model run 
#    testyears: a range of years to test
#    upgl_ref: where to apply upglacier thinning.  Default is 15km upstream, or top of glacier if flowline <15km
#    thinrate: thinning rate (constant) to apply at reference point
#    thinfunc: the option to define thinning as a function fit to obs (e.g. sinusoid) or as extreme climate scenario (e.g. exponential increase in thinning)
#    
#    returns list of dictionaries with model output
#    """
#    if thinvalues is None:  
#        thinvals = np.full(len(testyears), thinrate)
#    else:
#        thinvals = thinvalues
#    
#    modeldicts = [{} for j in range(len(linedicts))]
#    for j,d in enumerate(linedicts):
#        print 'Currently running line {}'.format(j)
#        sarr = d['Modelled'][0] #calling list of master glacier dicts for initialization before getting into modeldicts...
#        amax = sarr[-1] #can change if want to model shorter profile
#        refpt = min(amax, upgl_ref)
#        refht = d['Ref-profile-func'](refpt)
#        
#        bedf = d['bed']
#        sef = d['surface']
#        
#        tau_j = d['Best tau_0']
#        
#        dmodel = modeldicts[j]
#        dmodel['Termini'] = [L0*min(sarr)]
#        dmodel['Termrates'] = []
#        
#        for j, yr in enumerate(testyears):
#            #thinning = yr*thinrate
#            thinning = np.sum(thinvals[:j])
#            fwdmodel = PlasticProfile(bedf, tau_j, B_var, refpt, refht-thinning, 0, 25000, sef)
#            bkmodel = PlasticProfile(bedf, tau_j, B_var, refpt, refht-thinning, amax, 25000, sef)
#            modelterm = L0*min(fwdmodel[0]) #in m
#            dL = modelterm - dmodel['Termini'][-1]
#            #dmodel[yr] = fwdmodel #showing profile downstream of refpt
#            dmodel['Termini'].append(modelterm)
#            dmodel['Termrates'].append(dL) #dt = 1 by definition
#    
#    return modeldicts
#
#

##--------------------------
## VISUALIZATION - MAPPING
##--------------------------

def Greenland_map(service='ESRI_Imagery_World_2D', epsg=3413, xpixels=5000):
    """Function using Basemap to plot map for all of Greenland.
    Input:
        service: map appearance selected from ['World_Physical_Map', 'World_Shaded_Relief', 'World_Topo_Map', 'NatGeo_World_Map', 'ESRI_Imagery_World_2D', 'World_Street_Map', 'World_Imagery', 'Ocean_Basemap']
        epsg: identifier of specific map projection to use in plotting.  Default is 3413 (Polar Stereographic North).
    """
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=300, llcrnrlat=57, urcrnrlon=20, urcrnrlat=80, resolution='h')
    
    plt.figure()
    m.arcgisimage(service=service, xpixels=xpixels)
    plt.show()
    return m

##Convert coords into lat/lon so that Basemap can convert them back (don't know why this is necessary, but it works)
def flowline_latlon(coords, fromproj=pyproj.Proj("+init=epsg:3413"), toproj=pyproj.Proj("+init=EPSG:4326")):
    """Convert coords into lat/lon so that Basemap can convert them back for plotting (don't know why this is necessary, but it works)
    Defaults:
        fromproj = NSIDC Polar Stereographic North, EPSG 3413
        toproj = WGS84 lat-lon, EPSG 4326
    """
    xs = coords[:,0]
    ys = coords[:,1]
    x_lon, y_lat = pyproj.transform(fromproj, toproj, xs, ys)
    latlon_coords = np.asarray(zip(x_lon, y_lat))
    return latlon_coords

# set the colormap and centre the colorbar - from Joe Kington (StOv)
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

##--------------------------
## GREENLAND-SPECIFIC FILE I/O
##--------------------------

def read_termini(filename, year):
    """Make and return a dictionary of terminus positions, indexed by MEaSUREs ID.  These can then be plotted on a Greenland_map instance
    Input:
        filename = name of MEaSUREs terminus position shapefile to read
        year: year of the terminus position observations"""
    print 'Reading in MEaSUREs terminus positions for year ' + str(year)
    sf = shapefile.Reader(filename)
    fields = sf.fields[1:] #excluding the mute "DeletionFlag"
    field_names = [field[0] for field in fields]
    term_recs = sf.shapeRecords()
    termpts_dict = {}
    for r in term_recs:
        atr = dict(zip(field_names, r.record)) #dictionary of shapefile fields, so we can access GlacierID by name rather than index.  Index changes in later years.
        key = atr['GlacierID'] #MEaSUREs ID number for the glacier, found by name rather than index
        termpts_dict[key] = np.asarray(r.shape.points) #save points spanning terminus to dictionary
    return termpts_dict

# Finding intersection of terminus points with mainline--modified from hubbard-mainline-advance-v2.py
def projected_term_obs(termset, linestr):
    '''Given a termset from file input and LineString representation of a flowline, termline constructs a Shapely LineString of the terminus and returns the intersection of the two'''
    termarr = np.array(termset)
    termline = LineString(termarr)
    centrpt = termline.centroid
    arcdist = linestr.project(centrpt)
    if arcdist>0:
        return arcdist/1000
    else:
        near = linestr.distance(termline)  #in case terminus listed in MEaSUREs is farther advanced than max seaward extent of saved flowline
        return -near/1000
    
def advterm(termset, linestr):
    '''Given termset and LineString representation of a flowline, advterm finds which terminus position projects most advanced along central flowline and returns its arclength position'''
    x_term = termset[:, 0]  #Note need to change from [:, 1] to [:, 0] for x-coord, due to different data format for Hubbard
    y_term = termset[:, 1]
    projections = []
    for i in xrange(len(x_term)):
        proji = linestr.project(Point(x_term[i], y_term[i]))
        projections.append(proji)
    termmax = min(projections) #need minimum rather than max here because we are interested in the most advanced, i.e. lowest arc-length value projection of terminus
    return termmax/1000

def retterm(termset, linestr):
    '''Given termset (from file input above), retterm finds which terminus position projects most retreated (rel. 2007 terminus) along central flowline and returns its arclength position'''
    x_term = termset[:, 0]
    y_term = termset[:, 1]
    projections = []
    for i in xrange(len(x_term)):
        proji = linestr.project(Point(x_term[i], y_term[i]))
        projections.append(proji)
    termmin = max(projections) #max is the most retreated, i.e. highest arc-length value projection of terminus
    return termmin/1000



##Function to read MEaSUREs velocity GeoTIFFs
def read_velocities(filename, return_grid=True, return_proj=False):
    """Extract x, y, v from a MEaSUREs GeoTIFF.
    Input: 
        filename = GeoTIFF to be read
    Optional args:
        return_grid = whether to return x-y grid (default True) or only the velocity field (False)
        return_proj = whether to return the gdal projection parameters (default False)"""
    ds = gdal.Open(filename)
    #Get dimensions
    nc = ds.RasterXSize
    nr = ds.RasterYSize
    
    geotransform = ds.GetGeoTransform()
    xOrigin = geotransform[0]
    xPix = geotransform[1] #pixel width in x-direction
    yOrigin = geotransform[3]
    yPix = geotransform[5] #pixel height in y-direction
    
    lons = xOrigin + np.arange(0, nc)*xPix
    lats = yOrigin + np.arange(0, nr)*yPix
    
    x, y = np.meshgrid(lons, lats)
    
    vband = ds.GetRasterBand(1)
    varr = vband.ReadAsArray()
    
    if return_grid and return_proj:
        return x, y, varr, ds.GetProjection()
    elif return_grid:
        return x, y, varr
    else: 
        return varr

### Load-in functionality to read only terminus position and flux, lifted from Greenland-automated_summary_plots.py
def lightload(filename, glacier_name, output_dictionary):
    """Function to read only terminus position and flux from stored plastic model output.
    Input: 
        filename = the name of a pickle file with stored model output
        glacier_name = name or other identifier of the glacier to be read
        output_dictionary = name of an existing dictionary where we should put this output.
    Returns:
        output_dictionary modified to add the requested model output
    """
    output_dictionary[glacier_name] = {}
    
    with open(filename, 'rb') as handle:
        loadin = pickle.load(handle)
    
    N_Flowlines = loadin['N_Flowlines']
    mainline_termini = loadin['mainline_model_output']['Termini']
    mainline_flux = loadin['mainline_model_output']['Terminus_flux']
    output_dictionary[glacier_name][0] ={'Termini': mainline_termini, 'Terminus_flux': mainline_flux}
    
    if N_Flowlines >1:
        for n in range(N_Flowlines)[1::]:
            key_n = 'model_output_'+str(n)
            termini_n = loadin[key_n]['Termini']
            termflux_n = loadin[key_n]['Terminus_flux']
            output_dictionary[glacier_name][n] = {'Termini': termini_n, 'Terminus_flux': termflux_n}
    else:
        pass
        
    return output_dictionary

def scenario_cumulative_SLE(scenario_dictionary):
    sd = scenario_dictionary
    pernetwork_cumul_fx = []
    pernetwork_cumul_sle = []
    for j, gid in enumerate(glaciers_simulated):
        branch_fx = [np.nan_to_num(sd['GID{}'.format(gid)][k]['Terminus_flux']) for k in range(len(sd['GID{}'.format(gid)]))]
        total_fx = sum(branch_fx, axis=0)
        total_sle = (1E-12)*np.array(total_fx)/(361.8) #Gt ice/mm SLE conversion
        cumul_fx = np.cumsum(total_fx)
        cumul_sle = np.cumsum(total_sle)
        pernetwork_cumul_fx.append(cumul_fx)
        pernetwork_cumul_sle.append(cumul_sle)
    scenario_sle = np.cumsum(pernetwork_cumul_sle, axis=0)
    return scenario_sle
    
def compare_scenario_SLE(full_output_dictionary):
    """Calculate scenario_cumulative_SLE for all scenarios simulated, and compare them"""
    perscenario_SLE = []
    
    for s in full_output_dictionary.keys():
        print 'Scenario {}'.format(s)
        perscenario_SLE.append(scenario_cumulative_SLE(full_output_dictionary(s))[-1])
    
    return perscenario_SLE
    

    