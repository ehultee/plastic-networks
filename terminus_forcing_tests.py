## Testing upstream breakage & other necessary investigations to switch forcing to terminus
## 23 Apr 2018  EHU
##Edited 8 May 2018

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
from plastic_utilities_v2 import *
from GL_model_tools import *
from flowline_class_hierarchy import *

##-------------------
### READING IN BED
### COMMENT OUT IF DATA IS ALREADY READ IN TO YOUR SESSION
##-------------------

print 'Reading in surface topography'
gl_bed_path ='Documents/1. Research/2. Flowline networks/Model/Data/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
b_raw = fh.variables['bed'][:].copy() # bed topo
thick_mask = fh.variables['mask'][:].copy()
ss = np.ma.masked_where(thick_mask !=2, s_raw)#mask values: 0=ocean, 1=ice-free land, 2=grounded ice, 3=floating ice, 4=non-Greenland land
hh = np.ma.masked_where(thick_mask !=2, h_raw) 
bb = np.ma.masked_where(thick_mask !=2, b_raw)
## Down-sampling
X = xx[::2]
Y = yy[::2]
S = ss[::2, ::2]
H = hh[::2, ::2] 
B = bb[::2, ::2]
## Not down-sampling
#X = xx
#Y = yy
#S = ss
fh.close()

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], S.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], B.T[::, ::-1])

##-------------------
### INITIALIZING JAKOBSHAVN
##-------------------

jakcoords_main = Flowline_CSV('Documents/GitHub/plastic-networks/jakobshavn-mainline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_0 = Flowline(coords=jakcoords_main, index=0, name='Jak mainline', has_width=True)
Jakobshavn_main = PlasticNetwork(name='Jakobshavn Isbrae [main/south]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_main[0])
Jakobshavn_main.load_network(filename='Jakobshavn_main-w_35yr_model_output.pickle', load_mainline_output=True)

#jakcoords_sec = Flowline_CSV('Jakobshavn_secondary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
#jak_1 = Flowline(coords=jakcoords_sec, index=1, name='Jak secondary branch', has_width=True)
#Jakobshavn_sec = PlasticNetwork(name='Jakobshavn Isbrae [secondary/central]', init_type='Flowline', branches=(jak_1), main_terminus=jakcoords_sec[0])
#Jakobshavn_sec.load_network(filename='Jakobshavn_sec.pickle', load_mainline_output=False)
#
#jakcoords_tert = Flowline_CSV('Jakobshavn_tertiary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
#jak_2 = Flowline(coords=jakcoords_tert, index=2, name='Jak tertiary branch', has_width=True)
#Jakobshavn_tert = PlasticNetwork(name='Jakobshavn Isbrae [tertiary/north]', init_type='Flowline', branches=(jak_2), main_terminus=jakcoords_tert[0])
#Jakobshavn_tert.load_network(filename='Jakobshavn_tert.pickle', load_mainline_output=False)

#for gl in (Jakobshavn_main, Jakobshavn_sec, Jakobshavn_tert):
for gl in (Jakobshavn_main,):
    gl.process_full_lines(B_interp, S_interp, H_interp)
    for fln in gl.flowlines:
            fln.yield_type = gl.network_yield_type
            fln.optimal_tau = gl.network_tau
    gl.network_ref_profiles()


def check_upstream(network, which_years):
    """Checks upstream for multiple points of cliff failure
    Input:
        network: a PlasticNetwork instance with model output
        which_years: Array of which years of simulation to check
    """
    yieldtest = {}
    mainline = network.flowlines[0]
    #arc = ArcArray(mainline.coords)
    
    for yr in which_years:
        xarr = network.model_output[0][yr][0]
        thickness_array = np.array(network.model_output[0][yr][1]) - np.array(network.model_output[0][yr][2])
        thickness_function = interpolate.interp1d(xarr, thickness_array)
        yieldtest[yr] = []
        for pt in xarr:
            bed = mainline.bed_function(pt)/H0
            model_thickness = thickness_function(pt)
            B = mainline.Bingham_num(bed, model_thickness)
            balancethick = BalanceThick(bed, B)
            yieldtest[yr].append(model_thickness - balancethick)
        breakpoint1 = argwhere(sign(yieldtest[yr])>0)[0]
        breakpoint2 = argwhere(sign(yieldtest[yr])>0)[-1]
        print yr
        print breakpoint1, breakpoint2
    
    return yieldtest
        

def find_dHdL(flowline, profile, dL=None):
    """Function to compute successive profiles of length L-dL, L, L+dL to calculate dHdL over glacier flowline.
    Input: 
        profile: a plastic profile output from Flowline.plasticprofile of length L
        dL: spatial step to use in calculating dHdL.  Default 5 meters
    """
    if dL is None:
        dL = 5/flowline.L0 #nondimensional
    
    xmin = min(profile[0])
    xmax = max(profile[0])
    L_init = xmax-xmin
    
    #Terminus quantities
    SE_terminus = profile[1][0] #NEED TO CONFIRM that terminus is at [0] and not [-1]
    Bed_terminus = profile[2][0]
    H_terminus = SE_terminus - Bed_terminus 
    print H_terminus
    Bghm_terminus = flowline.Bingham_num(Bed_terminus, H_terminus)
    
    #Profile advanced by dL - note coord system means xmin-dL is more advanced, as x=0 is at initial terminus position
    bed_mindL = (flowline.bed_function(xmin-dL))/flowline.H0
    s_mindL = BalanceThick(bed_mindL, Bghm_terminus) + bed_mindL
    profile_mindL = flowline.plastic_profile(startpoint=xmin-dL, hinit = s_mindL, endpoint = xmax, surf = flowline.surface_function)
    H_mindL = np.array(profile_mindL[1]) - np.array(profile_mindL[2]) #array of ice thickness from profile
    Hx_mindL = interpolate.interp1d(profile_mindL[0], H_mindL, bounds_error=False, fill_value=0)
    print Hx_mindL(xmin)
    
    #Profile retreated by dL
    bed_plusdL = (flowline.bed_function(xmin+dL))/flowline.H0
    s_plusdL = BalanceThick(bed_plusdL, Bghm_terminus) + bed_plusdL
    profile_plusdL = flowline.plastic_profile(startpoint = xmin+dL, hinit = s_plusdL, endpoint = xmax, surf=flowline.surface_function)
    H_plusdL = np.array(profile_plusdL[1]) - np.array(profile_plusdL[2]) #array of ice thickness
    Hx_plusdL = interpolate.interp1d(profile_plusdL[0], H_plusdL, bounds_error=False, fill_value=0)
    print Hx_plusdL(xmin)
    
    dHdLx = lambda x: (Hx_plusdL(x) - Hx_mindL(x))/(2*dL)
    
    return dHdLx


#def terminus_dUdx(flowline, profile, A_factor=3.5E-25):
#    """Function to evaluate terminus dU/dx, needed for dLdt.
#    Input:
#        profile: a plastic profile output from Flowline.plasticprofile of the current time step
#        A: flow rate factor, assumed 3.5x10^(-25) Pa^-3 s^-1 for T=-10C based on Cuffey & Paterson
#    """
#    #Glen's flow law
#    dUdx_invseconds = A_factor * (flowline.optimal_tau)**3
#    s_per_annum = 31557600 #unit conversion to make calculation agree with accumulation
#    
#    dUdx_invannum = dUdx * s_per_annum
#    
#    return dUdx_invannum
#    

def balance_adot(network, V_field):
    """Function to compute spatially-averaged accumulation rate that balances observed terminus velocity
    Input:
        V_field: 2d-interpolated function to report magnitude of velocity (i.e. ice flow speed in m/a) given (x,y) coordinates
    """
    terminus = network.flowlines[0].coords[0][0:2] #All lines of a network should share a terminus at initial state
    terminus_speed = V_field(terminus) #dimensional in m/a
    terminus_width = network.flowlines[0].width_function(0) #dimensional in m
    terminus_thickness = network.H0*network.flowlines[0].ref_profile(0) - network.flowlines[0].bed_function(0) #dimensional in m
    termflux = terminus_speed * terminus_width * terminus_thickness
    
    sa = []
    for fl in network.flowlines:
        L = fl.length
        surface_area_nondim = quad(fl.width_function, 0, L)[0]
        surface_area = network.L0 * surface_area_nondim
        sa.append(surface_area)
    catchment_area = sum(sa)
    
    balance_a = termflux/catchment_area
    

    

def dLdt(flowline, profile, a_dot, rate_factor=3.5E-25, dL=None):
    """Function to compute terminus rate of advance/retreat given a mass balance forcing, a_dot.
    Input:
        profile: a plastic profile output from Flowline.plasticprofile of the current time step
        a_dot: net rate of ice accumulation/loss (in m/a)--spatially averaged over whole catchment for now
        rate_factor: flow rate factor A, assumed 3.5x10^(-25) Pa^-3 s^-1 for T=-10C based on Cuffey & Paterson
    """      
    xmin = min(profile[0])
    xmax = max(profile[0])
    L = xmax-xmin #length of the current profile, nondimensional
    print L
    
    if dL is None:
        dL=5/flowline.L0
    
    dHdL = find_dHdL(flowline, profile, dL)
    
    #Terminus quantities
    SE_terminus = profile[1][0] #terminus at [0], not [-1]--may return errors if running from head downstream, but this is for terminus forcing anyway
    print 'SE_terminus={}'.format(SE_terminus)
    Bed_terminus = profile[2][0]
    print 'Bed_terminus={}'.format(Bed_terminus)
    H_terminus = SE_terminus - Bed_terminus 
    Bghm_terminus = flowline.Bingham_num(Bed_terminus, H_terminus)
    Hy_terminus = BalanceThick(Bed_terminus, Bghm_terminus)
    print 'Hy_terminus={}'.format(Hy_terminus)

    #Quantities at adjacent grid point
    SE_adj = profile[1][1]
    Bed_adj = profile[2][1]
    H_adj = SE_adj - Bed_adj
    Bghm_adj = flowline.Bingham_num(Bed_adj, H_adj)
    Hy_adj = BalanceThick(Bed_adj, Bghm_adj)
    
    #Diffs
    dx_term = abs(profile[0][1] - profile[0][0]) #should be ~2m in physical units
    print 'dx_term={}'.format(dx_term)
    dHdx = (H_adj-H_terminus)/dx_term
    dHydx = (Hy_adj-Hy_terminus)/dx_term
    s_per_annum = 31557600 #unit conversion to make calculation agree with accumulation
    tau = flowline.Bingham_num(Bed_terminus, H_terminus) * (flowline.rho_ice * flowline.g * flowline.H0**2 / flowline.L0) #using Bingham_num handles whether tau_y constant or variable for selected flowline
    dUdx_terminus = s_per_annum * rate_factor * tau**3 

    Area_int = quad(dHdL, xmin, xmax)[0]
    print 'Area_int={}'.format(Area_int)
    #print 'dH/dL at terminus = {}'.format(dHdL(xmin))
    
    denom = dHydx - dHdx - (H_terminus**(-1))*dHdx*Area_int
    numerator = a_dot - dUdx_terminus*H_terminus - (a_dot*L*dHdx/H_terminus)
    
    result = numerator/denom
    
    return result

def terminus_time_evolve(network, testyears=arange(100), ref_branch_index=0, a_dot, a_dot_variable=None, upstream_limits=None, use_mainline_tau=True):
    """Time evolution on a network of Flowlines, forced from terminus.  Lines should be already optimised and include reference profiles from network_ref_profiles
    Arguments:
        testyears: a range of years to test, indexed by years from nominal date of ref profile (i.e. not calendar years)
        ref_branch_index: which branch to use for forcing.  Default is main branch ("0") but can be changed
        a_dot: spatially averaged accumulation rate (forcing)
    Optional args:   
        a_dot_variable: array of the same length as testyears with a different a_dot forcing to use in each year
        Offers the option to define thinning as persistence of obs or other nonlinear function.
        upstream_limits: array determining where to cut off modelling on each flowline, ordered by index.  Default is full length of lines.
        use_mainline_tau=False will force use of each line's own yield strength & type
    
        returns model output as dictionary for each flowline 
    """

    #Fixing default values
    if upstream_limits is None:
        upstream_limits=[fl.length for fl in network.flowlines]
    
    if a_dot_variable is None:  
        a_dot_vals = np.full(len(testyears), a_dot)
    else:
        a_dot_vals = a_dot_variable
    
    dt = mean(diff(testyears)) #size of time step
    
    model_output_dicts = [{'Termini': [0],
    'Terminus_heights': [fl.surface_function(0)],
    'Termrates': [],
    'Terminus_flux': []
    } for fl in network.flowlines]
    
    #Mainline reference
    ref_line = network.flowlines[ref_branch_index]
    ref_amax = upstream_limits[ref_branch_index]
    ref_surface = ref_line.surface_function
    #refpt = min(ref_amax, upgl_ref) #apply forcing at top of branch if shorter than reference distance.  In general would expect forcing farther downstream to give weaker response
    #refht = ref_line.ref_profile(refpt)
    if use_mainline_tau:
        ref_line.optimal_tau = network.network_tau
        ref_line.yield_type = network.network_yield_type
    else:
        pass
    refdict = model_output_dicts[ref_branch_index]
    refdict[0] = ref_line.ref_profile #initial condition for time evolution - needed to calculate calving flux at first timestep

    
    #Assume same terminus
    for k, yr in enumerate(testyears):
        a_dot_k = a_dot_vals[k]
        
        dLdt_annum = dLdt(flowline=ref_line, profile=refdict[k-1], a_dot=a_dot_k)
        #Ref branch
        new_termpos = network.L0*(refdict['Termini'][-1]+(dLdt_annum*dt)) #Multiply by dt in case dt!=1 annum
        new_term_bed = ref_line.bed_function(new_termpos/network.L0)
        previous_bed = ref_line.bed_function(refdict['Termini'][-1]/network.L0)
        previous_thickness = (refdict['Terminus_heights'][-1] - previous_bed)/network.H0 #nondimensional thickness for use in Bingham number
        new_termheight = BalanceThick(new_term_bed/network.H0, ref_line.Bingham_num(previous_bed/network.H0, previous_thickness)) + (new_term_bed/network.H0)
        new_profile = ref_line.plastic_profile(startpoint=new_termpos/network.L0, hinit=new_termheight, endpoint=ref_amax, surf=ref_surface)
        if yr>dt:
            termflux = ref_line.icediff(profile1=refdict[yr-dt], profile2=new_profile)
        else:
            termflux = np.nan
            
        refdict[yr] = new_profile
        refdict['Terminus_flux'].append(termflux)
        refdict['Termini'].append(new_termpos)
        refdict['Terminus_heights'].append(new_termheight)
        refdict['Termrates'].append(dLdt_annum*dt)
        
        #Other branches, incl. branch splitting
        for j, fl in enumerate(network.flowlines):
            out_dict = model_output_dicts[j]
            fl_amax = upstream_limits[j]
            separation_distance = ArcArray(fl.coords)[fl.intersections[1]] #where line separates from mainline
            if j==ref_branch_index:
                continue
            else:
                if use_mainline_tau:
                    fl.optimal_tau = network.network_tau
                    fl.yield_type = network.network_yield_type #this is probably too powerful, but unclear how else to exploit Bingham_number functionality
                else:
                    pass
                
                if out_dict['Termini'][-1]/network.L0 <= separation_distance : ## Below is only applicable while branches share single terminus 
                    dLdt_branch = dLdt_annum
                    branch_terminus = new_termpos
                    branch_termheight = new_termheight
                else: ##if branches have split, find new terminus quantities
                    dLdt_branch = dLdt(flowline=fl, profile=out_dict[k-1], a_dot=a_dot_k)
                    branch_terminus = network.L0*(out_dict['Termini'][-1] +(dLdt_branch*dt))
                    branch_term_bed = fl.bed_function(branch_terminus/network.L0)
                    previous_branch_bed = fl.bed_function(out_dict['Termini'][-1]/network.L0)
                    previous_branch_thickness = (out_dict['Terminus_heights'][-1] - previous_branch_bed)/network.H0
                    branch_termheight = BalanceThick(branch_term_bed/network.H0, fl.Bingham_num(previous_branch_bed/network.H0, previous_branch_thickness)) + (branch_term_bed/network.H0)
                    
                branchmodel = fl.plastic_profile(startpoint=branch_terminus/network.L0, hinit=branch_termheight, endpoint=fl_amax, surf=fl.surface_function)
                out_dict[yr] = branchmodel
                out_dict['Termini'].append(branch_terminus)
                out_dict['Terminus_heights'].append(branch_termheight)
                out_dict['Termrates'].append(dLdt_branch*dt)
                if yr > dt:
                    out_dict['Terminus_flux'].append(termflux)
                else:
                    out_dict['Terminus_flux'].append(np.nan)
    
    network.model_output = model_output_dicts
