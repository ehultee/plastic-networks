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

##Nondimensional time scaling
T_0 = 1 #annum, taken as the characteristic time scale
s_per_annum = 31557600 #seconds per annum, to convert time unit where necessary


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


##-------------------
### INITIALIZING JAKOBSHAVN
##-------------------

jakcoords_main = Flowline_CSV('Documents/GitHub/plastic-networks/jakobshavn-mainline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_0 = Flowline(coords=jakcoords_main, index=0, name='Jak mainline', has_width=True)
Jakobshavn_main = PlasticNetwork(name='Jakobshavn Isbrae [main/south]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_main[0])
Jakobshavn_main.load_network(filename='Jakobshavn_main-w_35yr_model_output.pickle', load_mainline_output=False)

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
        

def find_dHdL(flowline, profile, dL=None, debug_mode=False):
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
    
    if xmin-dL > 0:
        x_fwd = xmin-dL #note coord system has x=0 at the terminus, so xmin-dL is a more advanced position
        x_bk = xmin+dL
    else:
        x_fwd = xmin
        x_bk = xmin + 2*dL
    
    
    #Terminus quantities
    SE_terminus = profile[1][0] #CONFIRM that terminus is at [0] and not [-1]
    Bed_terminus = profile[2][0]
    H_terminus = SE_terminus - Bed_terminus 
    Bghm_terminus = flowline.Bingham_num(Bed_terminus, H_terminus)
    
    #Profile advanced by dL - note coord system means xmin-dL is more advanced, as x=0 is at initial terminus position
    bed_mindL = (flowline.bed_function(x_fwd))/flowline.H0
    s_mindL = BalanceThick(bed_mindL, Bghm_terminus) + bed_mindL
    profile_mindL = flowline.plastic_profile(startpoint=x_fwd, hinit = s_mindL, endpoint = xmax, surf = flowline.surface_function)
    H_mindL = np.array(profile_mindL[1]) - np.array(profile_mindL[2]) #array of ice thickness from profile
    Hx_mindL = interpolate.interp1d(profile_mindL[0], H_mindL, bounds_error=False, fill_value=0)
    
    #Profile retreated by dL
    bed_plusdL = (flowline.bed_function(x_bk))/flowline.H0
    s_plusdL = BalanceThick(bed_plusdL, Bghm_terminus) + bed_plusdL
    profile_plusdL = flowline.plastic_profile(startpoint = x_bk, hinit = s_plusdL, endpoint = xmax, surf=flowline.surface_function)
    H_plusdL = np.array(profile_plusdL[1]) - np.array(profile_plusdL[2]) #array of ice thickness
    Hx_plusdL = interpolate.interp1d(profile_plusdL[0], H_plusdL, bounds_error=False, fill_value=0)
    
    dHdLx = lambda x: (Hx_plusdL(x) - Hx_mindL(x))/(2*dL)
    
    if debug_mode:
        print 'Debugging dHdL.  Inspect:'
        print 'H_terminus={}'.format(H_terminus)
        print 'Hx_mindL={}'.format(Hx_mindL(xmin))
        print 'Hx_plusdL={}'.format(Hx_plusdL(xmin))
    else:
        pass
    
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

def balance_adot(network, V_field, use_width=False, L_limit=6.0):
    """Function to compute spatially-averaged accumulation rate that balances observed terminus velocity
    Input:
        V_field: 2d-interpolated function to report magnitude of velocity (i.e. ice flow speed in m/a) given (x,y) coordinates
        use_width: whether to consider upstream width when calculating balance accumulation.  Default is no.
        L_limit: nondimensional upstream distance that indicates how far we trust our model.  Default is 6.0 (60 km in dimensional units).
    """
    terminus = network.flowlines[0].coords[0][0:2] #All lines of a network should share a terminus at initial state
    terminus_speed = V_field(terminus[0], terminus[1]) #dimensional in m/a
    terminus_width = network.flowlines[0].width_function(0) #dimensional in m
    terminus_thickness = network.H0*network.flowlines[0].ref_profile(0) - network.flowlines[0].bed_function(0) #dimensional in m
    balance_thickness = network.H0*BalanceThick(network.flowlines[0].bed_function(0)/network.H0, network.flowlines[0].Bingham_num(network.flowlines[0].bed_function(0)/network.H0, terminus_thickness))
    termflux = terminus_speed * terminus_width * balance_thickness
    

    sa = []
    total_L = []
    for fl in network.flowlines:
        L = min(fl.length, L_limit)
        surface_area_nondim = quad(fl.width_function, 0, L)[0]
        surface_area = network.L0 * surface_area_nondim
        sa.append(surface_area)
        total_L.append(network.L0 * L)
    catchment_area = sum(sa)
    total_ice_length = sum(total_L)
    
    if use_width:
        balance_a = termflux/catchment_area
    else:
        balance_a = terminus_speed*balance_thickness/total_ice_length
    
    return balance_a
    

    

def dLdt(flowline, profile, a_dot, rate_factor=3.5E-25, dL=None, debug_mode=False):
    """Function to compute terminus rate of advance/retreat given a mass balance forcing, a_dot.
    Input:
        profile: a plastic profile output from Flowline.plasticprofile of the current time step
        a_dot: net rate of ice accumulation/loss.  Should be expressed in m/a /H0. Spatially averaged over whole catchment for now
        rate_factor: flow rate factor A, assumed 3.5x10^(-25) Pa^-3 s^-1 for T=-10C based on Cuffey & Paterson
    Returns dLdt in nondimensional units.  Multiply by L0 to get units of m/a (while T0=1a).
    """      
    xmin = min(profile[0])
    xmax = max(profile[0])
    L = xmax-xmin #length of the current profile, nondimensional
    
    if dL is None:
        dL=5/flowline.L0 #step over which to test dHdL profiles
    
    dHdL = find_dHdL(flowline, profile, dL)
    
    #Nondimensionalising rate factor
    inverse_units_of_A = T_0 * (flowline.rho_ice **3)*(flowline.g **3) * (flowline.H0 **6) / (flowline.L0 **3)
    #units_of_A = (flowline.L0 **3)/ (T_0*(flowline.rho_ice **3)*(flowline.g **3) *(flowline.H0 **6))
    #nondim_A = rate_factor * inverse_units_of_A
    
    #Terminus quantities
    SE_terminus = profile[1][0] #terminus at [0], not [-1]--may return errors if running from head downstream, but this is for terminus forcing anyway
    Bed_terminus = profile[2][0]
    H_terminus = SE_terminus - Bed_terminus 
    Bghm_terminus = flowline.Bingham_num(Bed_terminus, H_terminus)
    Hy_terminus = BalanceThick(Bed_terminus, Bghm_terminus)

    #Quantities at adjacent grid point
    SE_adj = profile[1][1]
    Bed_adj = profile[2][1]
    H_adj = SE_adj - Bed_adj
    Bghm_adj = flowline.Bingham_num(Bed_adj, H_adj)
    Hy_adj = BalanceThick(Bed_adj, Bghm_adj)
    
    #Diffs
    dx_term = abs(profile[0][1] - profile[0][0]) #should be ~2m in physical units
    dHdx = (H_adj-H_terminus)/dx_term
    dHydx = (Hy_adj-Hy_terminus)/dx_term
    tau = flowline.Bingham_num(Bed_terminus, H_terminus) * (flowline.rho_ice * flowline.g * flowline.H0**2 / flowline.L0) #using Bingham_num handles whether tau_y constant or variable for selected flowline
    dUdx_terminus = -1 * rate_factor * tau**3 #-1 due to sign convention with x increasing upstream from terminus
    nondim_dUdx_terminus = dUdx_terminus * inverse_units_of_A / (flowline.rho_ice * flowline.g * flowline.H0**2 / flowline.L0) #divide out units to get nondimensional quantity

    Area_int = quad(dHdL, xmin, xmax)[0]
    #print 'dH/dL at terminus = {}'.format(dHdL(xmin))
    
    
    denom = dHydx - dHdx* (1 + (Area_int/H_terminus))
    numerator = a_dot - nondim_dUdx_terminus*H_terminus + (a_dot*L*dHdx/H_terminus)
    
    result = numerator/denom
    
    if debug_mode:
        print 'For inspection on debugging:'
        print 'L={}'.format(L)
        print 'SE_terminus={}'.format(SE_terminus)
        print 'Bed_terminus={}'.format(Bed_terminus)
        print 'Hy_terminus={}'.format(Hy_terminus)
        print 'dx_term={}'.format(dx_term)
        print 'Area_int={}'.format(Area_int)
        print 'Checking dLdt: a_dot = {}. \n H dUdx = {}. \n Ub dHdx = {}.'.format(a_dot, nondim_dUdx_terminus*H_terminus, a_dot*L*dHdx/H_terminus) 
        print 'Denom = {}'.format(denom)
    else:
        pass

    return result

def terminus_time_evolve(network, testyears=arange(100), ref_branch_index=0, a_dot=None, a_dot_variable=None, upstream_limits=None, use_mainline_tau=True, debug_mode=False):
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
    
    if a_dot is None:
        a_dot = 0.2/network.H0
    
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
    ref_surface = ref_line.ref_profile
    #refpt = min(ref_amax, upgl_ref) #apply forcing at top of branch if shorter than reference distance.  In general would expect forcing farther downstream to give weaker response
    #refht = ref_line.ref_profile(refpt)
    if use_mainline_tau:
        ref_line.optimal_tau = network.network_tau
        ref_line.yield_type = network.network_yield_type
    else:
        pass
    refdict = model_output_dicts[ref_branch_index]
    refdict[0] = ref_line.plastic_profile(startpoint=0, hinit=ref_surface(0), endpoint=ref_amax, surf=ref_surface) #initial condition for time evolution - needed to calculate calving flux at first timestep

    
    #Assume same terminus
    for k, yr in enumerate(testyears):
        a_dot_k = a_dot_vals[k]
        
        if k<1:
            dLdt_annum = dLdt(flowline=ref_line, profile=refdict[0], a_dot=a_dot_k) * network.L0
        else:
            dLdt_annum = dLdt(flowline=ref_line, profile=refdict[k-1], a_dot=a_dot_k) * network.L0
        #Ref branch

        new_termpos_raw = refdict['Termini'][-1]-(dLdt_annum*dt) #Multiply by dt in case dt!=1 annum
        new_termpos = max(0, new_termpos_raw)
        if debug_mode:
            print 'dLdt_annum = {}'.format(dLdt_annum)
            print 'New terminus position = {}'.format(new_termpos)
        else:
            pass
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
        refdict['Terminus_heights'].append(new_termheight*network.H0)
        refdict['Termrates'].append(dLdt_annum*dt)
        
        #Other branches, incl. branch splitting
        for j, fl in enumerate(network.flowlines):
            out_dict = model_output_dicts[j]
            fl_amax = upstream_limits[j]
            if j==ref_branch_index:
                continue
            else:
                separation_distance = ArcArray(fl.coords)[fl.intersections[1]] #where line separates from mainline
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
                    dLdt_branch = dLdt(flowline=fl, profile=out_dict[k-1], a_dot=a_dot_k) * network.L0
                    branch_terminus = out_dict['Termini'][-1] -(dLdt_branch*dt)
                    branch_term_bed = fl.bed_function(branch_terminus/network.L0)
                    previous_branch_bed = fl.bed_function(out_dict['Termini'][-1]/network.L0)
                    previous_branch_thickness = (out_dict['Terminus_heights'][-1] - previous_branch_bed)/network.H0
                    branch_termheight = BalanceThick(branch_term_bed/network.H0, fl.Bingham_num(previous_branch_bed/network.H0, previous_branch_thickness)) + (branch_term_bed/network.H0)
                    
                branchmodel = fl.plastic_profile(startpoint=branch_terminus/network.L0, hinit=branch_termheight, endpoint=fl_amax, surf=fl.surface_function)
                if yr>dt:
                    branch_termflux = fl.icediff(profile1=out_dict[yr-dt], profile2=branchmodel)
                else:
                    branch_termflux = np.nan
                    
                out_dict[yr] = branchmodel
                out_dict['Termini'].append(branch_terminus)
                out_dict['Terminus_heights'].append(branch_termheight*network.H0)
                out_dict['Termrates'].append(dLdt_branch*dt)
                if yr > dt:
                    out_dict['Terminus_flux'].append(branch_termflux)
                else:
                    out_dict['Terminus_flux'].append(np.nan)
    
    network.model_output = model_output_dicts


#def balanceprofile(bedfunction, samplepoints, Hinit):
#    Bghm_init = flowline.Bingham_num(bedfunction(samplepoints[0]), Hinit)
#    Hy_init = BalanceThick(bedfunction(samplepoints[0]), Bghm_init)
#    Hy_list = [Hy_init]
#    
#    for j, pt in enumerate(samplepts):
#        Bghm = flowline.Bingham_num(bedfunction(pt), Hy_list[j]) #Hy_array initialized already at j=0
#        Hy = BalanceThick(bedfunction(pt), Bghm)
#        Hy_list.append(Hy)
#    
#    Hy_array = np.asarray(Hy_list)[1::]
#    return Hy_array
#    
#
#def solveL(flowline, profile, terminus_a_dot, Hy_profile, dxHy=None): #PICK UP HERE 1 JUNE
#    """Hy_profile should be an array resulting from balanceprofile
#    """
#    
#    xmin = min(profile[0])
#    xmax = max(profile[0])
#    dx = mean(diff(profile[0]))
#    L = xmax-xmin
#    
#    if dxHy is None:
#        dxHy = dx
#    
#    if dL is None:
#        dL = 5/flowline.L0
#            
#    dHdL = flowline.find_dHdL(profile, dL)
#    
#    if has_smb and terminus_balance is not None:
#        terminus_a_dot = terminus_balance
#    else:
#        terminus_a_dot = alpha_dot #setting value of surface mass balance at terminus = spatially averaged value.  Note this is likely to give strongly underestimated rates of advance/retreat
#
#    
#    #Diffs - using built-ins
#    dHdx_terminus = np.gradient(profile[1], varargs=dx)[0]
#    dHydx_terminus = np.gradient(Hy_profile, varargs=dxHy)[0]
#    
#    scipy.integrate.solve_ivp()

#
#def analyse_dLdt(flowline, profile, alpha_dot, rate_factor=1.7E-24, dL=None, debug_mode=False, has_smb=False, terminus_balance=None, submarine_melt=0):
#    """Function to compute terminus rate of advance/retreat given a mass balance forcing, a_dot.
#    Input:
#        profile: a plastic profile output from Flowline.plasticprofile of the current time step
#        alpha_dot: net rate of ice accumulation/loss.  Should be expressed in m/a /H0. Spatially averaged over whole catchment for now
#        rate_factor: flow rate factor A, assumed 3.5x10^(-25) Pa^-3 s^-1 for T=-10C based on Cuffey & Paterson
#    Optional inputs:
#        dL: passed to Flowline.find_dHdL as length of step over which to test dHdL profiles.  Default 5 m.
#        debug_mode: if True, turns on output for inspection with debugging.  Default False.
#        has_smb: describes whether we have surface mass balance data for this flowline.  Default False.  If true, will run with terminus_balance input in next argument
#        terminus_balance: value of surface mass balance (m/a /H0) at terminus, if available.  For now this is constant at the initial value.
#        submarine_melt: default 0 m/a (/L0).  If input is a value >0, dLdt applies that value annual-average submarine melt at the terminus, scaled by the portion of the terminus in contact with ocean water.
#    Returns dLdt in nondimensional units.  Multiply by L0 to get units of m/a (while T0=1a).
#    """      
#    xmin = min(profile[0])
#    xmax = max(profile[0])
#    L = xmax-xmin #length of the current profile, nondimensional
#    
#    if dL is None:
#        dL=5/flowline.L0 #step over which to test dHdL profiles
#    
#    dHdL = flowline.find_dHdL(profile, dL)
#    
#    if has_smb and terminus_balance is not None:
#        terminus_a_dot = terminus_balance
#    else:
#        terminus_a_dot = alpha_dot #setting value of surface mass balance at terminus = spatially averaged value.  Note this is likely to give strongly underestimated rates of advance/retreat
#        
#    
#    #Nondimensionalising rate factor
#    #inverse_units_of_A = flowline.T_0 * (flowline.rho_ice **3)*(flowline.g **3) * (flowline.H0 **6) / (flowline.L0 **3)
#    #s_per_annum = 31557600 #convert to time scale T0 = 1 a
#    #units_of_A = (flowline.L0 **3)/ (T_0*(flowline.rho_ice **3)*(flowline.g **3) *(flowline.H0 **6))
#    #nondim_A = rate_factor * inverse_units_of_A
#    
#    #Terminus quantities
#    SE_terminus = profile[1][0] #terminus at [0], not [-1]--may return errors if running from head downstream, but this is for terminus forcing anyway
#    Bed_terminus = profile[2][0]
#    print 'Bed_terminus={}'.format(Bed_terminus)
#    H_terminus = SE_terminus - Bed_terminus 
#    Bghm_terminus = flowline.Bingham_num(Bed_terminus, H_terminus)
#    Hy_terminus = BalanceThick(Bed_terminus, Bghm_terminus)
#
#    #Quantities at adjacent grid point
#    SE_adj = profile[1][1]
#    Bed_adj = profile[2][1]
#    H_adj = SE_adj - Bed_adj
#    Bghm_adj = flowline.Bingham_num(Bed_adj, H_adj)
#    Hy_adj = BalanceThick(Bed_adj, Bghm_adj)
#    
#    #Diffs
#    dx_term = abs(profile[0][1] - profile[0][0]) #should be ~2m in physical units
#    dHdx = (H_adj-H_terminus)/dx_term
#    dHydx = (Hy_adj-Hy_terminus)/dx_term
#    tau = flowline.Bingham_num(Bed_terminus, H_terminus) * (flowline.rho_ice * flowline.g * flowline.H0**2 / flowline.L0) #using Bingham_num handles whether tau_y constant or variable for selected flowline
#    dUdx_terminus = rate_factor * tau**3 #+1 due to sign convention with x increasing upstream from terminus and L increasing downstream
#    #nondim_dUdx_terminus = dUdx_terminus * inverse_units_of_A / ((flowline.rho_ice * flowline.g * flowline.H0**2 / flowline.L0)**3) #divide out units to get nondimensional quantity
#    nondim_dUdx_terminus = flowline.T_0 * dUdx_terminus #multiply by T_0 for nondimensionalisation
#    #nondim_dUdx_terminus = 0 #testing quasi-rigid hypothesis
#
#    Area_int = quad(dHdL, xmin, xmax)[0]
#    #print 'dH/dL at terminus = {}'.format(dHdL(xmin))
#    multiplier = 1 - (Area_int/H_terminus)
#    
#    denom = dHydx - dHdx* (1 - (Area_int/H_terminus))
#    numerator = terminus_a_dot - nondim_dUdx_terminus*H_terminus + (alpha_dot*L*dHdx/H_terminus)
#    
#    dLdt_viscoplastic = numerator/denom
#    
#    if Bed_terminus <0:
#        waterdepth = abs(Bed_terminus)
#    else:
#        waterdepth = 0
#    submarine_iceloss = submarine_melt * (waterdepth/H_terminus)
#    
#    result = dLdt_viscoplastic - submarine_iceloss
#    dimensional_result = (ref_line.L0 / ref_line.T_0) * result
#    print 'Dimensional dLdt = {} m/a'.format(dimensional_result)
#    
#    if debug_mode:
#        print 'For inspection on debugging:'
#        print 'L={}'.format(L)
#        print 'SE_terminus={}'.format(SE_terminus)
#        print 'Bed_terminus={}'.format(Bed_terminus)
#        print 'Hy_terminus={}'.format(Hy_terminus)
#        print 'dx_term={}'.format(dx_term)
#        print 'Area_int={}'.format(Area_int)
#        print 'Checking dLdt: terminus_a_dot = {}. \n H dUdx = {}. \n Ub dHdx = {}.'.format(terminus_a_dot, nondim_dUdx_terminus*H_terminus, alpha_dot*L*dHdx/H_terminus) 
#        print 'Denom: dHydx = {} \n dHdx = {} \n (1-Area_int/H_terminus) = {}'.format(dHydx, dHdx, multiplier)
#        print 'Viscoplastic dLdt={}'.format(dLdt_viscoplastic)
#        print 'Submarine ice loss = {}'.format(submarine_iceloss)
#    else:
#        pass
#
#    return result, dLdt_viscoplastic, numerator, dHydx, dHdx, multiplier, submarine_iceloss, nondim_dUdx_terminus*H_terminus
#    
##analyse dLdt in sample run for Jakobshavn
#ref_line = Jakobshavn_main.flowlines[0]
#ref_surface = ref_line.ref_profile
#ref_amax = 0.5*ref_line.length
#testpositions = linspace(0, ref_amax, num=20)
#terminus_SE = [ref_surface(0)]
#
#total_dLdt = []
#viscoplastic_dLdt = []
#numerators = []
#dHydxs = []
#dHdxs = []
#denomfactor = []
#submarine_loss = []
#dUdxH = []
#
#print 'Testing with absolute value of balance velocity, changed dUdx scaling, denser sampling'   
#for j,x in enumerate(testpositions):
#    new_term_bed = ref_line.bed_function(x/ref_line.L0)
#    previous_bed = ref_line.bed_function(testpositions[j-1]/ref_line.L0)
#    previous_thickness = (terminus_SE[-1] - previous_bed)/ref_line.H0 #nondimensional thickness for use in Bingham number
#    new_termheight = BalanceThick(new_term_bed/Jakobshavn_main.H0, ref_line.Bingham_num(previous_bed/ref_line.H0, previous_thickness)) + (new_term_bed/ref_line.H0)
#    prof = ref_line.plastic_profile(startpoint=x, hinit = new_termheight, endpoint=ref_amax, surf=ref_surface)
#    
#    dLdt_quantities = analyse_dLdt(flowline=ref_line, profile=prof, alpha_dot=Jakobshavn_main.sec_alphadot*ref_line.T_0/ref_line.H0, dL=1/L0, has_smb=True, terminus_balance=Jakobshavn_main.terminus_sec*ref_line.T_0/ref_line.H0, submarine_melt=100*ref_line.T_0/ref_line.L0, debug_mode=True)
#    
#    total_dLdt.append(dLdt_quantities[0])
#    viscoplastic_dLdt.append(dLdt_quantities[1])
#    numerators.append(dLdt_quantities[2])
#    dHydxs.append(dLdt_quantities[3])
#    dHdxs.append(dLdt_quantities[4])
#    denomfactor.append(dLdt_quantities[5])
#    submarine_loss.append(dLdt_quantities[6])
#    dUdxH.append(dLdt_quantities[7])
#    terminus_SE.append(new_termheight)
#
#total_dLdt_ma = (ref_line.L0 /ref_line.T_0) * np.array(total_dLdt)


def analyse_dLdt_dimensional(flowline, profile, alpha_dot, rate_factor=3.5E-25, dL=None, debug_mode=False, has_smb=False, terminus_balance=None, submarine_melt=0):
    """Function to compute terminus rate of advance/retreat given a mass balance forcing, a_dot.
    Input:
        profile: a plastic profile output from Flowline.plasticprofile of the current time step
        alpha_dot: net rate of ice accumulation/loss.  Should be expressed in m/a /H0. Spatially averaged over whole catchment for now
        rate_factor: flow rate factor A, assumed 3.5x10^(-25) Pa^-3 s^-1 for T=-10C based on Cuffey & Paterson
    Optional inputs:
        dL: passed to Flowline.find_dHdL as length of step over which to test dHdL profiles.  Default 5 m.
        debug_mode: if True, turns on output for inspection with debugging.  Default False.
        has_smb: describes whether we have surface mass balance data for this flowline.  Default False.  If true, will run with terminus_balance input in next argument
        terminus_balance: value of surface mass balance (m/a /H0) at terminus, if available.  For now this is constant at the initial value.
        submarine_melt: default 0 m/a (/L0).  If input is a value >0, dLdt applies that value annual-average submarine melt at the terminus, scaled by the portion of the terminus in contact with ocean water.
    Returns dLdt in nondimensional units.  Multiply by L0 to get units of m/a (while T0=1a).
    """      
    xmin = min(profile[0])*flowline.L0
    xmax = max(profile[0])*flowline.L0
    L = xmax-xmin #length of the current profile, nondimensional
    
    if dL is None:
        dL=5/flowline.L0 #step over which to test dHdL profiles
    
    dHdL = flowline.find_dHdL(profile, dL)
    
    if has_smb and terminus_balance is not None:
        terminus_a_dot = terminus_balance
    else:
        terminus_a_dot = alpha_dot #setting value of surface mass balance at terminus = spatially averaged value.  Note this is likely to give strongly underestimated rates of advance/retreat
        
    
    #Nondimensionalising rate factor
    #inverse_units_of_A = flowline.T_0 * (flowline.rho_ice **3)*(flowline.g **3) * (flowline.H0 **6) / (flowline.L0 **3)
    s_per_annum = 31557600 #convert to time scale T0 = 1 a
    #units_of_A = (flowline.L0 **3)/ (T_0*(flowline.rho_ice **3)*(flowline.g **3) *(flowline.H0 **6))
    #nondim_A = rate_factor * inverse_units_of_A
    
    #Terminus quantities
    SE_terminus = profile[1][0] *flowline.H0 #terminus at [0], not [-1]--may return errors if running from head downstream, but this is for terminus forcing anyway
    Bed_terminus = profile[2][0] *flowline.H0
    print 'Bed_terminus={}'.format(Bed_terminus)
    H_terminus = SE_terminus - Bed_terminus 
    Bghm_terminus = flowline.Bingham_num(Bed_terminus/flowline.H0, H_terminus/flowline.H0)
    Hy_terminus = BalanceThick(Bed_terminus/flowline.H0, Bghm_terminus) *flowline.H0

    #Quantities at adjacent grid point
    SE_adj = profile[1][1] *flowline.H0
    Bed_adj = profile[2][1] * flowline.H0
    H_adj = SE_adj - Bed_adj
    Bghm_adj = flowline.Bingham_num(Bed_adj/flowline.H0, H_adj/flowline.H0)
    Hy_adj = BalanceThick(Bed_adj/flowline.H0, Bghm_adj) * flowline.H0
    
    #Diffs
    dx_term = abs(profile[0][1] - profile[0][0]) *flowline.L0 #should be ~2m in physical units
    dHdx = (H_adj-H_terminus)/dx_term
    dHydx = (Hy_adj-Hy_terminus)/dx_term
    tau = Bghm_terminus * (flowline.rho_ice * flowline.g * flowline.H0**2 / flowline.L0) #using Bingham_num handles whether tau_y constant or variable for selected flowline
    dUdx_terminus = rate_factor * tau**3 # +1 due to sign convention with x increasing upstream from terminus and L increasing downstream
    dUdx_terminus_perannum = s_per_annum * dUdx_terminus #multiply by s_per_annum for consistent units

    Area_int = quad(dHdL, xmin, xmax)[0] * flowline.H0 # H0/L0 for dHdL, times L0 for integrating dx
    #print 'dH/dL at terminus = {}'.format(dHdL(xmin))
    multiplier = 1 - (Area_int/H_terminus)
    
    denom = dHydx - dHdx* (1 - (Area_int/H_terminus))
    numerator = terminus_a_dot - dUdx_terminus_perannum*H_terminus + (alpha_dot*L*dHdx/H_terminus)
    
    dLdt_viscoplastic = numerator/denom
    
    if Bed_terminus <0:
        waterdepth = abs(Bed_terminus)
    else:
        waterdepth = 0
    submarine_iceloss = submarine_melt * (waterdepth/H_terminus)
    
    result = dLdt_viscoplastic - submarine_iceloss
    nondimensional_result = result / flowline.L0
    print 'Dimensional dLdt = {} m/a'.format(result)
    print 'Nondimensional dLdt = {}'.format(nondimensional_result)
    
    if debug_mode:
        print 'For inspection on debugging - all should be DIMENSIONAL (m/a):'
        print 'L={}'.format(L)
        print 'SE_terminus={}'.format(SE_terminus)
        print 'Bed_terminus={}'.format(Bed_terminus)
        print 'Hy_terminus={}'.format(Hy_terminus)
        print 'dx_term={}'.format(dx_term)
        print 'Area_int={}'.format(Area_int)
        print 'Checking dLdt: terminus_a_dot = {}. \n H dUdx = {}. \n Ub dHdx = {}.'.format(terminus_a_dot, dUdx_terminus_perannum*H_terminus, alpha_dot*L*dHdx/H_terminus) 
        print 'Denom: dHydx = {} \n dHdx = {} \n (1-Area_int/H_terminus) = {}'.format(dHydx, dHdx, multiplier)
        print 'Viscoplastic dLdt={}'.format(dLdt_viscoplastic)
        print 'Submarine ice loss = {}'.format(submarine_iceloss)
    else:
        pass

    return result, dLdt_viscoplastic, terminus_a_dot, dUdx_terminus_perannum*H_terminus, alpha_dot*L*dHdx/H_terminus, dHydx, dHdx, multiplier, submarine_iceloss

ref_line = Jakobshavn_main.flowlines[0]
ref_surface = ref_line.ref_profile
ref_amax = 0.5*ref_line.length
testpositions = linspace(0, ref_amax, num=50)
terminus_SE = [ref_surface(0)]

total_dLdt = []
viscoplastic_dLdt = []
terminus_adots = []
dUdxH = []
UbdHdx = []
dHydxs = []
dHdxs = []
denomfactor = []
submarine_loss = []



print 'Testing with negative balance velocity, dimensional units'   
for j,x in enumerate(testpositions):
    new_term_bed = ref_line.bed_function(x/ref_line.L0)
    previous_bed = ref_line.bed_function(testpositions[j-1]/ref_line.L0)
    previous_thickness = (terminus_SE[-1] - previous_bed)/ref_line.H0 #nondimensional thickness for use in Bingham number
    new_termheight = BalanceThick(new_term_bed/Jakobshavn_main.H0, ref_line.Bingham_num(previous_bed/ref_line.H0, previous_thickness)) + (new_term_bed/ref_line.H0)
    prof = ref_line.plastic_profile(startpoint=x, hinit = new_termheight, endpoint=ref_amax, surf=ref_surface)
    
    dLdt_quantities = analyse_dLdt_dimensional(flowline=ref_line, profile=prof, alpha_dot=Jakobshavn_main.sec_alphadot, dL=1/L0, has_smb=True, terminus_balance=Jakobshavn_main.terminus_sec, submarine_melt=100, debug_mode=True)
    
    total_dLdt.append(dLdt_quantities[0])
    viscoplastic_dLdt.append(dLdt_quantities[1])
    terminus_adots.append(dLdt_quantities[2])
    dUdxH.append(dLdt_quantities[3])
    UbdHdx.append(dLdt_quantities[4])
    dHydxs.append(dLdt_quantities[5])
    dHdxs.append(dLdt_quantities[6])
    denomfactor.append(dLdt_quantities[7])
    submarine_loss.append(dLdt_quantities[8])
    terminus_SE.append(new_termheight)
#
#total_dLdt_ma = (ref_line.L0 /ref_line.T_0) * np.array(total_dLdt)


