## Meta-class of plastic glacier, which should have flowlines as classes of their own (because this is what we run model on)
## 29 Jan 2017  EHU

import warnings
import cPickle as pickle
from scipy.integrate import quad
from plastic_utilities_v2 import *
from GL_model_tools import *

class Ice(object):
    """Holds constants used in plastic model, set to default values but adjustable.
    Default values:
        H0 = 1e3, characteristic height for nondimensionalisation (in m)
        L0 = 10e3, characteristic length for nondimensionalisation (in m)
        g = 9.8 accel. due to gravity (m/s^2)
        rho_ice = 920.0, ice density (kg/m^3)
        rho_sea = 1020.0, seawater density (kg/m^3)
        default_Ty = 150e3, reasonable yield strength of ice based on lab range, Greenland obs, and Alaska optimisation (Pa)
        default_T0 = 130e3, reasonable strength of ice for Mohr-Coulomb failure criterion, based on Alaska optimisation (Pa)
    """
    def __init__(self, H0=1e3, L0=10e3, T_0=31557600, g=9.8, rho_ice=920.0, rho_sea=1020.0, default_Ty=150e3, default_T0=130e3):
        self.H0 = H0 #characteristic height for nondimensionalisation 
        self.L0 = L0
        self.T_0 = T_0 #1 annum = 31557600 s, default characteristic time.  Not to be confused with default_T0, a default Mohr-Coulumb stress
        self.g = g
        self.rho_ice = rho_ice #kg/m^3
        self.rho_sea = rho_sea #kg/m^3
        self.default_Ty = default_Ty #Pa
        self.default_T0 = default_T0 #Pa

class Branch(Ice):
    """A glacier branch that doesn't necessarily reach all the way to the terminus.  Essentially a container for coordinates.
    
    Attributes:
        coords: A list of x-y coordinates describing the line, ordered from terminus (lowest point) to head, expressed in meters
        index: An integer that identifies which flowline this is in the glacier
        order: An integer describing whether this is 
            0-> a main branch (reaches terminus)
            1-> a branch that flows directly to main branch
            2-> a branch that flows into a 1st-order branch, etc.
        flows_to: name of Branch or Flowline instance that receives this branch, if applicable
        intersect: point of intersection, if known
        width: array of glacier widths, if 
    """
    def __init__(self, coords, index=None, order=0, flows_to=None, intersect=None, has_width=True):
        Ice.__init__(self)
        self.full_coords = np.asarray(coords) #keeping width with other coords so it can be used in make_full_lines
        self.coords = np.asarray(coords)[:,0:2] # (x,y) coordinates of triples (x,y,width) read in by Flowline_CSV
        self.index = index
        self.order = order
        self.flows_to = flows_to
        self.intersect = intersect 
        if has_width:
            self.width = np.asarray(coords)[:,2]
        else:
            self.width = None
     
                     

class Flowline(Ice):
    """Holds all the information for a given flowline.  Inherits default constants from Ice class.
    For use in PlasticNetwork, all Flowlines should reach all the way to the terminus.  
    If you have a *branch* (i.e. a list of coordinates that intersects with another flowline upstream) then you should run TopToTerm first to make full-length lines.
    
    Attributes:
        coords: A list of x-y coordinates describing the line, ordered from terminus to head, expressed in meters
        index: An integer that identifies which flowline this is in the glacier
        name: name of the glacier or the particular flowline, if applicable
        initial_terminus: x-y coordinates of initial terminus position
        
    A complete instance should also have: 
        Callable bed function (make with Flowline.process_bed)
        Callable surface function
        Callable ice thickness function
        Value and type of best-fit yield strength
    """
    
    def __init__(self, coords, index=None, name=None, initial_terminus=(0,0), intersections=None, has_width=True, width_array=None):
        Ice.__init__(self)
        self.full_coords = np.asarray(coords)
        self.coords = np.asarray(coords)[:,0:2]
        self.length = ArcArray(self.coords)[-1]
        if index is None:
            self.index = np.nan
        else:
            self.index = index
        if name is None:
            self.name = 'Glacier'
        else:
            self.name = name
        self.intersections = intersections #should be set by make_full_lines in PlasticNetwork
        if has_width:
            try:
                self.width = np.absolute(np.asarray(coords)[:,2])
            except IndexError:
                self.width = np.absolute(width_array)
        else:
            self.width = None
        #if flows_to is not None:
        #    self.flows_to = flows_to
            
    
    def process_bed(self, B_field, extra_smoothing=True, degree=3, flowline_spacing=10):
        """Make callable bed function for this flowline.  B_field should be 2d-interpolated.
        degree: degree of spline fit for Univariate_Spline bed interpolation
        flowline_spacing: how closely to sample data along flowline before creating function.  Default 10m
        """
        if extra_smoothing: #turn on extra smoothing to use spline interp instead of linear.  Useful if you have tributaries joining with big differences in bed
            
            bettercoords = Evenspace(self.coords, spacing=50) #using plastic_utilities_v2 to get evenly spaced coords before interp
            arc = ArcArray(bettercoords)
            vals = [B_field(coord[0], coord[1]) for coord in bettercoords]
            sqzd = np.squeeze(vals) #need to alter shape to get into 1d interpolation
            funct=interpolate.UnivariateSpline(x=arc, y=sqzd, k=degree)
            self.bed_function = funct  
        else:
            self.bed_function = FlowProcess(self.coords, B_field)
    
    def process_surface(self, S_field):
        """Make callable surface function for this flowline.  S_field should be 2d-interpolated."""
        self.surface_function = FlowProcess(self.coords, S_field)
    
    def process_thickness(self, H_field):
        """Make callable ice thickness function along this flowline.  H_field should be 2d-interpolated."""
        self.thickness_function = FlowProcess(self.coords, H_field)
    
    def process_width(self):
        self.width_function = interpolate.interp1d(ArcArray(self.coords), self.width)
    
    def optimize_yield_strength(self, testrange=np.arange(50e3, 500e3, 5e3), arcmax=None, use_balancethick=True, initial_termpos=0, allow_upstream_breakage=False):
        """Run optimization and set the result to be the optimal value for the flowline instance.  
        Arcmax over which to run optimization can be adjusted according to how high up we think plastic approximation should go.
        initial_termpos defaults to 0 (i.e. the end of the flowline we've saved) but if observations against which we optimize have a different terminus, we can specify where to start the test profiles
        allow_upstream_breakage: argument passed to plasticmodel_error to allow or suppress upstream yield-strength breakage while optimizing
        NOTE: arcmax default is None but will set to self.length.  self unavailable at function define-time
        """
        
        #Fixing default that could not be used in definition of arguments
        if arcmax is None:
            arcmax = self.length
        CV_const_arr = []
        CV_var_arr = []
        
        bedf = self.bed_function
        surf = self.surface_function
        
        for tau in testrange:
            if use_balancethick:
                hinit = BalanceThick(bedf(initial_termpos)/H0, tau/(rho_ice*g*H0**2/L0))+(bedf(initial_termpos)/H0)
            else:
                hinit = surf(initial_termpos)/H0
                
            ##CONSTANT YIELD
            model_const = plasticmodel_error(bedf, tau, B_const, initial_termpos, hinit, arcmax, 25000, surf, allow_upstream_breakage=allow_upstream_breakage) #prescribed terminus thickness
            CV_const = model_const[1]
            ##VARIABLE YIELD
            model_var = plasticmodel_error(bedf, tau, B_var, initial_termpos, hinit, arcmax, 25000, surf, allow_upstream_breakage=allow_upstream_breakage) #prescribed terminus thickness
            CV_var = model_var[1]
            
            CV_const_arr.append(CV_const)
            CV_var_arr.append(CV_var)
        
        ##Masking 0 and NaN so that broken runs do not appear to be "optimal"
        CV_const_ma = np.ma.masked_equal(CV_const_arr, 0.0, copy=False)
        CV_var_ma = np.ma.masked_equal(CV_var_arr, 0.0, copy=False)
        
        constopt_index = np.argmin(CV_const_ma)
        varopt_index = np.argmin(CV_var_ma)
        
        constopt = testrange[constopt_index]
        varopt = testrange[varopt_index]
        
        if np.min(CV_const_ma)<np.min(CV_var_ma):
            yieldtype = 'constant'
            t_opt = constopt
        elif np.min(CV_const_ma)>np.min(CV_var_ma):
            yieldtype = 'variable'
            t_opt = varopt
        #elif np.min(CV_const_arr)==np.min(CV_var_arr): # This is rare enough that it probably means something else is wrong.
        #    yieldtype = 'variable'
        #    t_opt = varopt
        else:
            warnings.warn('Did not find minimum CV_RMS; cannot define optimal yield type.')
            
        self.yield_type = yieldtype
        self.optimal_tau = t_opt
        self.tau_const = constopt
        self.tau_var = varopt

    def Bingham_num(self, elev, thick, pos=None, time=None):
        """Defining the correct Bingham number function for plastic profile along flowline
        elev: nondimensional bed elevation
        thick: nondimensional ice thickness"""
            
        if self.yield_type is 'variable':
            #B_var from plastic_utilities, accounts for water pressure
            if elev<0:
                D = -elev #Water depth D the nondim bed topography value when Z<0
            else: 
                D = 0
            N = rho_ice*g*thick*H0 - rho_sea*g*D*H0
            mu = 0.01 #coefficient between 0 and 1, scales influence of water pressure
            tau_y = self.optimal_tau + mu*N
            return tau_y/(rho_ice*g*H0**2/L0)
        
        else: #do this for constant and as default
            return self.optimal_tau/(rho_ice*g*H0**2/L0) #B_const from plastic_utilities
            
    def plastic_profile(self, bedf=None, startpoint=0, hinit=None, endpoint=None, Npoints=25000, surf=None):
        """
        This snippet translates functionality of plastic_utilities_v2.PlasticProfile
        """
        
        H0 = self.H0
        L0 = self.L0
        
        if bedf is None:
            bedf = self.bed_function
        if hinit is None:
            hinit = self.surface_function(0)/H0
        if endpoint is None:
            endpoint = self.length
        if surf is None:
            surf = self.surface_function
        
        horiz = linspace(startpoint, endpoint, Npoints)
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
        thickarr.append(hinit-(bedf(startpoint)/H0))
        basearr.append(bedf(startpoint)/H0)
        obsarr.append((surf(startpoint))/H0)
        
        for x in horiz[1::]:
            bed = bedf(x)/H0  # value of interpolated bed function
            obsheight = (surf(x))/H0
            modelthick = thickarr[-1]
            B = self.Bingham_num(bed, modelthick, None, None)
            #Break statements for thinning below yield/water balance, or flotation
            if modelthick<BalanceThick(bed,B) and dx<0:
                print 'Thinned below water balance at x=' + str(10*x)+'km'
                break
            if modelthick<FlotationThick(bed) and dx<0: #should only break for flotation when running downstream in simulations 
                ##CONFIRM THAT THERE ARE OTHER WAYS TO CATCH FLOATING TERMINI 
                print 'Thinned below flotation at x=' + str(10*x)+'km'
                break
            else:
                basearr.append(bed)
                SEarr.append(SEarr[-1]+(B/modelthick)*dx) 
                thickarr.append(SEarr[-1]-basearr[-1])
                obsarr.append(obsheight)
        
        return (horiz[0:len(SEarr)], SEarr, basearr, obsarr)        

    def set_ref_profile(self):
        """Creates 1d interpolated function of reference profile height along flowline.
        Assigns function to be self.ref_profile  
        Give argument of nondim arclength, returns nondim surface elevation.  Useful for initial condition in thinning."""
        #plasticref = PlasticProfile(self.bed_function, self.Bingham_num, 0, self.surface_function(0)/H0, self.length, 10000, self.surface_function)
        plasticref = self.plastic_profile()
        ref_prof_interpolated = interpolate.interp1d(plasticref[0], plasticref[1], kind='linear', copy=True)
        self.ref_profile = ref_prof_interpolated
    
    def icediff(self, profile1, profile2, upstream_lim=None, shape='trapezoid', separation_buffer=None):
        """Calculate net ice loss due to calving between two plastic profiles (interpret as successive timesteps)
        Inputs:
            profile1: an output from Flowline.plastic_profile
            profile2: an output from Flowline.plastic_profile (at a later time step)
            shape: 'trapezoid' ##Add others later e.g. parabolic bed
            separation_buffer: distance away from intersection point (0D) to start counting flux from tributaries.  Should be ~0.5*width of downstream branch to prevent double-counting of ice mass.  Default 5000 m/L0. Does not affect single-branch networks.
        Output:
            dM, ice mass change at the terminus (due to calving flux) between profile 1 and profile 2, in kg
        """
        try:
            interpolated_func1 = interpolate.interp1d(profile1[0], profile1[1], kind='linear', copy=True, bounds_error=False) #Creating interpolated surface elevation profile
            bed_function1 = interpolate.interp1d(profile1[0], profile1[2]) #useful to calculate full-thickness volume change from retreat/advance
        except ValueError: #this happens when terminus retreats past upstream forcing point
            interpolated_func1 = lambda x: np.nan
            bed_function1 = lambda x: np.nan
        try:
            interpolated_func2 = interpolate.interp1d(profile2[0], profile2[1], kind='linear', copy=True, bounds_error=False) #Creating interpolated profile if needed
        except ValueError: #this happens when terminus retreats past upstream forcing point
            interpolated_func2 = lambda x: np.nan #return NaN once no longer calculating meaningful changes
    
        arcmax = self.length
                
        if upstream_lim is None:
            upstream_lim = arcmax
        upstream_limit = min(upstream_lim, arcmax) #catches if a uniform upstream limit has been applied that exceeds length of a particular flowline
        
        if separation_buffer is None:
            separation_buffer = 5000/self.L0
        
        if self.index != 0:
            idx = self.intersections[1] #for now will only catch one intersection...can refine later for more complex networks
            downstream_limit = ArcArray(self.coords)[idx] + separation_buffer
        else:
            downstream_limit = 0
    
        ## Limits of integration
        x1 = min(profile1[0]) #initial terminus position, in nondimensional units
        x2 = min(profile2[0]) #new terminus position
        print x1,x2
        #print 'Downstream limit={}.'.format(downstream_limit)
        
        dX = (x2-x1)*self.L0 #in physical units of m
        print 'dX={} m'.format(dX)
        
        w1 = self.width_function(x1)
        w2 = self.width_function(x2) ## BRANCH SEPARATION?
        #dW = abs(w2-w1) #in physical units of m
        print 'Widths: {}, {}'.format(w1, w2)
        print 'dW={} m'.format(w2-w1)
        
        if dX < 0: #second profile is more advanced than first
            bed_function2 = interpolate.interp1d(profile2[0], profile2[2]) #bed function that extends to more advanced terminus
            full_thickness = lambda x: self.H0*interpolated_func2(x) - self.H0*bed_function2(x)
        else:
            full_thickness = lambda x: self.H0*interpolated_func1(x) - self.H0*bed_function1(x)
        
        if x1 < downstream_limit:
            if x2 <= downstream_limit:
                frontal_dH = 0
            if x2 > downstream_limit:
                frontal_dH = quad(full_thickness, downstream_limit, x2)[0]
                #print 'frontal_dH integrand: \n full_thickness at downstream end: {} \n upstream end: {}'.format(full_thickness(downstream_limit), full_thickness(x2))
                #print 'Integrated: frontal_dH = {}'.format(frontal_dH)
        else:    
            frontal_dH = quad(full_thickness, x1, x2)[0] #should handle signs correctly 
            #print 'frontal_dH integrand: \n full_thickness at downstream end: {} \n upstream end: {}'.format(full_thickness(x1), full_thickness(x2))
            #print 'Integrated: frontal_dH = {}'.format(frontal_dH)
        frontal_dV = frontal_dH * self.L0 * 0.5*(w1+w2)
        frontal_dM = frontal_dV * self.rho_ice
        
        upstream_dH = lambda x: self.H0*(interpolated_func1(x) - interpolated_func2(x))*self.width_function(x)
        if x2 < downstream_limit:
            upstream_dV_raw = quad(upstream_dH, downstream_limit, upstream_limit)[0]
            #print 'upstream_dH function: interpolated_func1(downstream_limit)={}, interpolated_func2(downstream_limit)={}, width_function(downstream_limit)={}'.format(interpolated_func1(downstream_limit), interpolated_func2(downstream_limit), self.width_function(downstream_limit))
            #print 'upstream_dV integrand: \n  upstream_dH at downstream end: {} \n upstream end: {}'.format(upstream_dH(downstream_limit), upstream_dH(upstream_limit))
            #print 'Integrated: upstream_dH = {}'.format(upstream_dV_raw)
        else:
            upstream_dV_raw = quad(upstream_dH, x2, upstream_limit)[0]
            #print 'upstream_dH function: interpolated_func1(x2)={}, interpolated_func2(x2)={}, width_function(x2)={}'.format(interpolated_func1(x2), interpolated_func2(x2), self.width_function(x2))
            #print 'upstream_dV integrand: \n  upstream_dH at downstream end: {} \n upstream end: {}'.format(upstream_dH(x2), upstream_dH(upstream_limit))
            #print 'Integrated: upstream_dH = {}'.format(upstream_dV_raw)
        upstream_dV = upstream_dV_raw * self.L0
        upstream_dM = upstream_dV * self.rho_ice
        
        dM = frontal_dM + upstream_dM
        print 'dM {} = frontal {} + upstream {}'.format(dM, frontal_dM, upstream_dM)
        return dM #in kg
        
    def icediff_SL(self, profile1, profile2, upstream_lim=None, shape='trapezoid', separation_buffer=None, nan_to_zero=True):
        """Calculate sea level contribution due to calving between two plastic profiles (interpret as successive timesteps).  Nearly identical to icediff, but excludes ice already below sea level.
        Inputs:
            profile1: an output from Flowline.plastic_profile
            profile2: an output from Flowline.plastic_profile (at a later time step)
            shape: 'trapezoid' ##Add others later e.g. parabolic bed
            separation_buffer: distance away from intersection point (0D) to start counting flux from tributaries.  Should be ~0.5*width of downstream branch to prevent double-counting of ice mass.  Default 5000 m/L0. Does not affect single-branch networks.
            nan_to_zero: determines whether to fill NaN values with 0 to support summation.  Default true.
        Output:
            dM, ice mass change at the terminus (due to calving flux) between profile 1 and profile 2, in kg
        """
        try:
            interpolated_func1 = interpolate.interp1d(profile1[0], np.asarray(profile1[1]).squeeze(axis=-1), kind='linear', copy=True, bounds_error=False) #Creating interpolated surface elevation profile
            bed_function1 = interpolate.interp1d(profile1[0], profile1[2]) #useful to calculate full-thickness volume change from retreat/advance
        except ValueError: #this happens when terminus retreats past upstream forcing point
            print 'There has been a ValueError calculating the old interpolated surface profile.'
            interpolated_func1 = lambda x: np.nan
            bed_function1 = lambda x: np.nan
        try:
            interpolated_func2 = interpolate.interp1d(profile2[0], np.asarray(profile2[1]).squeeze(axis=-1), kind='linear', copy=True, bounds_error=False) #Creating interpolated profile if needed
        except ValueError: #this happens when terminus retreats past upstream forcing point
            print 'There has been a ValueError calculating the new interpolated surface profile.'
            interpolated_func2 = lambda x: np.nan #return NaN once no longer calculating meaningful changes
    
        arcmax = self.length
                
        if upstream_lim is None:
            upstream_lim = arcmax
        upstream_limit = min(upstream_lim, arcmax) #catches if a uniform upstream limit has been applied that exceeds length of a particular flowline
        
        if separation_buffer is None:
            separation_buffer = 5000/self.L0
        
        if self.index != 0:
            idx = self.intersections[1] #for now will only catch one intersection...can refine later for more complex networks
            downstream_limit = ArcArray(self.coords)[idx] + separation_buffer
        else:
            downstream_limit = 0
    
        ## Limits of integration
        x1 = min(profile1[0]) #initial terminus position, in nondimensional units
        x2 = min(profile2[0]) #new terminus position
        print x1,x2
        #print 'Downstream limit={}.'.format(downstream_limit)
        
        dX = (x2-x1)*self.L0 #in physical units of m
        print 'dX={} m'.format(dX)
        
        w1 = self.width_function(x1)
        w2 = self.width_function(x2) ## BRANCH SEPARATION?
        #dW = abs(w2-w1) #in physical units of m
        print 'Widths: {}, {}'.format(w1, w2)
        print 'dW={} m'.format(w2-w1)
        
        if dX < 0: #second profile is more advanced than first
            bed_function2 = interpolate.interp1d(profile2[0], profile2[2]) #bed function that extends to more advanced terminus
            full_thickness = lambda x: self.H0*interpolated_func2(x) - max(0, self.H0*bed_function2(x)) #max(0,D) adjusts for ice that is already below sea level
        else:
            full_thickness = lambda x: self.H0*interpolated_func1(x) - max(0, self.H0*bed_function1(x))
        
        if x1 < downstream_limit:
            if x2 <= downstream_limit:
                frontal_dH = 0
            if x2 > downstream_limit:
                frontal_dH = quad(full_thickness, downstream_limit, x2)[0]
                print 'frontal_dH integrand: \n full_thickness at downstream end: {} \n upstream end: {}'.format(full_thickness(downstream_limit), full_thickness(x2))
                print 'Integrated: frontal_dH = {}'.format(frontal_dH)
        else:    
            frontal_dH = quad(full_thickness, x1, x2)[0] #should handle signs correctly 
            print 'frontal_dH integrand: \n full_thickness at downstream end: {} \n upstream end: {}'.format(full_thickness(x1), full_thickness(x2))
            print 'Integrated: frontal_dH = {}'.format(frontal_dH)
        frontal_dV = frontal_dH * self.L0 * 0.5*(w1+w2)
        frontal_dM = frontal_dV * self.rho_ice
        if nan_to_zero:
            frontal_dM = np.nan_to_num(frontal_dM)
        
        upstream_dH = lambda x: self.H0*(interpolated_func1(x) - interpolated_func2(x))*self.width_function(x)
        if x2 < downstream_limit:
            upstream_dV_raw = quad(upstream_dH, downstream_limit, upstream_limit)[0]
            #print 'upstream_dH function: interpolated_func1(downstream_limit)={}, interpolated_func2(downstream_limit)={}, width_function(downstream_limit)={}'.format(interpolated_func1(downstream_limit), interpolated_func2(downstream_limit), self.width_function(downstream_limit))
            print 'upstream_dV integrand: \n  upstream_dH at downstream end: {} \n upstream end: {}'.format(upstream_dH(downstream_limit), upstream_dH(upstream_limit))
            print 'Integrated: upstream_dH = {}'.format(upstream_dV_raw)
        else:
            upstream_dV_raw = quad(upstream_dH, x2, upstream_limit)[0]
            #print 'upstream_dH function: interpolated_func1(x2)={}, interpolated_func2(x2)={}, width_function(x2)={}'.format(interpolated_func1(x2), interpolated_func2(x2), self.width_function(x2))
            print 'upstream_dV integrand: \n  upstream_dH at downstream end: {} \n upstream end: {}'.format(upstream_dH(x2), upstream_dH(upstream_limit))
            print 'Integrated: upstream_dH = {}'.format(upstream_dV_raw)
        upstream_dV = upstream_dV_raw * self.L0
        upstream_dM = upstream_dV * self.rho_ice
        if nan_to_zero:
            upstream_dM = np.nan_to_num(upstream_dM)
        
        dM = frontal_dM + upstream_dM
        print 'dM {} = frontal {} + upstream {}'.format(dM, frontal_dM, upstream_dM)
        return dM #in kg
    
    def find_dHdL(self, profile, dL=None, debug_mode=False):
        """Function to compute successive profiles of length L-dL, L, L+dL to calculate dHdL over glacier flowline.
        Input: 
            profile: a plastic profile output from Flowline.plasticprofile of length L
            dL: spatial step to use in calculating dHdL.  Default 5 meters
        """
        if dL is None:
            dL = 5/self.L0 #nondimensional
        
        xmin = min(profile[0])
        xmax = max(profile[0])
        L_init = xmax-xmin
        
        if xmin-dL > 0:
            x_fwd = xmin-dL #note coord system has x=0 at the terminus, so xmin-dL is a more advanced position
            x_bk = xmin+dL
        else:
            x_fwd = xmin
            x_bk = xmin + 2*dL
        if xmin+dL > xmax:
            x_fwd = xmin - 2*dL
            x_bk = xmin
        
        
        #Terminus quantities
        SE_terminus = profile[1][0] #CONFIRM that terminus is at [0] and not [-1]
        Bed_terminus = profile[2][0]
        H_terminus = SE_terminus - Bed_terminus 
        Bghm_terminus = self.Bingham_num(Bed_terminus, H_terminus)
        
        #Profile advanced by dL - note coord system means xmin-dL is more advanced, as x=0 is at initial terminus position
        bed_mindL = (self.bed_function(x_fwd))/self.H0
        s_mindL = BalanceThick(bed_mindL, Bghm_terminus) + bed_mindL
        profile_mindL = self.plastic_profile(startpoint=x_fwd, hinit = s_mindL, endpoint = xmax, surf = self.surface_function)
        H_mindL = np.array(profile_mindL[1]) - np.array(profile_mindL[2]) #array of ice thickness from profile
        Hx_mindL = interpolate.interp1d(profile_mindL[0].squeeze(), H_mindL.squeeze(), bounds_error=False, fill_value=0)
        
        #Profile retreated by dL
        bed_plusdL = (self.bed_function(x_bk))/self.H0
        s_plusdL = BalanceThick(bed_plusdL, Bghm_terminus) + bed_plusdL
        profile_plusdL = self.plastic_profile(startpoint = x_bk, hinit = s_plusdL, endpoint = xmax, surf=self.surface_function)
        H_plusdL = np.array(profile_plusdL[1]) - np.array(profile_plusdL[2]) #array of ice thickness
        Hx_plusdL = interpolate.interp1d(profile_plusdL[0].squeeze(), H_plusdL.squeeze(), bounds_error=False, fill_value=0)
        
        dHdLx = lambda x: (Hx_mindL(x) - Hx_plusdL(x))/(2*dL)
        
        if debug_mode:
            print 'Debugging dHdL.  Inspect:'
            print 'H_terminus={}'.format(H_terminus)
            print 'Hx_mindL={}'.format(Hx_mindL(xmin))
            print 'Hx_plusdL={}'.format(Hx_plusdL(xmin))
        else:
            pass
        
        return dHdLx
    
    def dLdt_dimensional(self, profile, alpha_dot, rate_factor=1.7E-24, dL=None, debug_mode=False, has_smb=False, terminus_balance=None, submarine_melt=0):
        """Function to compute terminus rate of advance/retreat given a mass balance forcing, a_dot.
        Input:
            profile: a plastic profile output from Flowline.plasticprofile of the current time step
            alpha_dot: net rate of ice accumulation/loss.  Should be expressed in m/a. Spatially averaged over whole catchment for now
            rate_factor: flow rate factor A, assumed 3.5x10^(-25) Pa^-3 s^-1 for T=-10C, or 1.7E-24 for T=-2C, based on Cuffey & Paterson
        Optional inputs:
            dL: passed to Flowline.find_dHdL as length of step over which to test dHdL profiles.  Default 5 m.
            debug_mode: if True, turns on output for inspection with debugging.  Default False.
            has_smb: describes whether we have surface mass balance data for this flowline.  Default False.  If true, will run with terminus_balance input in next argument
            terminus_balance: value of surface mass balance (m/a, dimensional) at terminus, if available.  For now this is constant at the initial value.
            submarine_melt: default 0 m/a (dimensional).  If input is a value >0, dLdt applies that value annual-average submarine melt at the terminus, scaled by the portion of the terminus in contact with ocean water.
        Returns dLdt in nondimensional units.  Multiply by L0 to get units of m/a (while T0=1a).
        """      
        xmin = min(profile[0])*self.L0
        xmax = max(profile[0])*self.L0
        L = xmax-xmin #length of the current profile, nondimensional
        
        if dL is None:
            dL=5/self.L0 #step over which to test dHdL profiles
        
        dHdL = self.find_dHdL(profile, dL)
        
        if has_smb and terminus_balance is not None:
            terminus_a_dot = terminus_balance
        else:
            terminus_a_dot = alpha_dot #setting value of surface mass balance at terminus = spatially averaged value.  Note this is likely to give strongly underestimated rates of advance/retreat
            
        
        #Converting rate factor to time scale 1 a
        s_per_annum = 31557600 #convert to time scale T0 = 1 a

        
        #Terminus quantities
        SE_terminus = profile[1][0] *self.H0 #terminus at [0], not [-1]--may return errors if running from head downstream, but this is for terminus forcing anyway
        Bed_terminus = profile[2][0] * self.H0
        H_terminus = SE_terminus - Bed_terminus 
        Bghm_terminus = self.Bingham_num(Bed_terminus/self.H0, H_terminus/self.H0)
        Hy_terminus = BalanceThick(Bed_terminus/self.H0, Bghm_terminus) * self.H0
    
        #Quantities at adjacent grid point
        SE_adj = profile[1][1] *self.H0
        Bed_adj = profile[2][1] *self.H0
        H_adj = SE_adj - Bed_adj
        Bghm_adj = self.Bingham_num(Bed_adj/self.H0, H_adj/self.H0)
        Hy_adj = BalanceThick(Bed_adj/self.H0, Bghm_adj) * self.H0
        
        #Diffs
        dx_term = abs(profile[0][1] - profile[0][0]) *self.L0 #should be ~2m in physical units
        dHdx = (H_adj-H_terminus)/dx_term
        dHydx = (Hy_adj-Hy_terminus)/dx_term
        tau = Bghm_terminus * (self.rho_ice * self.g * self.H0**2 / self.L0) #using Bingham_num handles whether tau_y constant or variable for selected flowline
        dUdx_terminus = s_per_annum * rate_factor * tau**3 #sign convention with x increasing upstream from terminus
    
        Area_int = quad(dHdL, xmin, xmax)[0] * self.H0 # H0/L0 for dHdL, times L0 for integrating dx
        multiplier = 1 - (Area_int/H_terminus)
        
        denom = dHydx - (dHdx* multiplier)
        numerator = terminus_a_dot - dUdx_terminus*H_terminus + (alpha_dot*L*dHdx/H_terminus)
        
        dLdt_viscoplastic = numerator/denom
        
        if Bed_terminus <0:
            waterdepth = abs(Bed_terminus)
        else:
            waterdepth = 0
        submarine_iceloss = submarine_melt * (waterdepth/H_terminus)
        
        result = dLdt_viscoplastic + submarine_iceloss #should be + due to sign convention -- new terminus position will be x_{n+1} = x_n + dLdt as x increases away from terminus
        
        if debug_mode:
            print 'For inspection on debugging - all should be DIMENSIONAL (m/a):'
            print 'L={}'.format(L)
            print 'SE_terminus={}'.format(SE_terminus)
            print 'Bed_terminus={}'.format(Bed_terminus)
            print 'Hy_terminus={}'.format(Hy_terminus)
            print 'dx_term={}'.format(dx_term)
            print 'Area_int={}'.format(Area_int)
            print 'Checking dLdt: terminus_a_dot = {}. \n H dUdx = {}. \n Ub dHdx = {}.'.format(terminus_a_dot, dUdx_terminus*H_terminus, alpha_dot*L*dHdx/H_terminus) 
            print 'Denom: dHydx = {} \n dHdx = {} \n (1-Area_int/H_terminus) = {}'.format(dHydx, dHdx, multiplier)
            print 'Viscoplastic dLdt={}'.format(dLdt_viscoplastic)
            print 'Submarine ice loss = {}'.format(submarine_iceloss)
        else:
            pass
    
        return result
  
              
class PlasticNetwork(Ice):
    """A container for a named glacier on which we want to run the plastic model.  Can take instances of Branch as arguments and make them Flowlines.
    
    Attributes:
        name: A string with what we call the glacier for modelling/analysis
        init_type: Can be 'Branch' or 'Flowline'. If we already have Flowlines we don't need to call make_full_lines
        branches: A tuple of Flowline objects describing the branches of this glacier.  MUST USE (BRANCH,) SYNTAX IF ONLY ONE BRANCH
        main_terminus: Coordinates of the terminus of the "main branch" (index 0) in the initial state
    """
    
    def __init__(self, name='Glacier', init_type='Branch', branches=[], main_terminus=(0,0)):
        Ice.__init__(self)
        self.name = name
        if init_type is 'Branch':
            self.branches = branches
            self.flowlines = []
            print 'You have initialised {} with glacier Branches.  Please call make_full_lines and process_full_lines to proceed.'.format(name)
        elif init_type is 'Flowline':
            self.flowlines = [branches]
            if not isinstance(branches, tuple):
                branches = (branches,)
            self.branches = branches
            print 'You have initialised {} with glacier Flowlines.  Please call process_full_lines to proceed with simulations.'.format(name)
        else:
            if not isinstance(branches, tuple):
                branches = (branches,)
            self.branches = branches
            self.flowlines = []
            print 'Unrecognised initialisation type.'
        self.main_terminus = main_terminus
    
    def make_full_lines(self):
        """Runs spatial KDTree search to connect flowlines at nearest points. Creates Flowline objects that reach from each branch head to the terminus.
        Note: this will overwrite any flowlines already saved to the network.
        """
        self.flowlines = []
        
        #KDTree search connects at closest points...may still be wonky
        #project somehow with Shapely?
        mainline = self.branches[0].full_coords
        maintree = spatial.KDTree(mainline[:,0:2])
        
        full_lines = {}
        j = 0
        ix = {}
        

        while j<len(self.branches): #making full length flowline for each branch
            branch = self.branches[j]
            br = branch.full_coords 
            
            if branch.intersect is not None: #allowing shortcut if intersection point already output by FilterMainTributaries
                idx = branch.intersect
                #add flows_to shortcut functionality here for branches that flow to a branch other than main 
            else: #conduct KDTree search if we have to
                pt = br[0, 0:2] #xy coordinates of an (x,y,width) triple
                dist, idx = maintree.query(pt, distance_upper_bound=5000)
            #print dist, idx
            
            if idx==len(mainline): #if branch does not intersect with the main line
                print 'Branch {} does not intersect main line.  Searching nearest trib.'.format(j)
                tribtree = spatial.KDTree(full_lines[j-1][:,0:2]) #line of nearest trib
                dist_t, idx_t = tribtree.query(pt, distance_upper_bound=1000)
                if idx_t==len(full_lines[j-1]):
                    print 'Branch {} also does not intersect tributary {}.  Appending raw line.  Use with caution.'.format(j, j-1)
                    full_lines[j] = br
                    ix[j] = None
                else:
                    tribfrag = branchlist[j-1][:idx_t]
                    fullbranch = np.concatenate((tribfrag, br))
                    full_lines[j] = fullbranch
                    ix[j] =(j-1, idx_t) #noting location of intersection with nearest tributary (if applicable)
                j+=1
            else:
                print mainline[idx]
                mainfrag = mainline[:idx]
                fullbranch = np.concatenate((mainfrag, br)) #Shapely project here
                full_lines[j] = fullbranch
                ix[j] = (0, idx) #noting location of intersection with mainline
                j+=1
        
        k = 0
        while k<len(full_lines):
            coordinates = np.asarray(full_lines[k])
            self.flowlines.append(Flowline(coordinates, index=k, intersections=ix[k]))
            k+=1
    
    def process_full_lines(self, bed_field, surface_field, thickness_field):
        """Processes bed, surface, and thickness functions for new flowlines.  Use if network initialized with Branches.
        """
        for line in self.flowlines:
            if line.index !=0:
                line.process_bed(bed_field, extra_smoothing=True, degree=3)
            else:
                line.process_bed(bed_field, extra_smoothing=False)
            line.process_surface(surface_field)
            line.process_thickness(thickness_field)
            line.process_width()
            
            if line.thickness_function(0) < FlotationThick(line.bed_function(0)): #uses FlotationThick from plastic_utilities_v2
                warnings.warn('{} line {} may have a floating terminus.  Run remove_floating and re-initialise.'.format(self.name, line.index))
    
    def remove_floating(self, test_SE=True):
        """Checks mainline for floating terminus, removes offending coordinates.
        Default arg.:
            test_SE: when set to the default value of True, function performs an additional test to look for negative surface elevation.  
                Tests SE = bed + thickness to identify ice that is either floating or a data artefact.
        """
        fltest = []
        mainline = self.flowlines[0]
        arc = ArcArray(mainline.coords)
        for pt in arc:
            thick_init = mainline.thickness_function(pt)
            bed_init = mainline.bed_function(pt)
            flthick = FlotationThick(bed_init)
            fltest.append(thick_init - flthick)
        nonfloating = argwhere(sign(fltest)>0) #all points where thickness exceeds flotation
        cleanthick = nonfloating[0] #first point where thickness exceeds flotation
        
        if test_SE:
            surftest = []
            for pt in arc:
                thick_init = mainline.thickness_function(pt)
                bed_init = mainline.bed_function(pt)
                se_test = thick_init + bed_init
                surftest.append(se_test)
            nonneg = argwhere(sign(surftest)>0)
            cleansurf = nonneg[0]
            cleanfrom = max(cleanthick, cleansurf) #index where both tests have been passed
        else:
            cleanfrom = cleanthick #if not testing surface elevation, simply use result of first test
        
        cleaned_coords = mainline.coords[int(cleanfrom)::]
        trimmed_width = mainline.width[int(cleanfrom)::]
        
        self.flowlines[0].coords = cleaned_coords
        self.flowlines[0].width = trimmed_width
        self.branches[0].coords = cleaned_coords
        self.branches[0].width = trimmed_width
        self.flowlines[0].length = ArcArray(cleaned_coords)[-1]
        
        print 'Coordinates with ice thickness less than flotation removed.  Please re-run make_full_lines and process_full_lines for network {}.'.format(self.name)
            
    
    def optimize_network_yield(self, testrange=arange(50e3, 500e3, 5e3), arcmax_list = None, check_all=False, use_balancethick=True, nw_initial_termpos=0, allow_upstream_breakage=False):
        """Runs yield strength optimization on each of the network Flowline objects and sets the yield strength of the 'main' branch as the network's default.
        Could find a better way to do this...
        arcmax_list: list of how high up to optimise in order of flowline index; default is full length (NOTE: default of "None" is set to more sensible default in function.  Same problem with self unavailable at define-time)
        check_all=False/True decides whether to do full optimisation to find tau_y on each line (True) or just take the main line value (False)
        nw_initial_termpos: default to 0, sets location of terminus from which to optimize
        """
        #Fixing default that could not be used in definition of arguments
        if arcmax_list is None:
            arcmax_list = [self.flowlines[i].length for i in range(len(self.flowlines))]
        
        try:
            self.network_yield_type = self.flowlines[0].yield_type
            self.network_tau = self.flowlines[0].optimal_tau
        except AttributeError:
            self.flowlines[0].optimize_yield_strength(testrange, arcmax=arcmax_list[0], use_balancethick=use_balancethick, initial_termpos=nw_initial_termpos, allow_upstream_breakage=allow_upstream_breakage)
            self.network_yield_type = self.flowlines[0].yield_type
            self.network_tau = self.flowlines[0].optimal_tau
        
        if check_all:
            for j,line in enumerate(self.flowlines):
                line.optimize_yield_strength(testrange, arcmax=arcmax_list[j], use_balancethick=use_balancethick, initial_termpos=nw_initial_termpos, allow_upstream_breakage=allow_upstream_breakage) #run optimisation for each flowline
        else:
            pass

    def network_ref_profiles(self, use_mainline_tau=True):
        """1D interpolated profile of reference surface along each processed flowline of the network.  Uses set_ref_profile from Flowline class.
        By default, constructs these using yield strength of main branch (so that all agree at terminus).  
        use_mainline_tau=False will force use of each line's own yield strength & type
        """
        for line in self.flowlines:
            if use_mainline_tau:
                line.yield_type = self.network_yield_type
                line.optimal_tau = self.network_tau
            else:
                pass #use the optimal values found in optimization procedure.  Might crash if didn't look for optimal values on each line when doing optimisation
                #could also insert code here to use spatially variable tau_y
            line.set_ref_profile()
                

    def network_time_evolve(self, testyears=arange(100), ref_branch_index=0, upgl_ref=15000/L0, thinrate=10/H0, thinvalues=None, upstream_limits=None, use_mainline_tau=True):
        """Time evolution on a network of Flowlines.  Lines should be already optimised and include reference profiles from network_ref_profiles 
        Arguments:
            testyears: a range of years to test, indexed by years from nominal date of ref profile (i.e. not calendar years)
            ref_branch_index: which branch to use for forcing.  Default is main branch ("0") but can be changed
            upgl_ref: where to apply upglacier thinning.  Default is 15km upstream, or top of glacier if flowline <15km
            thinrate: thinning rate (in case of constant) to apply at reference point.  Default is 10m/a, based on large Greenland outlets
        Optional args:   
            thinvals: array of the same length as testyears with how much thinning to apply in each year.
            Offers the option to define thinning as persistence of obs or other nonlinear function.
            upstream_limits: array determining where to cut off modelling on each flowline, ordered by index.  Default is full length of lines.
            use_mainline_tau=False will force use of each line's own yield strength & type
    
        returns model output as dictionary for each flowline 
        """
        
        #Fixing default values
        if upstream_limits is None:
            upstream_limits=[fl.length for fl in self.flowlines]
        
        if thinvalues is None:  
            thinvals = np.full(len(testyears), thinrate)
        else:
            thinvals = thinvalues
        
        dt = mean(diff(testyears)) #size of time step
        
        model_output_dicts = [{'Termini': [0],
        'Terminus_heights': [fl.surface_function(0)],
        'Termrates': [],
        'Terminus_flux': []
        } for fl in self.flowlines]
        
        #Setting up reference branch (default = mainline)
        ref_line = self.flowlines[ref_branch_index]
        ref_amax = upstream_limits[ref_branch_index]
        ref_surface = ref_line.surface_function
        refpt = min(ref_amax, upgl_ref) #apply forcing at top of branch if shorter than reference distance.  In general would expect forcing farther downstream to give weaker response
        refht = ref_line.ref_profile(refpt)
        if use_mainline_tau:
            ref_line.optimal_tau = self.network_tau
            ref_line.yield_type = self.network_yield_type
        else:
            pass
        refdict = model_output_dicts[ref_branch_index]
        refdict[0] = ref_line.ref_profile #initial condition for time evolution - needed to calculate calving flux at first timestep
        
        
        #Create plastic profiles on all branches of network, following methodology described in Ultee & Bassis (2017), for each time step
        for k, yr in enumerate(testyears):
            thinning = np.sum(thinvals[:k])
            
            #Forcing on reference branch (default ref branch = mainline)
            fwdmodel = ref_line.plastic_profile(startpoint=refpt, hinit = refht - thinning, endpoint=0, surf=ref_surface)
            modelterm_pos = self.L0*min(fwdmodel[0]) #in m
            modelterm_height = self.H0*(fwdmodel[1][-1]) #in m - need to confirm this is the right value in output
            full_line_model = ref_line.plastic_profile(startpoint=modelterm_pos/self.L0, hinit=modelterm_height/self.H0, endpoint=ref_amax, surf=ref_surface)
            
            dL = modelterm_pos - refdict['Termini'][-1]
            if yr > dt:
                termflux = ref_line.icediff(profile1=refdict[yr-dt], profile2=fwdmodel)
                refdict['Terminus_flux'].append(termflux)
            else:
                pass
            refdict[yr] = full_line_model #saving full modelled profile. need to change this output type
            refdict['Termini'].append(modelterm_pos)
            refdict['Terminus_heights'].append(modelterm_height) 
            refdict['Termrates'].append(dL/dt)

            
            #Run from terminus upstream for non-ref branches
            for j, fl in enumerate(self.flowlines):
                out_dict = model_output_dicts[j]
                fl_amax = upstream_limits[j]
                if j==ref_branch_index:
                    continue
                else:
                    if use_mainline_tau:
                        fl.optimal_tau = self.network_tau
                        fl.yield_type = self.network_yield_type #this is probably too powerful, but unclear how else to exploit Bingham_number functionality
                    else:
                        pass
                    branchmodel = fl.plastic_profile(startpoint=modelterm_pos/self.L0, hinit = modelterm_height/self.H0, endpoint=fl_amax, surf = fl.surface_function)
                    out_dict[yr] = branchmodel
                    if yr > dt:
                        termflux = fl.icediff(profile1=out_dict[yr-dt], profile2=branchmodel) #NEED MORE SOPHISTICATED HANDLING OF COMMON CHANNEL HERE
                        out_dict['Terminus_flux'].append(termflux)
                    else:
                        pass

        self.model_output = model_output_dicts
        
    def balance_adot(self, V_field, use_width=False, L_limit=6.0):
        """Function to compute spatially-averaged accumulation rate that balances observed terminus velocity
        Input:
            V_field: 2d-interpolated function to report magnitude of velocity (i.e. ice flow speed in m/a) given (x,y) coordinates
            use_width: whether to consider upstream width when calculating balance accumulation.  Default is no.
            L_limit: nondimensional upstream distance that indicates how far we trust our model.  Default is 6.0 (60 km in dimensional units).
        """
        terminus = self.flowlines[0].coords[0][0:2] #All lines of a network should share a terminus at initial state
        terminus_speed = V_field(terminus[0], terminus[1]) #dimensional in m/a
        terminus_width = self.flowlines[0].width_function(0) #dimensional in m
        terminus_thickness = self.H0*self.flowlines[0].ref_profile(0) - self.flowlines[0].bed_function(0) #dimensional in m
        balance_thickness = self.H0*BalanceThick(self.flowlines[0].bed_function(0)/self.H0, self.flowlines[0].Bingham_num(self.flowlines[0].bed_function(0)/self.H0, terminus_thickness))
        termflux = terminus_speed * terminus_width * balance_thickness
        
    
        sa = []
        total_L = []
        for fl in self.flowlines:
            L = min(fl.length, L_limit)
            surface_area_nondim = quad(fl.width_function, 0, L)[0]
            surface_area = self.L0 * surface_area_nondim
            sa.append(surface_area)
            total_L.append(self.L0 * L)
        catchment_area = sum(sa)
        total_ice_length = sum(total_L)
        
        if use_width:
            balance_a = termflux/catchment_area
        else:
            balance_a = terminus_speed*balance_thickness/total_ice_length
        
        self.balance_forcing = float(balance_a) #save to this network instance
        return balance_a     
    
    def terminus_time_evolve(self, testyears=arange(100), ref_branch_index=0, initial_term=0, alpha_dot=None, alpha_dot_variable=None, upstream_limits=None, 
    separation_buffer=None, use_mainline_tau=True, debug_mode=False, dt_rounding=5, dL=None, has_smb=False, terminus_balance=None, 
    submarine_melt=0, rate_factor=1.7E-24, output_heavy=False, sl_contribution=True):
        """Time evolution on a network of Flowlines, forced from terminus.  Lines should be already optimised and include reference profiles from network_ref_profiles
        Arguments:
            testyears: a range of years to test, indexed by years from nominal date of ref profile (i.e. not calendar years)
            ref_branch_index: which branch to use for forcing.  Default is main branch ("0") but can be changed
            initial_term: where to initialize, in m relative to initial terminus of main flowline.  Default is 0 (start at initial terminus)
            alpha_dot: spatially averaged accumulation rate (forcing)
        Optional args:   
            alpha_dot_variable: array of the same length as testyears with a different alpha_dot forcing to use in each year
            Offers the option to define thinning as persistence of obs or other nonlinear function.
            upstream_limits: array determining where to cut off modelling on each flowline, ordered by index.  Default is full length of lines.
            separation_buffer: determines how far away from point of intersection tributaries should start their own simulation.  Default 5000 m/H0 - prevents double-counting in troughs ~10km wide.
            use_mainline_tau=False will force use of each line's own yield strength & type
            debug_mode=True will turn on interim output for inspection
            dt_rounding: determines how many digits of dt to keep--default 5.  If your time step is less than 1 annum, you may find numerical error.  dt_rounding takes care of it. 
            dL: passed to Flowline.find_dHdL as length of step over which to test dHdL profiles.  Default 5 m.
            has_smb: passed to Flowline.dLdt_dimensional to indicate whether we have spatially varying SMB/SEC
            terminus_balance: dimensional (m/a), passed to dLdt_dimensional to indicate thinning rate at terminus
            submarine_melt: allows inclusion of submarine melt rate - set to 0 for now as its implementation is inconsistent (11 Jun 2018)
            rate_fator: passed to Flowline.dLdt_dimensional. A in Glen's flow law (Pa^-3 s^-1)
            output_heavy: determines whether to keep all profiles or just Termini, Terminus_flux summary arrays for each branch
            returns model output as dictionary for each flowline 
        """
    
        #Fixing default values
        if upstream_limits is None:
            upstream_limits=[fl.length for fl in self.flowlines]
        
        if alpha_dot is None:
            alpha_dot = 0.2 #m/a, dimensional
        
        if alpha_dot_variable is None:  
            alpha_dot_vals = np.full(len(testyears), alpha_dot)
        else:
            alpha_dot_vals = alpha_dot_variable
        
        if dL is None:
            dL = 5/self.L0
        
        if separation_buffer is None:
            separation_buffer_default = 5000/self.L0
        
        dt = round(mean(diff(testyears)), dt_rounding) #size of time step, rounded to number of digits specified by dt_rounding
        
        model_output_dicts = [{'Termini': [0],
        'Terminus_heights': [fl.surface_function(0)],
        'Termrates': [],
        'Terminus_flux': [],
        'testyears': testyears
        } for fl in self.flowlines]
        
        #Mainline reference
        ref_line = self.flowlines[ref_branch_index]
        ref_amax = upstream_limits[ref_branch_index]
        ref_surface = ref_line.ref_profile
        #refpt = min(ref_amax, upgl_ref) #apply forcing at top of branch if shorter than reference distance.  In general would expect forcing farther downstream to give weaker response
        #refht = ref_line.ref_profile(refpt)
        if use_mainline_tau:
            ref_line.optimal_tau = self.network_tau
            ref_line.yield_type = self.network_yield_type
        else:
            pass
        refdict = model_output_dicts[ref_branch_index]
        if initial_term != 0:
            initial_term_bed = ref_line.bed_function(initial_term)
            term_balance_height = (initial_term_bed/self.H0) + BalanceThick(initial_term_bed/self.H0, ref_line.Bingham_num(initial_term_bed/self.H0, 0)) 
            refdict[0] = ref_line.plastic_profile(startpoint=initial_term, hinit=term_balance_height, endpoint=ref_amax, surf=ref_surface) #initial condition for time evolution - needed to calculate calving flux at first timestep
        else: #use initial obs
            refdict[0] = ref_line.plastic_profile(startpoint=0, hinit=ref_surface(0), endpoint=ref_amax, surf=ref_surface) #initial condition for time evolution - needed to calculate calving flux at first timestep
    
        
        #Assume same terminus
        for k, yr in enumerate(testyears):
            print yr
            alpha_dot_k = alpha_dot_vals[k]
            key = round(yr-dt, dt_rounding) #allows dictionary referencing when dt is < 1 a
            
            if k<1:
                dLdt_annum = ref_line.dLdt_dimensional(profile=refdict[0], alpha_dot=alpha_dot_k, debug_mode=debug_mode, dL=dL, has_smb=has_smb, terminus_balance=terminus_balance, submarine_melt=submarine_melt, rate_factor=rate_factor)
            else:
                dLdt_annum = ref_line.dLdt_dimensional(profile=refdict[key], alpha_dot=alpha_dot_k, debug_mode=debug_mode, dL=dL, has_smb=has_smb, terminus_balance=terminus_balance, submarine_melt=submarine_melt, rate_factor=rate_factor)
            if np.isnan(dLdt_annum): #happens if retreat hits edge of domain unexpectedly
                dLdt_annum = refdict['Termrates'][-1] / dt #replace with last non-nan value
            else:
                pass
            #Ref branch
    
            new_termpos_raw = refdict['Termini'][-1]+(dLdt_annum*dt) #Multiply by dt in case dt!=1 annum.  Multiply dLdt by L0 because its output is nondimensional
            new_termpos_posdef = max(0, new_termpos_raw)
            if new_termpos_posdef > (ref_amax * self.L0):
                print 'Terminus retreated past upstream limit. Resetting terminus position to = upstream limit.'
                new_termpos = refdict['Termini'][-1] #set to previous good position
                #new_termpos = ref_amax *self.L0 #glacier sits at edge of domain if terminus retreats beyond upstream limit
            else:
                new_termpos = new_termpos_posdef
            if debug_mode:
                print 'Looping. year={}.'.format(yr)
                print 'yr-dt={}.'.format(yr-dt)
                print 'key={}'.format(key)
                print 'dLdt_annum = {} m'.format(dLdt_annum)
                print 'New terminus position = {}'.format(new_termpos)
            else:
                pass
            previous_bed = ref_line.bed_function(refdict['Termini'][-1]/self.L0)
            try:
                new_term_bed = ref_line.bed_function(new_termpos/self.L0)
            except ValueError: #interpolated bed function sometimes complains that ref_amax * L0 is above the interpolation range.  Use last good value if this is the case.
                new_term_bed = previous_bed
            previous_thickness = (refdict['Terminus_heights'][-1] - previous_bed)/self.H0 #nondimensional thickness for use in Bingham number
            new_termheight = BalanceThick(new_term_bed/self.H0, ref_line.Bingham_num(previous_bed/self.H0, previous_thickness)) + (new_term_bed/self.H0)
            new_profile = ref_line.plastic_profile(startpoint=new_termpos/self.L0, hinit=new_termheight, endpoint=ref_amax, surf=ref_surface)
            if yr>dt:
                #key = round(yr-dt, dt_rounding)
                if sl_contribution:
                    termflux = ref_line.icediff_SL(profile1=refdict[key], profile2=new_profile)
                else:
                    termflux = ref_line.icediff(profile1=refdict[key], profile2=new_profile)
            else:
                termflux = np.nan
            
            new_key = round(yr, dt_rounding)    
            refdict[new_key] = new_profile
            refdict['Terminus_flux'].append(termflux)
            refdict['Termini'].append(new_termpos)
            refdict['Terminus_heights'].append(new_termheight*self.H0)
            refdict['Termrates'].append(dLdt_annum*dt)
            
            #Other branches, incl. branch splitting
            for j, fl in enumerate(self.flowlines):
                out_dict = model_output_dicts[j]
                fl_amax = upstream_limits[j]
                if j==ref_branch_index:
                    continue
                else:
                    #separation_distance = ArcArray(fl.coords)[fl.intersections[1]] + separation_buffer #where line separates from mainline
                    try:
                        separation_buffer = 0.5*ref_line.width_function(ArcArray(fl.coords)[fl.intersections[1]])/self.L0
                    except AttributeError:
                        separation_buffer = separation_buffer_default #set = default value given as optional argument if reference line has no width
                    separation_distance = ArcArray(fl.coords)[fl.intersections[1]] + separation_buffer #where line separates from mainline

                    if use_mainline_tau:
                        fl.optimal_tau = self.network_tau
                        fl.yield_type = self.network_yield_type #this is probably too powerful, but unclear how else to exploit Bingham_number functionality
                    else:
                        pass
                    
                    if out_dict['Termini'][-1]/self.L0 < ArcArray(fl.coords)[fl.intersections[1]] : ## Below is only applicable while branches share single terminus 
                        print 'Previous terminus {} < {}'.format(out_dict['Termini'][-1]/self.L0, ArcArray(fl.coords)[fl.intersections[1]])
                        dLdt_branch = dLdt_annum
                        branchmodel_full = fl.plastic_profile(startpoint=new_termpos/self.L0, hinit=new_termheight, endpoint=fl_amax, surf=fl.surface_function) #model profile only down to intersection
                        mainbranch_int_idx = (np.abs(branchmodel_full[0] - ArcArray(fl.coords)[fl.intersections[1]])).argmin()
                        mainbranch_tribheight = refdict[new_key][1][mainbranch_int_idx] #finding surface elevation across main trough based on modelling of main branch
                        cutoff_yieldheight = BalanceThick(fl.bed_function(separation_distance)/fl.H0, fl.Bingham_num(0,0)) + (fl.bed_function(separation_distance)/fl.H0)
                        ##Note for constant yield, Bingham_num(0,0) = Bingham_num everywhere.  Might introduce bug for Mohr-Coulomb case
                        
                        branch_termheight = max(mainbranch_tribheight, cutoff_yieldheight) #identify and handle if tributary has separated from main branch based on whether main surface elevation is below reasonable intersection with trib
                        print 'Main branch tributary height: {}, cutoff yield height {}'.format(mainbranch_tribheight, cutoff_yieldheight)
                        if branch_termheight==cutoff_yieldheight: #tributary has fully separated from downstream branch
                            print 'Tributary {} has separated due to thinning disconnection from main branch at time {} a'.format(j, key+dt)
                            branch_terminus = separation_distance * self.L0
                            branchmodel = fl.plastic_profile(startpoint=separation_distance, hinit=cutoff_yieldheight, endpoint=fl_amax, surf=fl.surface_function)
                        ##print 'cutoff_idx = {}. Distance ={}'.format(cutoff_idx, branchmodel_full[0][cutoff_idx])
                        #cutoff_height = branchmodel_full[1][cutoff_idx]
                        #print 'Bed around cutoff point = {}, {}, {}'.format(branchmodel_full[2][cutoff_idx-1], branchmodel_full[2][cutoff_idx], branchmodel_full[2][cutoff_idx+1])
                        #print 'Cutoff_height = {}.  Cutoff_yieldheight = {}.  Branch_termheight = {}'.format(cutoff_height, cutoff_yieldheight, branch_termheight)
                        else:
                            branch_terminus = new_termpos
                            cutoff_idx = (np.abs(branchmodel_full[0] - separation_distance)).argmin() #retrieve index of closest point in full profile to the tributary cutoff
                            branchmodel = tuple([branchmodel_full[j][cutoff_idx::] for j in range(len(branchmodel_full))]) #truncate branchmodel profile to be counted only above cutoff point

                    else: ##if branches have split, find new terminus quantities
                        print 'Tributary {} has split from main branch at time {} a'.format(j, key+dt)
                        if k<1:
                            out_dict[0] = fl.plastic_profile(startpoint=0, hinit=fl.surface_function(0), endpoint=fl_amax, surf=fl.surface_function) #initial condition for time evolution - needed to calculate calving flux at first timestep
                            dLdt_branch = fl.dLdt_dimensional(profile=out_dict[0], alpha_dot=alpha_dot_k, debug_mode=debug_mode, dL=dL, has_smb=has_smb, terminus_balance=terminus_balance, submarine_melt=submarine_melt, rate_factor=rate_factor)
                        else:
                            dLdt_branch = fl.dLdt_dimensional(profile=out_dict[key], alpha_dot=alpha_dot_k, debug_mode=debug_mode, dL=dL, has_smb=has_smb, terminus_balance=terminus_balance, submarine_melt=submarine_melt, rate_factor=rate_factor)
                        if np.isnan(dLdt_branch):
                            dLdt_branch = out_dict['Termrates'][-1] /dt
                        else:
                            pass
                        branch_terminus_raw = out_dict['Termini'][-1] + (dLdt_branch*dt)
                        branch_terminus_posdef = max(0, branch_terminus_raw) #catching the rare case when branches have separated but one branch suddenly readvances past initial terminus
                        if branch_terminus_posdef > (fl_amax * self.L0):
                            print 'Terminus retreated past upstream limit. Resetting terminus position to = upstream limit.'
                            branch_terminus = out_dict['Termini'][-1] #catching whether terminus has retreated beyond upstream limit, setting to previous good value if so
                        else:
                            branch_terminus = branch_terminus_posdef
                        branch_term_bed = fl.bed_function(branch_terminus/self.L0)
                        previous_branch_bed = fl.bed_function(out_dict['Termini'][-1]/self.L0)
                        previous_branch_thickness = (out_dict['Terminus_heights'][-1] - previous_branch_bed)/self.H0
                        branch_termheight = BalanceThick(branch_term_bed/self.H0, fl.Bingham_num(previous_branch_bed/self.H0, previous_branch_thickness)) + (branch_term_bed/self.H0)    
                        branchmodel = fl.plastic_profile(startpoint=branch_terminus/self.L0, hinit=branch_termheight, endpoint=fl_amax, surf=fl.surface_function)
                    
                    if yr>dt:
                        if sl_contribution:
                            branch_termflux = fl.icediff_SL(profile1=out_dict[key], profile2=branchmodel, separation_buffer=separation_buffer)
                        else:
                            branch_termflux = fl.icediff(profile1=out_dict[key], profile2=branchmodel, separation_buffer=separation_buffer)
                    else:
                        branch_termflux = np.nan
                    
                    #print 'Year {} tributary {} terminus = {}'.format(yr, fl.index, branch_terminus)
                    out_dict[new_key] = branchmodel
                    out_dict['Termini'].append(branch_terminus)
                    out_dict['Terminus_heights'].append(branch_termheight*self.H0)
                    out_dict['Termrates'].append(dLdt_branch*dt)
                    if yr > dt:
                        out_dict['Terminus_flux'].append(branch_termflux)
                    else:
                        out_dict['Terminus_flux'].append(np.nan)
        
        light_output = [{'testyears': model_output_dicts[j]['testyears'], 'Termini': model_output_dicts[j]['Termini'], 'Terminus_flux': model_output_dicts[j]['Terminus_flux']} for j in range(len(self.flowlines))]
        
        if output_heavy:
            self.model_output = model_output_dicts
        else:
            self.model_output = light_output
        
        #self.model_output = model_output_dicts
  
        
    def save_network(self, filename=None):
        """Write essential information about a PlasticNetwork instance to a pickle.
        filename: defaults to self.name
        
        Output:
            'network_name': name as initialized
            'N_Flowlines': how many flowlines the network has
            'network_tau': optimal value of tau_y
            'network_yield_type': constant or variable yield strength
            'mainline_model_output': full profiles for each year, plus Termini, Termrates, Terminus_heights, and Terminus_flux
        """
        if filename is None:
            fn = str(self.name)
            fn1 = fn.replace(" ", "")
            fn2 = fn1.replace("[", "-")
            fn3 = fn2.replace("/", "_")
            fn4 = fn3.replace("]", "")
            filename='{}.pickle'.format(fn4)
        
        N_Flowlines = len(self.flowlines)
        
        try:
            balance_forcing = self.balance_forcing
        except AttributeError: #likely will not have attribute if only saving runs forced from upstream
            balance_forcing = None
        
        output_dict = {
        'network_name': self.name,
        'N_Flowlines': N_Flowlines,
        'network_tau': self.network_tau,
        'network_yield_type': self.network_yield_type,
        'balance_forcing': balance_forcing,
        'mainline_model_output': self.model_output[0]
        }
        
        if N_Flowlines >1:
            for n in range(N_Flowlines):
                output_n = self.model_output[n]
                key_n = 'model_output_'+str(n)
                output_dict[key_n] = output_n
            else:
                pass
        
        with open(filename, 'wb') as handle:
            pickle.dump(output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def load_network(self, filename, load_mainline_output=False, load_tributary_output=False):
        """Loads in pickled information about a network previously run.  You will still have to run process_full_lines and network_ref_profiles if you want to do new model runs.
        load_model_output: default True loads in saved model output from saved run.  Useful for analysis without re-running.
        """
        with open(filename, 'rb') as handle:
            loaded_dict = pickle.load(handle)
        print loaded_dict['network_name']
        
        if loaded_dict['network_name'] != self.name:
            warnings.warn('Name of loaded network does not match name of this network instance.  Check that you have loaded the correct pickle.')
        
        N_Flowlines = loaded_dict['N_Flowlines']
        
        self.network_tau = loaded_dict['network_tau']
        self.network_yield_type = loaded_dict['network_yield_type']
        try:
            self.balance_forcing = loaded_dict['balance_forcing']
        except KeyError: #if balance forcing was not saved
            self.balance_forcing = None
        if load_mainline_output:
            self.model_output = {} #needs to be initialized
            self.model_output[0] = loaded_dict['mainline_model_output']
        if N_Flowlines >1 and load_tributary_output:
            for n in range(N_Flowlines)[1::]:
                key_n = 'model_output_'+str(n)
                load_in_output_n = loaded_dict[key_n]
                self.model_output[n] = load_in_output_n
        else:
            pass
            