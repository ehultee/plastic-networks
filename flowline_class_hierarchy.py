## Meta-class of plastic glacier, which should have flowlines as classes of their own (because this is what we run model on)
## 29 Jan 2017  EHU

import warnings
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
    def __init__(self, H0=1e3, L0=10e3, g=9.8, rho_ice=920.0, rho_sea=1020.0, default_Ty=150e3, default_T0=130e3):
        self.H0 = H0 #characteristic height for nondimensionalisation 
        self.L0 = L0
        self.g = g
        self.rho_ice = 920.0
        self.rho_sea = 1020.0
        self.default_Ty = 150e3
        self.default_T0 = 130e3

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
    """
    def __init__(self, coords, index=None, order=0, flows_to=None, intersect=None):
        Ice.__init__(self)
        self.coords = np.asarray(coords) 
        self.order = order
        self.flows_to = flows_to
        self.intersect = intersect 
    
    #def make_full_line(self):
         #TRANSFERRED THIS FUNCTIONALITY TO PLASTICNETWORK BECAUSE THAT'S THE CONTEXT WHERE IT WILL BE USED
 
                     

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
    
    def __init__(self, coords, index=None, name=None, initial_terminus=(0,0), intersections=None):
        Ice.__init__(self)
        self.coords = np.asarray(coords)
        self.length = ArcArray(self.coords)[-1]
        if index is None:
            self.index = np.nan
        else:
            self.index = index
        if name is None:
            self.name = 'Glacier'
        else:
            self.name = name
        #if flows_to is not None:
        #    self.flows_to = flows_to
            
    
    def process_bed(self, B_field):
        """Make callable bed function for this flowline.  B_field should be 2d-interpolated."""
        self.bed_function = FlowProcess(self.coords, B_field)
    
    def process_surface(self, S_field):
        """Make callable surface function for this flowline.  S_field should be 2d-interpolated."""
        self.surface_function = FlowProcess(self.coords, S_field)
    
    def process_thickness(self, H_field):
        """Make callable ice thickness function along this flowline.  H_field should be 2d-interpolated."""
        self.thickness_function = FlowProcess(self.coords, H_field)
    
    def optimize_yield_strength(self, testrange=range(50e3, 150e3, 5e3), arcmax=self.length):
        """Run optimization and set the result to be the optimal value for the flowline instance.  
        Arcmax over which to run optimization can be adjusted according to how high up we think plastic approximation should go."""
        CV_const_arr = []
        CV_var_arr = []
        
        bedf = self.bed_function
        surf = self.bed_function
        
        for tau in testrange:
                
            ##CONSTANT YIELD
            model_const = plasticmodel_error(bedf, tau, B_const, 0, surf(0)/H0, arcmax, 25000, surf) #prescribed terminus thickness
            CV_const = model_const[1]
            ##VARIABLE YIELD
            model_var = plasticmodel_error(bedf, tau, B_var, 0, surf(0)/H0, arcmax, 25000, surf) #prescribed terminus thickness
            CV_var = model_var[1]
            
            CV_const_arr.append(CV_const)
            CV_var_arr.append(CV_var)
        
        constopt_index = np.argmin(CV_const_arr)
        varopt_index = np.argmin(CV_var_arr)
        
        constopt = testrange[constopt_index]
        varopt = testrange[varopt_index]
        
        if np.min(CV_const_arr)<np.min(CV_var_arr):
            yieldtype = 'constant'
            t_opt = constopt
        elif np.min(CV_const_arr_)>np.min(CV_var_arr):
            yieldtype = 'variable'
            t_opt = varopt
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
        if self.yield_type is 'constant':
            return self.optimal_tau/(rho_ice*g*H0**2/L0) #B_const from plastic_utilities
            
        elif self.yield_type is 'variable':
            #B_var from plastic_utilities, accounts for water pressure
            if elev<0:
                D = -elev #Water depth D the nondim bed topography value when Z<0
            else: 
                D = 0
            N = rho_ice*g*thick*H0 - rho_sea*g*D*H0
            mu = 0.01 #coefficient between 0 and 1, scales influence of water pressure
            tau_y = self.optimal_tau/ + mu*N
            return tau_y/(rho_ice*g*H0**2/L0)
            
    def set_ref_profile(self):
        """Creates 1d interpolated function of reference profile height along flowline.  
        Give argument of nondim arclength, returns nondim surface elevation.  Useful for initial condition in thinning."""
        plasticref = PlasticProfile(self.bed_function, self.Bingham_num, 0, self.surface_function(0)/H0, self.length, 10000, self.surface_function)
        ref_prof_interpolated = interpolate.interp1d(plasticref[0], plasticref[1], kind='linear', copy=True)
        self.ref_profile = ref_prof_interpolated
    
  
              
class PlasticNetwork(Ice):
    """A container for a named glacier on which we want to run the plastic model.  Can take instances of Branch as arguments and make them Flowlines.
    
    Attributes:
        name: A string with what we call the glacier for modelling/analysis
        init_type: Can be 'Branch' or 'Flowline'. If we already have Flowlines we don't need to call make_full_lines
        branches: A list of Flowline objects describing the branches of this glacier
        main_terminus: Coordinates of the terminus of the "main branch" (index 0) in the initial state
    """
    
    def __init__(self, name='Glacier', init_type='Branch', branches=[], main_terminus=(0,0)):
        Ice.__init__(self)
        self.name = name
        if init_type is 'Branch':
            self.branches = branches
            self.flowlines = []
            print 'You have initialised with glacier Branches.  Please call make_full_lines and process_full_lines to proceed.'
        elif init_type is 'Flowline':
            self.flowlines = branches
            print 'You have initialised with glacier Flowlines.  Ready to proceed with simulations.'
        else:
            self.branches = branches
            self.flowlines = []
            print 'Unrecognised initialisation type.'
        self.main_terminus = main_terminus
    
    def make_full_lines(self):
        #KDTree search just connects at closest points...may still be wonky
        #project somehow with Shapely?
        mainline = self.branches[0].coords
        maintree = spatial.KDTree(mainline)
        
        full_lines = {}
        j = 0
        

        while j<len(self.branches): #making full length flowline for each branch
            branch = self.branches[j] 
            pt = branch[0]
            dist, idx = maintree.query(pt, distance_upper_bound=5000)
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
                fullbranch = np.concatenate((mainfrag, branch)) #Shapely project here
                full_lines[j] = fullbranch
                j+=1
        
        while k<len(full_lines):
            coordinates = np.asarray(full_lines[k])
            self.flowlines[k] = Flowline(coordinates, index=k)
    
    def process_full_lines(self, bed_field, surface_field, thickness_field):
        """Processes bed, surface, and thickness functions for new flowlines.  Use if network initialized with Branches.
        """
        for line in self.flowlines:
            line.process_bed(bed_field)
            line.process_surface(surface_field)
            line.process_thickness(thickness_field)

    def optimize_network_yield(self, testrange=range(50e3, 150e3, 5e3), arcmax_list = [self.flowlines[i].length for i in range(len(self.flowlines))], check_all=False):
        """Runs yield strength optimization on each of the network Flowline objects and sets the yield strength of the 'main' branch as the network's default.
        Could find a better way to do this...
        arcmax_list: list of how high up to optimise in order of flowline index; default is full length
        check_all=False/True decides whether to do full optimisation to find tau_y on each line (True) or just take the main line value (False)
        """
        
        if check_all:
            for j,line in enumerate(self.flowlines):
                line.optimize_yield_strength(testrange, arcmax=arcmax_list[j]) #run optimisation for each flowline
        else:
            try:
                self.network_yield_type = self.flowlines[0].yield_type
                self.network_tau = self.flowlines[0].optimal_tau
            except AttributeError:
                self.flowlines[0].optimize_yield_strength(testrange, arcmax=arcmax_list[0])
                self.network_yield_type = self.flowlines[0].yield_type
                self.network_tau = self.flowlines[0].optimal_tau

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
    
    def network_time_evolve(self, testyears=arange(100), ref_branch_index=0, upgl_ref=15000/L0, thinrate=10/H0, thinvalues=None, upstream_limits=[fl.length for fl in self.flowlines], use_mainline_tau=True):
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
    
        returns [TYPE] model output
        """
        if thinvalues is None:  
            thinvals = np.full(len(testyears), thinrate)
        else:
            thinvals = thinvalues
        
        dt = mean(diff(testyears))
        
        model_output_dicts = [{} for l in self.flowlines]
        for j,fl in enumerate(self.flowlines):
            #print 'Currently running {} line {}'.format(self.name, j)
            #sarr = fl.ref_profile[0] #arclength array of reference profile...but here this is an interpolated function, so might not work
            amax = upstream_limits[j]
            #amin = ArcArray(fl.coords) #probably always 0?
            refpt = min(amax, upgl_ref) #forces at top of branch if shorter than reference distance.  In general would expect forcing farther downstream to give weaker response
            refht = fl.ref_profile(refpt)
            
            bed = fl.bed_function
            surf = fl.surface_fuction
            if use_mainline_tau:
                tau_j = self.network_tau
                fl.yield_type = self.network_yield_type #this is probably too powerful, but unclear how else to exploit Bingham_number functionality
            else:
                tau_j = fl.optimal_tau #use mainline value and type, surely?
            
            dmodel = model_output_dicts[j]
            #dmodel['Termini'] = [(self.L0)*amin]
            dmodel['Termini'] = [0]
            dmodel['Termrates'] = []
            
            #NEED TO DO THIS DIFFERENTLY FOR MAINLINE THAN FOR OTHER LINES
            #TIME STEP OUTSIDE OF FLOWLINE FOR-LOOP
            for k, yr in enumerate(testyears):
                thinning = np.sum(thinvals[:k])
                fwdmodel = PlasticProfile(bed, tau_j, fl.Bingham_num, refpt, refht - thinning, 0, 25000, surf)
                bkmodel = PlasticProfile(bed, tau_j, fl.Bingham_num, refpt, refht - thinning, amax, 25000, surf)
                modelterm = L0*min(fwdmodel[0]) #in m
                dL = modelterm - dmodel['Termini'][-1]
                dmodel[yr] = fwdmodel
                dmodel['Termini'].append(modelterm)
                dmodel['Termrates'].append(dL/dt)
                
            