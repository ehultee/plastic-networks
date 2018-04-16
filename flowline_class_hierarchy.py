## Meta-class of plastic glacier, which should have flowlines as classes of their own (because this is what we run model on)
## 29 Jan 2017  EHU

import warnings
import pickle
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
        if has_width:
            try:
                self.width = np.asarray(coords)[:,2]
            except IndexError:
                self.width = width_array
        else:
            self.width = None
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
    
    def process_width(self):
        self.width_function = interpolate.interp1d(ArcArray(self.coords), self.width)
    
    def optimize_yield_strength(self, testrange=np.arange(50e3, 500e3, 5e3), arcmax=None):
        """Run optimization and set the result to be the optimal value for the flowline instance.  
        Arcmax over which to run optimization can be adjusted according to how high up we think plastic approximation should go.
        NOTE: arcmax default is None but will set to self.length.  self unavailable at function define-time"""
        
        #Fixing default that could not be used in definition of arguments
        if arcmax is None:
            arcmax = self.length
        CV_const_arr = []
        CV_var_arr = []
        
        bedf = self.bed_function
        surf = self.surface_function
        
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
        elif np.min(CV_const_arr)>np.min(CV_var_arr):
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
            tau_y = self.optimal_tau/ + mu*N
            return tau_y/(rho_ice*g*H0**2/L0)
        
        else: #do this for constant and as default
            return self.optimal_tau/(rho_ice*g*H0**2/L0) #B_const from plastic_utilities
            
    def plastic_profile(self, bedf=None, startpoint=0, hinit=None, endpoint=None, Npoints=25000, surf=None):
        """
        This snippet translates functionality of plastic_utilities_v2.PlasticProfile"""
        
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
            #Break statements for thinning below yield, water balance, or flotation
            if dx<0:
                if modelthick<BalanceThick(bed,B):
                    print 'Thinned below water balance at x=' + str(10*x)+'km'
                    break
            if modelthick<FlotationThick(bed) and dx<0: #should only break for flotation when running downstream in simulations 
                ##CONFIRM THAT THERE ARE OTHER WAYS TO CATCH FLOATING TERMINI 
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
        
        return (horiz[0:len(SEarr)], SEarr, basearr, obsarr)        

    def set_ref_profile(self):
        """Creates 1d interpolated function of reference profile height along flowline.
        Assigns function to be self.ref_profile  
        Give argument of nondim arclength, returns nondim surface elevation.  Useful for initial condition in thinning."""
        #plasticref = PlasticProfile(self.bed_function, self.Bingham_num, 0, self.surface_function(0)/H0, self.length, 10000, self.surface_function)
        plasticref = self.plastic_profile()
        ref_prof_interpolated = interpolate.interp1d(plasticref[0], plasticref[1], kind='linear', copy=True)
        self.ref_profile = ref_prof_interpolated
    
  
              
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
        
        """
        #KDTree search just connects at closest points...may still be wonky
        #project somehow with Shapely?
        mainline = self.branches[0].full_coords
        maintree = spatial.KDTree(mainline[:,0:2])
        
        full_lines = {}
        j = 0
        

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
                if idx==len(full_lines[j-1]):
                    print 'Branch {} also does not intersect tributary {}.  Appending raw line.  Use with caution.'.format(j, j-1)
                    full_lines[j] = br
                else:
                    tribfrag = branchlist[j-1][:idx_t]
                    fullbranch = np.concatenate((tribfrag, br))
                    full_lines[j] = fullbranch
                j+=1
            else:
                print mainline[idx]
                mainfrag = mainline[:idx]
                fullbranch = np.concatenate((mainfrag, br)) #Shapely project here
                full_lines[j] = fullbranch
                j+=1
        
        k = 0
        while k<len(full_lines):
            coordinates = np.asarray(full_lines[k])
            self.flowlines.append(Flowline(coordinates, index=k))
            k+=1
    
    def process_full_lines(self, bed_field, surface_field, thickness_field):
        """Processes bed, surface, and thickness functions for new flowlines.  Use if network initialized with Branches.
        """
        for line in self.flowlines:
            line.process_bed(bed_field)
            line.process_surface(surface_field)
            line.process_thickness(thickness_field)
            line.process_width()
            
            if line.thickness_function(0) < FlotationThick(line.bed_function(0)): #uses FlotationThick from plastic_utilities_v2
                warnings.warn('{} line {} may have a floating terminus.  Run remove_floating and re-initialise.'.format(self.name, line.index))
    
    def remove_floating(self):
        """Checks mainline for floating terminus, removes offending coordinates.
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
        cleanfrom = nonfloating[0] #first point where thickness exceeds flotation
        
        cleaned_coords = mainline.coords[cleanfrom::]
        trimmed_width = mainline.width[cleanfrom::]
        
        self.flowlines[0].coords = cleaned_coords
        self.flowlines[0].width = trimmed_width
        self.branches[0].coords = cleaned_coords
        self.branches[0].width = trimmed_width
        self.flowlines[0].length = ArcArray(cleaned_coords)[-1]
        
        print 'Coordinates with ice thickness less than flotation removed.  Please re-run make_full_lines and process_full_lines for network {}.'.format(self.name)
            
    
    def optimize_network_yield(self, testrange=arange(50e3, 500e3, 5e3), arcmax_list = None, check_all=False):
        """Runs yield strength optimization on each of the network Flowline objects and sets the yield strength of the 'main' branch as the network's default.
        Could find a better way to do this...
        arcmax_list: list of how high up to optimise in order of flowline index; default is full length (NOTE: default of "None" is set to more sensible default in function.  Same problem with self unavailable at define-time)
        check_all=False/True decides whether to do full optimisation to find tau_y on each line (True) or just take the main line value (False)
        """
        #Fixing default that could not be used in definition of arguments
        if arcmax_list is None:
            arcmax_list = [self.flowlines[i].length for i in range(len(self.flowlines))]
        
        try:
            self.network_yield_type = self.flowlines[0].yield_type
            self.network_tau = self.flowlines[0].optimal_tau
        except AttributeError:
            self.flowlines[0].optimize_yield_strength(testrange, arcmax=arcmax_list[0])
            self.network_yield_type = self.flowlines[0].yield_type
            self.network_tau = self.flowlines[0].optimal_tau
        
        if check_all:
            for j,line in enumerate(self.flowlines):
                line.optimize_yield_strength(testrange, arcmax=arcmax_list[j]) #run optimisation for each flowline
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
 
    def icediff(self, profile1, profile2, interpolated_func1=None, interpolated_func2=None, shape='frustum'):
        """Calculate net ice loss due to calving between two plastic profiles (interpret as successive timesteps)
        Inputs:
            profile1: an output from Flowline.plastic_profile
            profile2: an output from Flowline.plastic_profile (at a later time step)
            interpolated_func1: an interpolated surface elevation function describing profile1 in terms of arc length
            interpolated_func2: an interpolated surface elevation function describing profile2 in terms of arc length
            shape: 'frustum' ##Add others later e.g. parabolic bed
        Output:
            dM, ice mass change at the terminus (due to calving flux) between profile 1 and profile 2
        """
        if interpolated_func1 is None:
            interpolated_func1 = interpolate.interp1d(profile1[0], profile1[1], kind='linear', copy=True) #Creating interpolated surface elevation profile
        if interpolated_func2 is None:
            interpolated_func2 = interpolate.interp1d(profile2[0], profile2[1], kind='linear', copy=True) #Creating interpolated profile if needed
        
        ## Limits of integration
        x1 = min(profile1[0]) #initial terminus position, in nondimensional units
        x2 = min(profile2[0]) #new terminus position
        dX = (x2-x1)*self.L0 #in physical units of m
        print 'dX={} m'.format(dX)
        
        w1 = self.flowlines[0].width_function(x1)
        w2 = self.flowlines[0].width_function(x2) ## ADD BRANCH SEPARATION
        #dW = abs(w2-w1) #in physical units of m
        print 'dW={} m'.format(w2-w1)
        
        b1 = self.flowlines[0].bed_function(x1) ## ADD BRANCH SEPARATION
        b2 = self.flowlines[0].bed_function(x2) #physical units of m
        print 'Bed change = {}-{} m'.format(b2, b1)
        thickness1 = self.H0*(interpolated_func1(x1))-b1 
        thickness2 = self.H0*(interpolated_func1(x2))-b2 #determine if this comes from profile1 or profile2
        #dH = abs(thickness2-thickness1)
        print 'Thickness change = {} m'.format(thickness2-thickness1)
        
        dV_frustum = dX*((w2*thickness2) + (w2+w1)*(thickness2+thickness1) + (w1*thickness1))/6 #volume change, approximated as frustum of a pyramid
        print dV_frustum
        
        dM = dV_frustum * self.rho_ice #mass change
        
        return dM #will be negative in cases of ice loss
               

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
                termflux = self.icediff(profile1=refdict[yr-dt], profile2=fwdmodel)
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

        self.model_output = model_output_dicts
        
    
    def save_network_instance(self, filename=None):
        if filename is None:
            fn = str(self.name)
            fn1 = fn.replace(" ", "")
            fn2 = fn1.replace("[", "-")
            fn3 = fn2.replace("/", "_")
            fn4 = fn3.replace("]", "")
            filename='{}.pickle'.format(fn4)
        
        N_Flowlines = len(self.flowlines)
        
        output_dict = {
        'network_name': self.name,
        'N_Flowlines': N_Flowlines,
        'network_tau': self.network_tau,
        'network_yield_type': self.network_yield_type
        }
        
        with open(filename, 'wb') as handle:
            pickle.dump(output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def load_network_instance(self, filename):
        with open(filename, 'rb') as handle:
            loaded_dict = pickle.load(handle)
        
        print loaded_dict['network_name']
        
        self.network_tau = loaded_dict['network_tau']
        self.network_yield_type = loaded_dict['network_yield_type']
