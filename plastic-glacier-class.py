# Defining a class to run plastic model functions
# 16 Nov 17  EHU
from plastic_utilities_v2 import *

##Bingham stress functions-----------------
#def B_var(elev, thick, pos, time): #variation by altitude and ice thickness (effective pressure at base)...pos, time arguments required by plasticmodel
#    if elev<0:
#        D = -elev #Water depth D the nondim bed topography value when Z<0
#    else:
#        D = 0
#    N = rho_ice*g*thick*H0 - rho_sea*g*D*H0
#    mu = 0.01 #Coefficient between 0 and 1
#    tau_y = tau_0 + mu*N
#    return tau_y/(rho_ice*g*H0**2/L0) 
#    
#def B_const(elev, thick, pos, time):  #functional form of B if using a constant yield strength
#    return tau_yield/(rho_ice*g*H0**2/L0)

class PlasticGlacier:
    def __init__(self, glaciername, flowline, bedfunction, yield_type, tau_val, obs_SE_function, inityear, years=None, linenumber=None):
        self.name = glaciername
        self.xy = flowline
        self.bed = bedfunction
        self.tau = tau_val
        self.obs_SE = obs_SE_function
        self.inityear = inityear
        #self.year_range = years
        #if years is None:
        #    print 'Warning: you will need to enter a valid year range to model time evolution'
        self.yield_type = yield_type
        if linenumber is not None:
            self.linenumber = linenumber #to allow for networks with more than one line
        #if yield_type.lower() in ['v', 'var', 'variable', 'coulomb']:
        #    self.yieldfunc = B_var
        #    print '{} Glacier line successfully initialised with variable yield strength.'.format(glaciername)
        #elif yield_type.lower() in ['c', 'cons', 'const', 'constant']:
        #    self.yieldfunc = B_const
        #    print '{} Glacier line successfully initialised with constant yield strength.'.format(glaciername)
        #else:
        #    print 'Unrecognised yield type. Please select constant or variable (Coulomb) yield strength.'
    
    def yieldfunc(self, elev, thick, pos=None, time=None): #rewriting functions B_var and B_const to root out bug
        if self.yield_type.lower() in ['v', 'var', 'variable', 'coulomb']:
            if elev<0:
                D = -elev #Water depth D the nondim bed topography value when Z<0
            else:
                D = 0
            N = rho_ice*g*thick*H0 - rho_sea*g*D*H0
            mu = 0.01 #Coefficient between 0 and 1
            tau_y = self.tau + mu*N
            return tau_y/(rho_ice*g*H0**2/L0) 
        elif self.yield_type.lower() in ['c', 'cons', 'const', 'constant']:
            return self.tau/(rho_ice*g*H0**2/L0)
    
    def set_year_range(self, years):
        """Enter years when snapshots should be generated for time evolution."""
        self.year_range = years
            
    def length(self):
        a = ArcArray(self.xy)
        return (L0/1000)*a[-1]
        
    def Profile(self, startpoint, endpoint, Npoints, term_initialisation=None, h_in=None, is_initial=None):
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
        
        if term_initialisation.lower() in ['o', 'obs', 'observed']:
            hinit = self.obs_SE(startpoint)
        elif term_initialisation.lower() in ['b', 'balance', 'balancethick', 'model', 'modelled']:
            hinit = BalanceThick((self.bed(startpoint)/H0), self.yieldfunc(self.bed(startpoint)/H0, (self.obs_SE(startpoint)-self.bed(startpoint))/H0, 0, 0)) + (self.bed(startpoint)/H0)
        elif term_initialisation.lower() in ['midstream', 'no', 'false', 'not term', 'not terminus']:
            hinit = h_in
            
        
        SEarr.append(hinit)
        thickarr.append(hinit-(self.bed(startpoint)/H0))
        basearr.append(self.bed(startpoint)/H0)
        obsarr.append((self.obs_SE(startpoint))/H0)
        
        for x in horiz[1::]:
            bedval = self.bed(x)/H0  # value of interpolated bed function
            obsheight = (self.obs_SE(x))/H0
            modelthick = thickarr[-1]
            B = self.yieldfunc(bedval, modelthick, None, None)
            #Break statements for thinning below yield, water balance, or flotation
            if dx<0:
                if modelthick<BalanceThick(bedval,B):
                    print 'Thinned below water balance at x=' + str(10*x)+'km'
                    break
            if modelthick<FlotationThick(bedval):
                print 'Thinned below flotation at x=' + str(10*x)+'km'
                break
            if modelthick<4*B*H0/L0:
                print 'Thinned below yield at x=' +str(10*x)+'km'
                break
            else:
                basearr.append(bedval)
                SEarr.append(SEarr[-1]+(B/modelthick)*dx) 
                thickarr.append(SEarr[-1]-basearr[-1])
                obsarr.append(obsheight)
        
        error = np.sqrt(((np.array(SEarr)-np.array(obsarr))**2).mean())
        CVrms = error/mean(SEarr)
        #print 'RMS error: '+ str(error) +',   CV(RMSE): ' + str(CVrms)
        
        if is_initial is True:
            self.init_profile = (horiz[0:len(SEarr)], SEarr, basearr)
        
        return (horiz[0:len(SEarr)], SEarr, basearr, obsarr)
    
    def TimeEvol_simple(self, refpoint, dHfunc, simyears, ups_max=None):
        """Refpoint is the upstream point where thinning/thickening will be applied - express in m/L0
        dHfunc should be a function that takes a year as input and returns the change in ice thickness at the reference point"""
        
        ## integrate initial profile to make a function
        modelint = interpolate.interp1d(self.init_profile[0], self.init_profile[1], kind='linear', copy=True)
        refheight = modelint(refpoint)
        if ups_max is None:
            arcmax = (1000/L0)*self.length
        else:
            arcmax = ups_max
        
        dHvals = dHfunc(simyears)
        modeltermini = [0]
        termrates = []
        self.TimeEvol = {}
        
        for j, year in enumerate(simyears):
            dH = dHvals[j]
            fwdmodel = Profile(refpoint, 0, 25000, h_in=refheight+dH, is_initial=False)
            bkmodel = Profile(refpoint, arcmax, 25000, h_in=refheight+dH, is_initial=False)
            modelterm = L0*min(fwdmodel[0]) #in m
            
            model_x = np.concatenate((fwdmodel[0][::-1], bkmodel[0]))
            model_prof = np.concatenate((fwdmodel[1][::-1], bkmodel[1]))
            self.snapshot[year] = (model_x, model_prof) #can I save arrays of arrays to dictionaries?
            
            dL = modelterm - modeltermini[-1]
            modeltermini.append(modelterm)
            termrates.append(dL) #dt = 1 by definition
    
    #def TimeEvol_SMB(self, SMBfunction):
    ## change length
    ## make new profile

    