## Debugging Helheim hanging tributary
## 14 Jun 2018  EHU

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

#Smoothing bed to check effect on dLdt
unsmoothB = B
smoothB = gaussian_filter(B, 2)
#B_processed = np.ma.masked_where(thick_mask !=2, smoothB)

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], S.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])

## Initializing Helheim network
helcoords_0, helcoords_1, helcoords_2 = Flowline_CSV('Helheim-network-w_width.csv', 3, has_width=True, flip_order=False)
hel_0 = Branch(coords=helcoords_0, index=0, order=0)
hel_1 = Branch(coords=helcoords_1, index=1, order=1, flows_to=0)
hel_2 = Branch(coords=helcoords_2, index=2, order=1, flows_to=0)
hel_branches = (hel_0, hel_1, hel_2)
Helheim = PlasticNetwork(name='Helheim', init_type='Branch', branches=hel_branches, main_terminus=helcoords_0[0])
Helheim.make_full_lines()
Helheim.load_network(filename='Helheim.pickle')


#glacier_networks=(Helheim,)

Helheim.process_full_lines(B_interp, S_interp, H_interp)
#Helheim.optimize_network_yield(check_all=True, use_balancethick=True)
for fln in Helheim.flowlines:
    fln.yield_type  = Helheim.network_yield_type
    fln.optimal_tau = Helheim.network_tau
Helheim.network_ref_profiles()

###-------------------
#### FORWARD PROJECTION--FORCING FROM TERMINUS
###-------------------
testyears = arange(100, step=0.25)
db=False

Helheim.terminus_time_evolve(testyears=testyears, alpha_dot=0, dL=1/L0, separation_buffer=10000/L0, has_smb=False, terminus_balance=0, submarine_melt = 0, debug_mode=db, rate_factor=3.7E-26) 
Helheim_multibranch_flux = [Helheim.model_output[j]['Terminus_flux'] for j in range(len(Helheim.flowlines))]


#####-------------------
###### PROJECTION--HANGING TRIBUTARY ONLY
#####-------------------
#fl = Helheim.flowlines[1]
#fl_amax = fl.length
#separation_buffer = 10000/L0
#separation_distance = ArcArray(fl.coords)[fl.intersections[1]] + separation_buffer

#testyears = arange(20, step=0.5)
#dt = 0.5
#dt_rounding=3
#trib_dict = {'Termini': [separation_distance],
#        'Terminus_heights': [fl.surface_function(separation_distance)],
#        'Termrates': [],
#        'Terminus_flux': [],
#        'testyears': testyears
#        }
#trib_dict[0] = fl.plastic_profile(startpoint=separation_distance, hinit=fl.surface_function(separation_distance), endpoint=fl_amax, surf=fl.ref_profile)
#
#for k, yr in enumerate(testyears):
#    print yr
#    key = round(yr-dt, dt_rounding) #allows dictionary referencing when dt is < 1 a
#            
#    if k<1:
#        dLdt_annum = fl.dLdt_dimensional(profile=trib_dict[0], alpha_dot=0, debug_mode=False, dL=1/L0, has_smb=False, terminus_balance=0, submarine_melt=0, rate_factor=3.7E-26)
#    else:
#        dLdt_annum = fl.dLdt_dimensional(profile=trib_dict[key], alpha_dot=0, debug_mode=False, dL=1/L0, has_smb=False, terminus_balance=0, submarine_melt=0, rate_factor=3.7E-26)
#    if np.isnan(dLdt_annum): #happens if retreat hits edge of domain unexpectedly
#        dLdt_annum = trib_dict['Termrates'][-1] / dt #replace with last non-nan value
#    else:
#        pass
#    
#    new_termpos_raw = trib_dict['Termini'][-1]+(dLdt_annum*dt) #Multiply by dt in case dt!=1 annum.  Multiply dLdt by L0 because its output is nondimensional
#    new_termpos_posdef = max(0, new_termpos_raw)
#    if new_termpos_posdef > (fl_amax * fl.L0):
#        print 'Terminus retreated past upstream limit. Resetting terminus position to = upstream limit.'
#        new_termpos = trib_dict['Termini'][-1] #set to previous good position
#        #new_termpos = ref_amax *self.L0 #glacier sits at edge of domain if terminus retreats beyond upstream limit
#    else:
#        new_termpos = new_termpos_posdef
#    previous_bed = fl.bed_function(trib_dict['Termini'][-1]/fl.L0)
#    try:
#        new_term_bed = fl.bed_function(new_termpos/fl.L0)
#    except ValueError: #interpolated bed function sometimes complains that ref_amax * L0 is above the interpolation range.  Use last good value if this is the case.
#        new_term_bed = previous_bed
#    previous_thickness = (trib_dict['Terminus_heights'][-1] - previous_bed)/fl.H0 #nondimensional thickness for use in Bingham number
#    new_termheight = BalanceThick(new_term_bed/fl.H0, fl.Bingham_num(previous_bed/fl.H0, previous_thickness)) + (new_term_bed/fl.H0)
#    new_profile = fl.plastic_profile(startpoint=new_termpos/fl.L0, hinit=new_termheight, endpoint=fl_amax, surf=fl.ref_profile)
#    if yr>dt:
#        #key = round(yr-dt, dt_rounding)
#        termflux = fl.icediff(profile1=trib_dict[key], profile2=new_profile)
#    else:
#        termflux = np.nan
#    
#    new_key = round(yr, dt_rounding)    
#    trib_dict[new_key] = new_profile
#    trib_dict['Terminus_flux'].append(termflux)
#    trib_dict['Termini'].append(new_termpos)
#    trib_dict['Terminus_heights'].append(new_termheight*fl.H0)
#    trib_dict['Termrates'].append(dLdt_annum*dt)



####-------------------
##### SUMMARY PLOTTING
####-------------------

##projections = [Jakobshavn_main.model_output, Jakobshavn_sec.model_output, Jakobshavn_tert.model_output, KogeBugt.model_output, Helheim.model_output, Kanger.model_output]
names = ['Helheim [main]', 'Helheim [line 1]', 'Helheim [line 2]']
#rates = ['{0:.2f} m/a'.format(gl.sec_alphadot) for gl in glacier_networks]
##styles = [':', '-.', '--', '-']
markers = ['o', '.', ',', '^', 'd', '*']
styles = ['-', ':', '-.', '-', '-', '-']
cmap = matplotlib.cm.get_cmap('winter')
colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
alt_colors = matplotlib.cm.get_cmap('Greys')([0.3, 0.5, 0.7, 0.9])


plt.figure(figsize=(12,8))
for j in range(len(Helheim_multibranch_flux)):
    plt.plot(testyears, 1E-12*np.array(Helheim_multibranch_flux[j]), linewidth=4, linestyle=styles[j], color=colors[j], label=names[j])
    plt.plot(testyears[::50], 1E-12*np.array(Helheim_multibranch_flux[j][::50]), linewidth=0, marker=markers[j], ms=10, color=colors[j])
    plt.fill_between(testyears, y1=1E-12*np.array(Helheim_multibranch_flux[j]), y2=0, color=colors[j], alpha=0.5)    
#plt.legend(loc='upper right')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Terminus ice flux [Gt/a]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0,15.5)
#plt.axes().set_xticks([0,5,10, 15])
#plt.axes().set_ylim(0, 18)
#plt.axes().set_yticks([0, 5, 10, 15])
plt.show()

plt.figure('Helheim terminus change')
for k, mo in enumerate(Helheim.model_output): #for each branch j
    #colork = matplotlib.cm.get_cmap('viridis')(k/len(Helheim.model_output))
    #markerk = (k+2, mod(k+1, 3), 0)
    lsk = (':', '-.', '--')
    plt.plot(testyears[::], -0.001*np.array(mo['Termini'])[:-1:], linewidth=4, color='k', ls=lsk[k], label='Helheim line {}'.format(k))
    #plt.plot(testyears[::10], -0.001*np.array(mo['Termini'])[:-1:10], linewidth=0, color='k', marker=markerk, ms=10, label='Helheim line {}'.format(k))
#plt.legend(loc='upper right')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Terminus change [km]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0, 20)
#plt.axes().set_xticks([0, 5, 10, 15, 20])
#plt.axes().set_xticklabels(['0', '', '', '', '20'])
#plt.axes().set_ylim(-16, 1)
#plt.axes().set_yticks([-15, -10, -5, 0])
plt.show()