# Forward projections with calving flux on Greenland glaciers
# 9 Apr 2018  EHU

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
### FINDING GLACIERS
##-------------------

jakcoords_main = Flowline_CSV('Documents/GitHub/plastic-networks/jakobshavn-mainline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_0 = Flowline(coords=jakcoords_main, index=0, name='Jak mainline', has_width=True)
Jakobshavn_main = PlasticNetwork(name='Jakobshavn Isbrae [main/south]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_main[0])
Jakobshavn_main.load_network(filename='JakobshavnIsbrae-main_south.pickle')

jakcoords_sec = Flowline_CSV('Jakobshavn_secondary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_1 = Flowline(coords=jakcoords_sec, index=1, name='Jak secondary branch', has_width=True)
Jakobshavn_sec = PlasticNetwork(name='Jakobshavn Isbrae [secondary/central]', init_type='Flowline', branches=(jak_1), main_terminus=jakcoords_sec[0])
Jakobshavn_sec.load_network(filename='Jakobshavn_sec.pickle')

jakcoords_tert = Flowline_CSV('Jakobshavn_tertiary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_2 = Flowline(coords=jakcoords_tert, index=2, name='Jak tertiary branch', has_width=True)
Jakobshavn_tert = PlasticNetwork(name='Jakobshavn Isbrae [tertiary/north]', init_type='Flowline', branches=(jak_2), main_terminus=jakcoords_tert[0])
Jakobshavn_tert.load_network(filename='Jakobshavn_tert.pickle')

kbcoords = Flowline_CSV('KogeBugt-mainline-w_width.csv', 1, has_width=True, flip_order=True)[0]
kb_line = Flowline(coords=kbcoords, index=0, name='Koge Bugt mainline', has_width=True)
KogeBugt = PlasticNetwork(name='Koge Bugt', init_type='Flowline', branches=(kb_line), main_terminus=kbcoords[0])
KogeBugt.load_network(filename='KogeBugt.pickle')

#
#### INTERSECTING LINES
helcoords_0, helcoords_1, helcoords_2 = Flowline_CSV('Helheim-network-w_width.csv', 3, has_width=True, flip_order=False)
hel_0 = Branch(coords=helcoords_0, index=0, order=0)
hel_1 = Branch(coords=helcoords_1, index=1, order=1, flows_to=0)
hel_2 = Branch(coords=helcoords_2, index=2, order=1, flows_to=0)
hel_branches = (hel_0, hel_1, hel_2)
Helheim = PlasticNetwork(name='Helheim', init_type='Branch', branches=hel_branches, main_terminus=helcoords_0[0])
Helheim.make_full_lines()
Helheim.load_network(filename='Helheim.pickle')

kangercoords_0, kangercoords_1, kangercoords_2, kangercoords_3, kangercoords_4 = Flowline_CSV('Documents/GitHub/plastic-networks/kangerlussuaq-network-w_width.csv', 5, has_width=True, flip_order=False)
kanger_0 = Branch(coords=kangercoords_0, index=0, order=0)
kanger_1 = Branch(coords=kangercoords_1, index=1, order=1, flows_to=0, intersect=174)
kanger_2 = Branch(coords=kangercoords_2, index=2, order=1, flows_to=0, intersect=191) #DIFFERENT FROM PREVIOUS BRANCH 2.  NEW FLOWLINE SET AS OF 31 MAR 2018
kanger_3 = Branch(coords=kangercoords_3, index=3, order=1, flows_to=0, intersect=146)
kanger_4 = Branch(coords=kangercoords_4, index=4, order=1, flows_to=0, intersect=61)
kanger_branches = (kanger_0, kanger_1, kanger_3, kanger_4)
Kanger = PlasticNetwork(name='Kangerlussuaq', init_type='Branch', branches=kanger_branches, main_terminus=kangercoords_0[0])
Kanger.make_full_lines()
Kanger.load_network(filename='Kangerlussuaq.pickle')


##-------------------
### PROCESSING LINE FUNCTIONS + OPTIMIZING YIELD
##-------------------

glacier_networks = (Jakobshavn_main, Jakobshavn_sec, Jakobshavn_tert, KogeBugt, Helheim, Kanger) #list which glaciers we're handling

for gl in glacier_networks:
    print gl.name
    gl.process_full_lines(B_interp, S_interp, H_interp)
    if gl in (KogeBugt, Kanger): # Add more sophisticated code to catch warnings?
        gl.remove_floating()
        if gl in (Kanger,):
            gl.make_full_lines()
        else:
            pass
        gl.process_full_lines(B_interp, S_interp, H_interp)
    #gl.optimize_network_yield(check_all=False)
    for fln in gl.flowlines:
        fln.yield_type = gl.network_yield_type
        fln.optimal_tau = gl.network_tau
    gl.network_ref_profiles()


###-------------------
#### FORWARD PROJECTION
###-------------------

testyears = arange(20)

jak_rate = 10/H0 #m/a, mid/low end from Csatho 2008
kb_rate = 2/H0
hel_rate = 5/H0
kanger_rate = 5/H0
thinrates = [jak_rate, kb_rate, hel_rate, kanger_rate]

Jakobshavn_main.network_time_evolve(testyears=testyears, thinrate=jak_rate)

KogeBugt.network_time_evolve(testyears=testyears, thinrate=kb_rate)

Helheim.network_time_evolve(testyears=testyears, thinrate=hel_rate)

Kanger.network_time_evolve(testyears=testyears, thinrate=kanger_rate)


   
###-------------------
#### REALISTIC FORCING
###-------------------

#
#obs_ish_projections = []
#for j, procglac in enumerate(processed_glaciers):
#    fwdmodel_dicts = PlasticEvol(procglac, testyears, thinrate = thinrates[j])
#    obs_ish_projections.append(fwdmodel_dicts)
#
#linear_increase_thinning = [np.linspace(tr, 2*tr, num=100) for tr in thinrates]
#lin_inc_projections = []
#for j, procglac in enumerate(processed_glaciers):
#    fwdmodel_dicts = PlasticEvol(procglac, arange(100), thinvalues=linear_increase_thinning[j])
#    lin_inc_projections.append(fwdmodel_dicts)
#
#   
###-------------------
#### SUMMARY PLOTTING
###-------------------
#
#names = ['Jakobshavn', 'Koge Bugt', 'Helheim', 'Kangerlussuaq']
##rates = ['10 m/a', '2 m/a', '5 m/a', '5 m/a'] #for constant 'obs-ish' projections
#rates = ['10-20 m/a', '2-4 m/a', '5-10 m/a', '5-10 m/a']
##styles = [':', '-.', '--', '-']
#markers = ['o', '^', 'd', '*']
#cmap = matplotlib.cm.get_cmap('winter')
#colors = cmap([0.3, 0.5, 0.7, 0.9])
#
#plt.figure()
##for j,pr in enumerate(lowend_projection_storage):
###for j,pr in enumerate(obs_ish_projections):
#for j,pr in enumerate(lin_inc_projections):
#    plt.plot(arange(100), -0.001*np.array(pr[0]['Termini'][1:]), linewidth=4, color=colors[j], label='{}, {}'.format(names[j], rates[j]))
#    plt.plot(arange(100)[::5], -0.001*np.array(pr[0]['Termini'][1:])[::5], linewidth=0, marker=markers[j], ms=10, color=colors[j])
#plt.legend(loc='lower left')
#plt.axes().set_xlabel('Year of simulation', size=30)
#plt.axes().set_ylabel('Terminus change [km]', size=30)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_ylim(-16, 1)
#plt.axes().set_yticks([-15, -10, -5, 0])
#plt.show()
