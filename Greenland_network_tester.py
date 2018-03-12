# Testing Flowline, Branch, and PlasticNetwork classes for comparison with function-based AGU results
# 9 Mar 2018  EHU

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
#gl_bed_path ='Documents/1. Research/2. Flowline networks/Model/Data/BedMachine-Greenland/MCdataset-2015-04-27.nc'
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

jakcoords_0, jakcoords_1, jakcoords_2 = Flowline_CSV('Documents/1. Research/2. Flowline networks/Model/Greenland_tests/Flowline_sets/jakobshavn-thead_lines.csv', 3)
jak_0 = Flowline(coords=jakcoords_0, index=0, name='Jak mainline')
jak_1 = Flowline(coords=jakcoords_1, index=1, name='Jak line 1')
jak_2 = Flowline(coords=jakcoords_2, index=2, name='Jak line 2')
Jakobshavn_main = PlasticNetwork(name='Jakobshavn Isbrae [main]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_0[0])
Jakobshavn_1 = PlasticNetwork(name='Jakobshavn Isbrae 1', init_type='Flowline', branches=(jak_1), main_terminus=jakcoords_1[0])
Jakobshavn_2 = PlasticNetwork(name='Jakobshavn Isbrae 2', init_type='Flowline', branches=(jak_2), main_terminus=jakcoords_2[0])

kbcoords = Flowline_CSV('Documents/1. Research/2. Flowline networks/Model/Greenland_tests/Flowline_sets/koge_bugt-thead_lines.csv', 1)[0]
kb_line = Flowline(coords=kbcoords, index=0, name='Koge Bugt mainline')
KogeBugt = PlasticNetwork(name='Koge Bugt', init_type='Flowline', branches=kb_line, main_terminus=kbcoords[0])


### INTERSECTING LINES - ADDITIONAL PROCESSING NEEDED
helcoords_0, helcoords_1, helcoords_2 = Flowline_CSV('Documents/1. Research/2. Flowline networks/Model/Greenland_tests/Flowline_sets/helheim-branch_lines.csv', 3)
hel_0 = Branch(coords=helcoords_0, index=0, order=0)
hel_1 = Branch(coords=helcoords_1, index=1, order=1, flows_to=0)
hel_2 = Branch(coords=helcoords_2, index=2, order=1, flows_to=0)
hel_branches = (hel_0, hel_1, hel_2)
Helheim = PlasticNetwork(name='Helheim', init_type='Branch', branches=hel_branches, main_terminus=helcoords_0[0])
Helheim.make_full_lines()


kangercoords_0, kangercoords_1, kangercoords_2, kangercoords_3, kangercoords_4, kangercoords_5 = Flowline_CSV('Documents/1. Research/2. Flowline networks/Model/Greenland_tests/Flowline_sets/kangerlussuaq-branch_lines.csv', 6)
kanger_0 = Branch(coords=kangercoords_0, index=0, order=0)
kanger_1 = Branch(coords=kangercoords_1, index=1, order=1, flows_to=0)
#kanger_2 = Branch(coords=kangercoords_2, index=2, order=2, flows_to=1) #branch 2 did not connect under previous methodology.  Keep low expectations here
kanger_3 = Branch(coords=kangercoords_3, index=3, order=1, flows_to=0)
kanger_4 = Branch(coords=kangercoords_4, index=4, order=1, flows_to=0)
kanger_branches = (kanger_0, kanger_1, kanger_3, kanger_4)
Kanger = PlasticNetwork(name='Kangerlussuaq', init_type='Branch', branches=kanger_branches, main_terminus=kangercoords_0[0])
Kanger.make_full_lines()

#separate network for kanger_5, which has separate terminus
kanger_secondary = Flowline(coords=kangercoords_5, index=0)
Kanger_s = PlasticNetwork(name='Kangerlussuaq [secondary]', init_type='Flowline', branches=kanger_secondary, main_terminus=kangercoords_5[0])


##-------------------
### PROCESSING LINE FUNCTIONS
##-------------------
for gl in (Jakobshavn_main, Jakobshavn_1, Jakobshavn_2, KogeBugt, Helheim, Kanger, Kanger_s):
    gl.process_full_lines(B_interp, S_interp, H_interp)
    if gl in (Kanger, Kanger_s): # Add more sophisticated code to catch warnings?
        gl.remove_floating()
        gl.make_full_lines()
        gl.process_full_lines(B_interp, S_interp, H_interp)
#    gl.optimize_network_yield(check_all=True)


##debugging Kangerlussuaq
#for gl in (Kanger, Kanger_s):
#    gl.process_full_lines(B_interp, S_interp, H_interp)
#    gl.remove_floating()
#    gl.make_full_lines()
#    gl.process_full_lines(B_interp, S_interp, H_interp)
#    gl.optimize_network_yield(check_all=False)

#fltest = []
#kanger_arcmain = ArcArray(Kanger.flowlines[0].coords)
#for pt in kanger_arcmain:
#    gimpthick = Kanger.flowlines[0].thickness_function(pt)
#    bmbed = Kanger.flowlines[0].bed_function(pt)
#    flthick = FlotationThick(bmbed)
#    fltest.append(gimpthick-flthick)
#nonfloating = argwhere(sign(fltest)>0) #tops flotation after pt 18--hard-coding for now (Dec 13 2017)



##-------------------
### TY OPTIMIZATION
##-------------------



###-------------------
#### OPTIMAL TY - HARD CODED FROM SEPARATE OPTIMIZATION FOR NOW
###-------------------
#bestfit_jak = [(325000.0, 265000.0), (255000.0, 220000.0), (415000.0, 390000.0)]
#bestfit_kb = [(340e3, 295e3)]
#bestfit_hel = [(330000.0, 270000.0), (150000.0, 65000.0), (265000.0, 205000.0)]
##bestfit_kanger = [(315000.0, 260000.0),
## (245000.0, 195000.0),
## (220000.0, 135000.0), #Starts at intersection rather than terminus--use with caution
## (235000.0, 180000.0),
## (210000.0, 170000.0),
## (50000.0, 50000.0)] #buggy - thins below flotation even though GIMP thickness exceeds flotation thickness
#bestfit_kanger_limited =  [(315000.0, 260000.0),
# (245000.0, 195000.0),
# (235000.0, 180000.0),
# (210000.0, 170000.0)]
#
#all_tyfits = [bestfit_jak, bestfit_kb, bestfit_hel, bestfit_kanger_limited]
#
#
###-------------------
#### PROCESSING GLACIERS
###-------------------
#keylist = ('bed', 'surface', 'thickness')
#fieldlist = (B_interp, S_interp, H_interp)
#
#all_glaciers = [jak_master, kb_master, hel_master, kanger_limited] #ADD IV IF TIME
#processed_glaciers = []
#for i, g in enumerate(all_glaciers):
#    bf_taus = all_tyfits[i]
#    processed = ProcessDicts(g, keylist, fieldlist, bf_taus)
#    processed_glaciers.append(processed)
#
###-------------------
#### FORWARD PROJECTION
###-------------------
#
#testyears = arange(100)
#
##jak_proj = PlasticEvol(processed_glaciers[0], testyears)
##kb_proj = PlasticEvol(processed_glaciers[1], testyears)
##hel_proj = PlasticEvol(processed_glaciers[2], testyears)
##kanger_proj = PlasticEvol(processed_glaciers[3], testyears)
#
#
#projection_storage = []
#for procglac in processed_glaciers:
#    fwdmodel_dicts = PlasticEvol(procglac, testyears)
#    projection_storage.append(fwdmodel_dicts)
#
#   
###-------------------
#### REALISTIC FORCING
###-------------------
#jak_rate = 10/H0 #m/a, mid/low end from Csatho 2008
#kb_rate = 2/H0
#hel_rate = 5/H0
#kanger_rate = 5/H0
#thinrates = [jak_rate, kb_rate, hel_rate, kanger_rate]
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
