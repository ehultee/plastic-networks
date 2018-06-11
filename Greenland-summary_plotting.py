## Loading in output from scenarios to make summary plots
## 9 Jun 2018  EHU

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

#Smoothing bed to check effect on dLdt
unsmoothB = B
smoothB = gaussian_filter(B, 2)
#B_processed = np.ma.masked_where(thick_mask !=2, smoothB)

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], S.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])

print 'Reading in 5-year surface elevation change'
gl_sec_path ='Documents/GitHub/plastic-networks/Data/CS2-SEC_5yr.nc'
#gl_sec_path ='Documents/GitHub/plastic-networks/Data/CS2-SEC_2yr.nc'
fh3 = Dataset(gl_sec_path, mode='r')
x_sec = fh3.variables['x'][:].copy() #x-coord (polar stereo)
y_sec = fh3.variables['y'][:].copy() #y-coord (polar stereo)
t_sec = fh3.variables['t'][:].copy() #average time of slice (days since 1 JAN 2000)
sec_raw = fh3.variables['SEC'][:].copy()
fh3.close()

sec_i_masked = np.ma.masked_greater(sec_raw[:,:,0], 9000)
sec_i_excludemasked = np.ma.filled(sec_i_masked, fill_value=np.mean(sec_i_masked))
#sec_i_regrid = interpolate.griddata((x_sec.ravel(), y_sec.ravel()), sec_i_masked.ravel(), (Xmat, Ymat), method='nearest')
SEC_i = interpolate.RectBivariateSpline(x_sec, y_sec, sec_i_excludemasked.T)


##-------------------
### READING IN GLACIERS
##-------------------


jakcoords_main = Flowline_CSV('Documents/GitHub/plastic-networks/Data/jakobshavn-mainline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_0 = Flowline(coords=jakcoords_main, index=0, name='Jak mainline', has_width=True)
Jak_main_warmbaseline = PlasticNetwork(name='Jakobshavn Isbrae [main/south]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_main[0])
Jak_main_warmbaseline.load_network(filename='JakobshavnIsbrae-main_south-9Jun18-0forcing-warmice-dt025a.pickle', load_mainline_output=True)
Jak_main_coldbaseline = PlasticNetwork(name='Jakobshavn Isbrae [main/south]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_main[0])
Jak_main_coldbaseline.load_network(filename='JakobshavnIsbrae-main_south-8Jun18-0forcing-coldice-dt025a.pickle', load_mainline_output=True)

jakcoords_sec = Flowline_CSV('Jakobshavn_secondary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_1 = Flowline(coords=jakcoords_sec, index=0, name='Jak secondary branch', has_width=True)
Jak_sec_warmbaseline = PlasticNetwork(name='Jakobshavn Isbrae [secondary/central]', init_type='Flowline', branches=(jak_1), main_terminus=jakcoords_sec[0])
Jak_sec_warmbaseline.load_network(filename='JakobshavnIsbrae-secondary_central-9Jun18-0forcing-warmice-dt025a.pickle', load_mainline_output=True)
Jak_sec_coldbaseline = PlasticNetwork(name='Jakobshavn Isbrae [secondary/central]', init_type='Flowline', branches=(jak_1), main_terminus=jakcoords_sec[0])
Jak_sec_coldbaseline.load_network(filename='JakobshavnIsbrae-secondary_central-8Jun18-0forcing-coldice-dt025a.pickle', load_mainline_output=True)

jakcoords_tert = Flowline_CSV('Jakobshavn_tertiary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_2 = Flowline(coords=jakcoords_tert, index=0, name='Jak tertiary branch', has_width=True)
Jak_tert_warmbaseline = PlasticNetwork(name='Jakobshavn Isbrae [tertiary/north]', init_type='Flowline', branches=(jak_2), main_terminus=jakcoords_tert[0])
Jak_tert_warmbaseline.load_network(filename='JakobshavnIsbrae-tertiary_north-9Jun18-0forcing-warmice-dt025a.pickle', load_mainline_output=True)
Jak_tert_coldbaseline = PlasticNetwork(name='Jakobshavn Isbrae [tertiary/north]', init_type='Flowline', branches=(jak_2), main_terminus=jakcoords_tert[0])
Jak_tert_coldbaseline.load_network(filename='JakobshavnIsbrae-tertiary_north-8Jun18-0forcing-coldice-dt025a.pickle', load_mainline_output=True)

kbcoords = Flowline_CSV('KogeBugt-mainline-w_width.csv', 1, has_width=True, flip_order=True)[0]
kb_line = Flowline(coords=kbcoords, index=0, name='Koge Bugt mainline', has_width=True)
KogeBugt_warmbaseline = PlasticNetwork(name='Koge Bugt', init_type='Flowline', branches=(kb_line), main_terminus=kbcoords[0])
KogeBugt_warmbaseline.load_network(filename='KogeBugt-9Jun18-0forcing-warmice-dt025a.pickle', load_mainline_output=True)
KogeBugt_coldbaseline = PlasticNetwork(name='Koge Bugt', init_type='Flowline', branches=(kb_line), main_terminus=kbcoords[0])
KogeBugt_coldbaseline.load_network(filename='KogeBugt-8Jun18-0forcing-coldice-dt025a.pickle', load_mainline_output=True)


### INTERSECTING LINES
helcoords_0, helcoords_1, helcoords_2 = Flowline_CSV('Helheim-network-w_width.csv', 3, has_width=True, flip_order=False)
hel_0 = Branch(coords=helcoords_0, index=0, order=0)
hel_1 = Branch(coords=helcoords_1, index=1, order=1, flows_to=0)
hel_2 = Branch(coords=helcoords_2, index=2, order=1, flows_to=0)
hel_branches = (hel_0, hel_1, hel_2)
#Helheim.make_full_lines()
Helheim_warmbaseline = PlasticNetwork(name='Helheim', init_type='Branch', branches=hel_branches, main_terminus=helcoords_0[0])
Helheim_warmbaseline.load_network(filename='Helheim-9Jun18-0forcing-warmice-dt025a.pickle', load_mainline_output=True, load_tributary_output=True)
Helheim_coldbaseline = PlasticNetwork(name='Helheim', init_type='Branch', branches=hel_branches, main_terminus=helcoords_0[0])
Helheim_coldbaseline.load_network(filename='Helheim-8Jun18-0forcing-coldice-dt025a.pickle', load_mainline_output=True, load_tributary_output=True)

kangercoords_0, kangercoords_1, kangercoords_2, kangercoords_3, kangercoords_4 = Flowline_CSV('Documents/GitHub/plastic-networks/Data/kangerlussuaq-network-w_width.csv', 5, has_width=True, flip_order=False)
kanger_0 = Branch(coords=kangercoords_0, index=0, order=0)
kanger_1 = Branch(coords=kangercoords_1, index=1, order=1, flows_to=0, intersect=174)
kanger_2 = Branch(coords=kangercoords_2, index=2, order=1, flows_to=0, intersect=191) #DIFFERENT FROM PREVIOUS BRANCH 2.  NEW FLOWLINE SET AS OF 31 MAR 2018
kanger_3 = Branch(coords=kangercoords_3, index=3, order=1, flows_to=0, intersect=146)
kanger_4 = Branch(coords=kangercoords_4, index=4, order=1, flows_to=0, intersect=61)
kanger_branches = (kanger_0, kanger_1, kanger_3, kanger_4)
#Kanger.make_full_lines()
Kanger_coldbaseline = PlasticNetwork(name='Kangerlussuaq', init_type='Branch', branches=kanger_branches, main_terminus=kangercoords_0[0])
Kanger_coldbaseline.load_network(filename='Kangerlussuaq-8Jun18-0forcing-coldice-dt025a.pickle', load_mainline_output=True, load_tributary_output=True)
Kanger_warmbaseline = PlasticNetwork(name='Kangerlussuaq', init_type='Branch', branches=kanger_branches, main_terminus=kangercoords_0[0])
Kanger_warmbaseline.load_network(filename='Kangerlussuaq-9Jun18-0forcing-warmice-dt025a.pickle', load_mainline_output=True, load_tributary_output=True)


###-------------------
#### PROJECTION SUMMARY
###-------------------
testyears = arange(50, step=0.25)

#Warm baseline fluxes
#Kanger_multibranch_flux = [Kanger_warmbaseline.model_output[j]['Terminus_flux'] for j in range(len(Kanger.flowlines))]
#Kanger_total_flux = sum(Kanger_multibranch_flux, axis = 0) #note that Kanger_multibranch_flux is multidimensional, needs care in summing
#Helheim_multibranch_flux = [Helheim_warmbaseline.model_output[j]['Terminus_flux'] for j in range(len(Helheim.flowlines))]
#Helheim_total_flux = sum(Helheim_multibranch_flux, axis=0)
#Jakobshavn_multibranch_flux = [Jak_main_warmbaseline.model_output[0]['Terminus_flux'], Jak_sec_warmbaseline.model_output[0]['Terminus_flux'], Jak_tert_warmbaseline.model_output[0]['Terminus_flux']]
#Jakobshavn_total_flux = sum(Jakobshavn_multibranch_flux, axis=0)
#fluxes = [Jak_main_warmbaseline.model_output[0]['Terminus_flux'], Jak_sec_warmbaseline.model_output[0]['Terminus_flux'], Jak_tert_warmbaseline.model_output[0]['Terminus_flux'], KogeBugt_warmbaseline.model_output[0]['Terminus_flux'], Helheim_warmbaseline.model_output[0]['Terminus_flux'], Kanger_warmbaseline.model_output[0]['Terminus_flux']]
#total_fluxes = [Jakobshavn_total_flux, KogeBugt_warmbaseline.model_output[0]['Terminus_flux'], Helheim_total_flux, Kanger_total_flux]
#NOTE NEED TO INCLUDE FLUXES FROM ALL LINES OF HELHEIM, KANGER


#Cold baseline fluxes for plot
Kanger_multibranch_flux_cold = [np.absolute(Kanger_coldbaseline.model_output[j]['Terminus_flux']) for j in range(len(Kanger.flowlines))]
Kanger_total_flux_cold = sum(Kanger_multibranch_flux_cold, axis = 0) #note that Kanger_multibranch_flux is multidimensional, needs care in summing
Helheim_multibranch_flux_cold = [np.absolute(Helheim_coldbaseline.model_output[j]['Terminus_flux']) for j in range(len(Helheim.flowlines))]
Helheim_total_flux_cold = sum(Helheim_multibranch_flux_cold, axis=0)
total_fluxes_split_CB = [Jak_main_coldbaseline.model_output[0]['Terminus_flux'], Jak_sec_coldbaseline.model_output[0]['Terminus_flux'], Jak_tert_coldbaseline.model_output[0]['Terminus_flux'], KogeBugt_coldbaseline.model_output[0]['Terminus_flux'], Helheim_total_flux_cold, Kanger_total_flux_cold] #totals from all termini, with Jak split by branch (network)

fluxes_cleaned = []
sle = [] #will be array of annual sea level contributions
for flux in total_fluxes_split_CB:
    flux_c = np.nan_to_num(flux)
    flux_abs = np.absolute(flux_c) #quick fix 9 Jun for some lines that have negative width upstream.  Fix more substantially later.
    fluxes_cleaned.append(flux_abs)
    sleq = (1E-12)*np.array(flux_abs)/(361.8) #Gt ice / mm sea level equiv conversion
    sle.append(sleq)
cumul_sle_pernetwork = []
total_sle = []
for sl in sle:
    c = np.cumsum(sl)
    cumul_sle_pernetwork.append(c)
total_sle = np.cumsum(cumul_sle_pernetwork, axis=0)

projections_WB = [Jak_main_warmbaseline.model_output, Jak_sec_warmbaseline.model_output, Jak_tert_warmbaseline.model_output, KogeBugt_warmbaseline.model_output, Helheim_warmbaseline.model_output, Kanger_warmbaseline.model_output]
projections_CB = [Jak_main_coldbaseline.model_output, Jak_sec_coldbaseline.model_output, Jak_tert_coldbaseline.model_output, KogeBugt_coldbaseline.model_output, Helheim_coldbaseline.model_output, Kanger_coldbaseline.model_output]
names = ['Sermeq Kujalleq [main]', 'Sermeq Kujalleq [central]', 'Sermeq Kujalleq [north]', 'Koge Bugt', 'Helheim', 'Kangerlussuaq']
combined_networks = ['Sermeq Kujalleq', 'Koge Bugt', 'Helheim', 'Kangerlussuaq']
##styles = [':', '-.', '--', '-']
markers = ['o', '.', ',', '^', 'd', '*']
styles = ['-', ':', '-.', '-', '-', '-']
cmap = matplotlib.cm.get_cmap('winter')
colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
alt_colors = matplotlib.cm.get_cmap('Greys')([0.3, 0.5, 0.7, 0.9])

#terminus
plt.figure()
for j, pr in enumerate(projections_WB):
    print j
    plt.plot(testyears, -0.001*np.array(pr[0]['Termini'][1::]), linewidth=4, color=colors[j], linestyle=styles[j], label='{}'.format(names[j]))
    plt.plot(testyears[::20], -0.001*np.array(pr[0]['Termini'][1::])[::20], linewidth=0, marker=markers[j], ms=10, color=colors[j])
plt.legend(loc='lower left')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Terminus change [km]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
plt.axes().set_ylim(-100, 1)
plt.axes().set_yticks([-75, -50, -25, 0])
plt.title('50 a (dt=0.25 a) dynamic change only - Gaussian sigma=2 - A=1.7E-24', fontsize=20)
plt.show()


##SINGLE NETWORK - splitting termini
plt.figure('Kangerlussuaq terminus change')
for k in range(len(Kanger_warmbaseline.model_output)): #for each branch j
    #colork = matplotlib.cm.get_cmap('viridis')(k/len(Kanger_warmbaseline.model_output))
    mo = Kanger_warmbaseline.model_output[k] #for some reason enumerate doesn't work with loaded-in output, so we're stuck with this
    markerk = (k+2, mod(k+1, 3), 0)
    plt.plot(testyears, -0.001*np.array(mo['Termini'][0:-1:]), linewidth=4, color='k')
    plt.plot(testyears[::10], -0.001*np.array(mo['Termini'][0:-1:])[::10], linewidth=0, color='k', marker=markerk, ms=10)
#plt.legend(loc='lower left')
plt.axes().set_xlabel('Year of simulation', size=30)
plt.axes().set_ylabel('Terminus change [km]', size=30)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
plt.axes().set_ylim(-50, 1)
plt.axes().set_yticks([-50, -25, 0])
plt.show()

###Flux
plt.figure()
for j in range(len(names)):
    plt.plot(testyears, 1E-12*np.array(fluxes[j]), linewidth=4, linestyle=styles[j], color=colors[j], label=names[j])
    plt.plot(testyears[::50], 1E-12*np.array(fluxes[j][::50]), linewidth=0, marker=markers[j], ms=10, color=colors[j])
    plt.fill_between(testyears, y1=1E-12*np.array(fluxes[j]), y2=0, color=colors[j], alpha=0.5)    
plt.legend(loc='upper right')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Terminus ice flux [Gt/a]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
plt.axes().set_xlim(0,50)
plt.axes().set_xticks([0,10,20, 30, 40, 50])
plt.axes().set_ylim(0, 401)
plt.axes().set_yticks([0, 50, 100, 150, 200, 250, 300, 350, 400])
plt.show()

####Sea level equivalent
plt.figure()
for j in range(len(names)):
    plt.plot(1+testyears[::], total_sle[j], linewidth=4, color=colors[::][j], label=names[j])
    plt.plot(1+testyears[::5], total_sle[j][::5], linewidth=0, marker=markers[j], ms=10, color=colors[::][j])
    if j==0:
        plt.fill_between(1+testyears[::], y1=total_sle[j], y2=0, color=colors[::][j], alpha=0.7)  
    else:
        plt.fill_between(1+testyears[::], y1=total_sle[j], y2=total_sle[j-1], color=colors[::][j], alpha=0.7)      
plt.legend(loc='upper left')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
plt.axes().set_xlim(0, 50)
plt.axes().set_xticks([0, 10, 20, 30, 40, 50])
#plt.axes().set_ylim(-16, 1)
#plt.axes().set_yticks([0, 1, 2, 3, 4])
plt.show()