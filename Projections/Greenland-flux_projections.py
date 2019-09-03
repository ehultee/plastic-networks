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
## Special import for SERMeQ modules
import sys
sys.path.insert(0, 'Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *


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

#Smoothing bed
unsmoothB = B
smoothB = gaussian_filter(B, 2)
smoothS = gaussian_filter(S, 2) #17 Jan 19 - smoothing S as well for consistency with auto-selected networks
#B_processed = np.ma.masked_where(thick_mask !=2, smoothB)

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothS.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])

### Reading in SENTINEL velocity map
#print 'Now reading in (vector) velocity map'
#v_path = 'Documents/1. Research/2. Flowline networks/Model/Data/ESA-Greenland/greenland_iv_500m_s1_20161223_20170227_v1_0.nc'
#fh2 = Dataset(v_path, mode='r')
#xv = fh2.variables['x'][:].copy()
#yv = fh2.variables['y'][:].copy()
##yv = yv_flipped[::-1]
#v_raw = fh2.variables['land_ice_surface_velocity_magnitude'][:].copy() #this is v(y, x)
#vx_raw = fh2.variables['land_ice_surface_easting_velocity'][:].copy()
#vy_raw =fh2.variables['land_ice_surface_northing_velocity'][:].copy()
#v_upper = np.ma.masked_greater(v_raw, 10000)
#vx_upper = np.ma.masked_greater(vx_raw, 10000)
#vy_upper = np.ma.masked_greater(vy_raw, 10000)
#fh2.close()
#
### Interpolate SENTINEL and sample at BedMachine points
#print 'Now interpolating to same grid'
#vf_x = interpolate.interp2d(xv, yv[::-1], vx_upper[::-1,::])
#vf_y = interpolate.interp2d(xv, yv[::-1], vy_upper[::-1,::])
#vf = interpolate.interp2d(xv, yv[::-1], v_upper[::-1, ::])

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
### FINDING GLACIERS
##-------------------

jakcoords_main = Flowline_CSV('Documents/GitHub/plastic-networks/Data/jakobshavn-mainline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_0 = Flowline(coords=jakcoords_main, index=0, name='Jak mainline', has_width=True)
Jakobshavn_main = PlasticNetwork(name='Jakobshavn Isbrae [main/south]', init_type='Flowline', branches=(jak_0), main_terminus=jakcoords_main[0])
Jakobshavn_main.load_network(filename='JakobshavnIsbrae-main_south.pickle')

jakcoords_sec = Flowline_CSV('Jakobshavn_secondary-flowline-w_width.csv', 1, has_width=True, flip_order=True)[0]
jak_1 = Flowline(coords=jakcoords_sec, index=0, name='Jak secondary branch', has_width=True)
Jakobshavn_sec = PlasticNetwork(name='Jakobshavn Isbrae [secondary/central]', init_type='Flowline', branches=(jak_1), main_terminus=jakcoords_sec[0])
Jakobshavn_sec.load_network(filename='Jakobshavn_sec.pickle')

jakcoords_tert = Flowline_CSV('Jakobshavn_tertiary-flowline-w_width.csv', 1, has_width=True, flip_order=True)[0]
jak_2 = Flowline(coords=jakcoords_tert, index=0, name='Jak tertiary branch', has_width=True)
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

kangercoords_0, kangercoords_1, kangercoords_2, kangercoords_3, kangercoords_4 = Flowline_CSV('Documents/GitHub/plastic-networks/Data/kangerlussuaq-network-w_width.csv', 5, has_width=True, flip_order=False)
kanger_0 = Branch(coords=kangercoords_0, index=0, order=0)
kanger_1 = Branch(coords=kangercoords_1, index=1, order=1, flows_to=0, intersect=174)
#kanger_2 = Branch(coords=kangercoords_2, index=2, order=1, flows_to=0, intersect=191) #DIFFERENT FROM PREVIOUS BRANCH 2.  NEW FLOWLINE SET AS OF 31 MAR 2018
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
#glacier_networks=(Jakobshavn_main,)

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
        fln.yield_type  = gl.network_yield_type
        fln.optimal_tau = gl.network_tau
    gl.network_ref_profiles()
#
#
####-------------------
##### FORWARD PROJECTION--FORCING FROM TERMINUS
####-------------------
##
testyears = arange(20, step=0.25)
db=False
##
#testglac = (Jakobshavn_main,)
testglac = glacier_networks #to test all
###testglac = (Kanger,)
##
#Finding SEC rates and making persistence projection
for gl in testglac:
    print gl.name
    gl.sec_mainline = np.asarray([SEC_i(gl.flowlines[0].coords[i,0], gl.flowlines[0].coords[i,1]) for i in range(len(gl.flowlines[0].coords))])
    away_from_edge = np.argmin(gl.sec_mainline)
    gl.sec_alphadot = np.mean(gl.sec_mainline[away_from_edge::])
    variable_forcing = linspace(start=gl.sec_alphadot, stop=2*gl.sec_alphadot, num=len(testyears))
    gl.terminus_sec = float(min(gl.sec_mainline.flatten()))#using min because values close to edge get disrupted by mask interpolation
    gl.terminus_time_evolve(testyears=testyears, alpha_dot=gl.sec_alphadot, dL=1/L0, separation_buffer=10000/L0, has_smb=True, terminus_balance=gl.terminus_sec, submarine_melt = 0, debug_mode=db, rate_factor=3.7E-26) 
#    
#    print 'Saving output for {}'.format(gl.name)
#    fn = str(gl.name)
#    fn1 = fn.replace(" ", "")
#    fn2 = fn1.replace("[", "-")
#    fn3 = fn2.replace("/", "_")
#    fn4 = fn3.replace("]", "")
#    fn5 = fn4+'-2Jul18-persistence-coldice-100a_dt025a.pickle'
#    gl.save_network(filename=fn5)

Kanger_multibranch_flux = [Kanger.model_output[j]['Terminus_flux'] for j in range(len(Kanger.flowlines))]
Kanger_total_flux = sum(Kanger_multibranch_flux, axis = 0) #note that Kanger_multibranch_flux is multidimensional, needs care in summing
Helheim_multibranch_flux = [Helheim.model_output[j]['Terminus_flux'] for j in range(len(Helheim.flowlines))]
Helheim_total_flux = sum(Helheim_multibranch_flux, axis=0)
Jakobshavn_multibranch_flux = [Jakobshavn_main.model_output[0]['Terminus_flux'], Jakobshavn_sec.model_output[0]['Terminus_flux'], Jakobshavn_tert.model_output[0]['Terminus_flux']]
Jakobshavn_total_flux = sum(Jakobshavn_multibranch_flux, axis=0)
main_fluxes = [Jakobshavn_main.model_output[0]['Terminus_flux'], Jakobshavn_sec.model_output[0]['Terminus_flux'], Jakobshavn_tert.model_output[0]['Terminus_flux'], KogeBugt.model_output[0]['Terminus_flux'], Helheim.model_output[0]['Terminus_flux'], Kanger.model_output[0]['Terminus_flux']] #main calving termini only
total_fluxes = [Jakobshavn_total_flux, KogeBugt.model_output[0]['Terminus_flux'], Helheim_total_flux, Kanger_total_flux]
total_fluxes_split = [Jakobshavn_main.model_output[0]['Terminus_flux'], Jakobshavn_sec.model_output[0]['Terminus_flux'], Jakobshavn_tert.model_output[0]['Terminus_flux'], KogeBugt.model_output[0]['Terminus_flux'], Helheim_total_flux, Kanger_total_flux] #totals from all termini, with Jak split by branch (network)

fluxes_cleaned = []
sle = [] #will be array of annual sea level contributions
for flux in total_fluxes_split:
    flux_c = np.nan_to_num(flux)
    #flux_a = np.absolute(flux_c)
    fluxes_cleaned.append(flux_c)
    sleq = (1E-12)*np.array(flux_c)/(361.8) #Gt ice / mm sea level equiv conversion
    sle.append(sleq)
cumul_sle_pernetwork = []
#total_sle = []
for sl in sle:
    c = np.cumsum(sl)
    cumul_sle_pernetwork.append(c)
total_sle = np.cumsum(cumul_sle_pernetwork, axis=0)

#### Splitting by branch for Jakobshavn
##fluxes_cleaned = []
##sle = []
##for flux in total_fluxes_split:
##    flux_c = np.nan_to_num(flux)
##    flux_a = np.absolute(flux_c)
##    fluxes_cleaned.append(flux_a)
##    sleq = (1E-12)*np.array(flux_a)/(361.8)
##    sle.append(sleq)
##cumul_sle_split = []
##for sl in sle:
##    c = np.cumsum(sl)
##    cumul_sle_split.append(c)
##total_sle_split = np.cumsum(cumul_sle_split, axis=0)
###
#######-------------------
######## SUMMARY PLOTTING
#######-------------------
###
#projections = [Jakobshavn_main.model_output, Jakobshavn_sec.model_output, Jakobshavn_tert.model_output, KogeBugt.model_output, Helheim.model_output, Kanger.model_output]
#names = ['Sermeq Kujalleq [main]', 'Sermeq Kujalleq [central]', 'Sermeq Kujalleq [north]', 'Koge Bugt', 'Helheim', 'Kangerlussuaq']
#combined_networks = ['Sermeq Kujalleq', 'Koge Bugt', 'Helheim', 'Kangerlussuaq']
##rates = ['{0:.2f} m/a'.format(gl.sec_alphadot) for gl in glacier_networks]
###styles = [':', '-.', '--', '-']
#markers = ['o', '.', ',', '^', 'd', '*']
#styles = ['-', ':', '-.', '-', '-', '-']
#cmap = cm.get_cmap('winter')
#colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
#alt_colors = matplotlib.cm.get_cmap('Greys')([0.3, 0.5, 0.7, 0.9])
##
####
#####terminus
#plt.figure('Terminus retreat, 100 a -30C ice, dynamic drawdown', figsize=(12,8))
#for j, pr in enumerate(projections):
#    print j
#    plt.plot(testyears, -0.001*np.array(pr[0]['Termini'][1::]), linewidth=4, color=colors[j], linestyle=styles[j], label='{}'.format(names[j]))
#    plt.plot(testyears[::20], -0.001*np.array(pr[0]['Termini'][1::])[::20], linewidth=0, marker=markers[j], ms=10, color=colors[j])
##plt.legend(loc='lower left')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Terminus change [km]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_xlim(0, 82)
##plt.axes().set_xticks([7, 32, 57, 82])
##plt.axes().set_xticklabels(['2025', '2050', '2075', '2100']) #label projection out to calendar year 2100
##plt.axes().set_xlim(0, 82)
#plt.axes().set_xticks([0, 25, 50, 75, 100])
##plt.axes().set_ylim(-100, 1)
#plt.axes().set_yticks([-40, -30, -20, -10, 0])
#plt.title('100 yr dynamic drawdown, -30C ice', fontsize=26)
#plt.show()
####
#######JAKOBSHAVN ONLY - checking monotonicity
######plt.figure('Koge Bugt terminus change')
#######for k, mo in enumerate(Helheim.model_output): #for each branch j
######mo = KogeBugt.model_output[0]
######plt.plot(0.1*arange(len(mo['Termini'])), -0.001*np.array(mo['Termini']), linewidth=4, color='k', label='Koge Bugt mainline')
######plt.plot(0.1*arange(len(mo['Termini']))[::5], -0.001*np.array(mo['Termini'])[::5], linewidth=0, color='k', marker='.', ms=10)
######plt.legend(loc='lower left')
######plt.axes().set_xlabel('Year of simulation', size=30)
######plt.axes().set_ylabel('Terminus change [km]', size=30)
######plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#######plt.axes().set_ylim(-16, 1)
#######plt.axes().set_yticks([-15, -10, -5, 0])
######plt.show()
#####
####SINGLE NETWORK - splitting termini
#plt.figure('Kanger terminus change')
#for k, mo in enumerate(Kanger.model_output): #for each branch j
#    #colork = matplotlib.cm.get_cmap('viridis')(k/len(Helheim.model_output))
#    #markerk = (k+2, mod(k+1, 3), 0)
#    lsk = (':', '-.', '--', '-')
#    plt.plot(testyears[::], -0.001*np.array(mo['Termini'])[:-1:], linewidth=4, color='k', ls=lsk[k], label='Kanger line {}'.format(k))
#    #plt.plot(testyears[::10], -0.001*np.array(mo['Termini'])[:-1:10], linewidth=0, color='k', marker=markerk, ms=10, label='Helheim line {}'.format(k))
##plt.legend(loc='upper right')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Terminus change [km]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0, 20)
#plt.axes().set_xticks([0, 5, 10, 15, 20])
#plt.axes().set_xticklabels(['0', '', '', '', '20'])
##plt.axes().set_ylim(-16, 1)
##plt.axes().set_yticks([-15, -10, -5, 0])
#plt.show()
##
####Flux
##plt.figure(figsize=(12,8))
##for j in range(len(names)):
##    plt.plot(testyears, 1E-12*np.array(fluxes_cleaned[j]), linewidth=4, linestyle=styles[j], color=colors[j], label=names[j])
##    plt.plot(testyears[::50], 1E-12*np.array(fluxes_cleaned[j][::50]), linewidth=0, marker=markers[j], ms=10, color=colors[j])
##    plt.fill_between(testyears, y1=1E-12*np.array(fluxes_cleaned[j]), y2=0, color=colors[j], alpha=0.5)    
###plt.legend(loc='upper right')
##plt.axes().set_xlabel('Year of simulation', size=20)
##plt.axes().set_ylabel('Terminus ice flux [Gt/a]', size=20)
##plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
###plt.axes().set_xlim(0,15.5)
###plt.axes().set_xticks([0,5,10, 15])
###plt.axes().set_ylim(0, 18)
###plt.axes().set_yticks([0, 5, 10, 15])
##plt.show()
########
#############Sea level equivalent
#plt.figure(figsize=(12,8))
#for j in range(len(names)):
#    plt.plot(testyears[::], total_sle[j], linewidth=4, color=colors[::][j], label=names[j])
#    plt.plot(testyears[::5], total_sle[j][::5], linewidth=0, marker=markers[j], ms=10, color=colors[::][j])
#    if j==0:
#        plt.fill_between(testyears[::], y1=total_sle[j], y2=0, color=colors[::][j], alpha=0.7)  
#    else:
#        plt.fill_between(testyears[::], y1=total_sle[j], y2=total_sle[j-1], color=colors[::][j], alpha=0.7)      
##plt.legend(loc='upper left')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0, 100)
#plt.axes().set_xticks([0, 25, 50, 75, 100])
##plt.axes().set_ylim(-16, 1)
##plt.axes().set_yticks([0, 1, 2, 3, 4])
#plt.show()