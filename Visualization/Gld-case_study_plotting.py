## Plotting individual glacier results for case studies
## 24 Apr 2019  EHU

import numpy as np
import matplotlib.pyplot as plt
import csv
import shapefile
#import collections
import glob
#from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib.patches as mpatches
from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
## Special import for SERMeQ modules
import sys
sys.path.insert(0, 'Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *


### Topography needed to remove floating points from saved coords
###
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
#bb = np.ma.masked_where(thick_mask !=2, b_raw)
bb = b_raw #don't mask, to allow bed sampling from modern bathymetry (was subglacial in ~2006)
## Down-sampling
X = xx[::2]
Y = yy[::2]
S = ss[::2, ::2]
H = hh[::2, ::2] 
B = bb[::2, ::2]
M = thick_mask[::2,::2]
## Not down-sampling
#X = xx
#Y = yy
#S = ss
fh.close()

#Smoothing bed and surface
unsmoothB = B
smoothB = gaussian_filter(B, 2)
smoothS = gaussian_filter(S, 2)
#B_processed = np.ma.masked_where(thick_mask !=2, smoothB)

#Replacing interpolated surface with bed+thickness
S_new = np.add(B, H)

S_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothS.T[::, ::-1])
H_interp = interpolate.RectBivariateSpline(X, Y[::-1], H.T[::, ::-1])
B_interp = interpolate.RectBivariateSpline(X, Y[::-1], smoothB.T[::, ::-1])



###--------------------------------------
#### GLACIERS TO PLOT
###--------------------------------------
## Which glaciers are available
glacier_ids = range(1,195) #MEaSUREs glacier IDs to process.
not_present = (93, 94, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 169) #glacier IDs missing from set
## which ones get special treatment
added_jan19 = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)
seaward_projected = (61, 64, 82, 83, 99, 130, 132, 139, 140, 141, 156, 157, 158, 161, 167, 170, 178, 179, 180, 184) 
special_treatment = np.concatenate((added_jan19, seaward_projected))
errors = (5, 18, 19, 29, 71, 92, 95, 97, 101, 107, 108, 120, 134) #glacier IDs that crashed in hindcasting 12 Mar 2019
rmv = np.concatenate((not_present, errors))
for n in rmv:
    try:
        glacier_ids.remove(n)
    except ValueError:
        pass
#glaciers_to_plot=np.copy(glacier_ids).tolist()
#for m in special_treatment:
#    try:
#        glaciers_to_plot.remove(m)
#    except ValueError:
#        pass

glaciers_to_plot = [g for g in glacier_ids if g in (3, 105, 137, 175)]


testyears = arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence',)

#datemarker = '2019-02-08' #markers on filename to indicate date run
tempmarker = 'min10Cice' #and temperature of ice
timestepmarker = '8a_dt025a' #and total time and timestep

full_output_dicts = {}

for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for gid in glaciers_to_plot:
        fn = glob.glob('Documents/GitHub/Data_unsynced/Hindcasted_networks/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] #using glob * to select files of multiple run dates
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    #for i, gid in enumerate(glaciers_to_plot):
    #    fn = 'GID{}-{}-{}-{}-{}.pickle'.format(gid, datemarker, s, tempmarker, timestepmarker)
    #    lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

perscenario_fluxes = []
perscenario_SLE = []

for s in full_output_dicts.keys():
    print 'Scenario {}'.format(s)
    out = full_output_dicts[s]
    pernetwork_cumul_fx = []
    pernetwork_cumul_sle = []
    for j, gid in enumerate(glaciers_to_plot):
        branch_fx = [np.nan_to_num(out['GID{}'.format(gid)][k]['Terminus_flux']) for k in range(len(out['GID{}'.format(gid)]))]
        total_fx = sum(branch_fx, axis=0)
        total_sle = (1E-12)*np.array(total_fx)/(361.8) #Gt ice/mm SLE conversion
        cumul_fx = np.cumsum(total_fx)
        cumul_sle = np.cumsum(total_sle)
        pernetwork_cumul_fx.append(cumul_fx)
        pernetwork_cumul_sle.append(cumul_sle)
    scenario_flux = np.cumsum(pernetwork_cumul_fx, axis=0)
    perscenario_fluxes.append(scenario_flux[-1])
    scenario_sle = np.cumsum(pernetwork_cumul_sle, axis=0)
    print max(scenario_sle[-1])
    print(scenario_sle[-1][-1])
    perscenario_SLE.append(scenario_sle[-1])

gl_termpos_fldr = 'Documents/GitHub/Data_unsynced/MEaSUREs-termini'
basefiles = ['/termini_0607_v01_2', '/termini_0708_v01_2', '/termini_0809_v01_2', '/termini_1213_v01_2', '/termini_1415_v01_2', '/termini_1516_v01_2']
obs_years = [2006, 2007, 2008, 2012, 2014, 2015] #compare with term of hindcast, 2006-2014

termini = {}
for i,b in enumerate(basefiles):
    yr = obs_years[i]
    fn = gl_termpos_fldr+b
    termini[yr] = read_termini(fn, yr) #creating dictionary for each year
    print len(termini[yr])

nw_base_fpath = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/Gld-autonetwork-GID'
projected_termini = {gid: [] for gid in glaciers_to_plot}
seaward_coords_fpath = 'Documents/GitHub/Data_unsynced/Auto_selected-networks/Seaward_coords/Gld-advnetwork-GID' 
termpos_corrections = {gid: 0 for gid in glacier_ids}


for gid in glaciers_to_plot:
    print 'Reading in glacier ID: '+str(gid)
    #if gid in added_jan19:
    #    filename = nw_base_fpath+str(gid)+'-date_2019-01-10.csv'
    filename = glob.glob(nw_base_fpath+'{}-date_*.csv'.format(gid))[0] #using glob * to select files of different run dates
    
    coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
    if gid in seaward_projected:
        seaward_fn = seaward_coords_fpath+'{}-fwd_2000_m.csv'.format(gid)
        seaward_coords = Flowline_CSV(seaward_fn, has_width=True, flip_order=True)[0]
        branch_0 = Branch(coords=np.concatenate((seaward_coords, coords_list[0])), index=0, order=0) #saving extended central branch as main
        termpos_correction = 10*max(ArcArray(seaward_coords)) #how much glacier length has been added to initial line, i.e. how much terminus shifted in coordinate system, in km
        print termpos_correction
    else:
        branch_0 = Branch(coords=coords_list[0], index=0, order=0) #saving central branch as main
        termpos_correction = 0
    termpos_corrections[gid] = termpos_correction
    branch_list = [branch_0]

    nw = PlasticNetwork(name='GID'+str(gid), init_type='Branch', branches=branch_list, main_terminus=branch_0.coords[0])
    nw.make_full_lines()
    if gid not in seaward_projected:  #remove floating, but not from lines that have been artificially extended
        print 'Removing floating points from glacier ID: '+str(gid)
        nw.process_full_lines(B_interp, S_interp, H_interp)
        nw.remove_floating()
    mainline = LineString(nw.flowlines[0].coords)
    
    for yr in obs_years:
        try:
            termpts = termini[yr][gid] #get terminus points for each year
            t = projected_term_obs(termpts, mainline) #project onto main flowline
            r = retterm(termpts, mainline) #find most retreated point
            a = advterm(termpts, mainline) #find most advanced point
            print 'GID {}, t={}, r={}, a={}'.format(gid, t, r, a)
            projected_termini[gid].append((-1*termpos_correction)+np.asarray((a, t, r))) #add these to dictionary of projected termini per glacier
        except KeyError:
            print 'No terminus found in year {} for GID {}.'.format(yr, gid)
            projected_termini[gid].append((0, np.nan, 0))    
    

###--------------------------------------
#### PLOTTING
###--------------------------------------   
  
## Settings for plots   
#labels = [str(g[1]) for g in gids_by_name if g[0] in glaciers_to_plot] #set what the glaciers will be called in plotting.  Default is simply their MEaSUREs ID
labels = [str(g) for g in glaciers_to_plot]
markers = ['o', '.', ',', '^', 'd', '*']
styles = ['-', ':', '-.', '-', '-', '-']
cmap = cm.get_cmap('winter')
scenario_colors = cm.get_cmap('Blues')([0.1, 0.3, 0.5, 0.7, 0.9])
#colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
colors = cmap(linspace(0.1, 0.9, num=len(glaciers_to_plot)))
alt_colors = cm.get_cmap('Greys')([0.2, 0.3, 0.5, 0.7, 0.9])
plot_years = 2006+np.array(testyears)


#
#####terminus
#plt.figure()
#for j, gid in enumerate(glaciers_to_plot):
#    print gid
#    term_positions = full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::]
#    ms_selection = mod(j, len(styles))
#    plt.plot(testyears, -0.001*np.array(term_positions), linewidth=2, color='Gainsboro', linestyle=styles[ms_selection], label=labels[j])
#    #plt.plot(testyears[::4], -0.001*np.array(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::])[::4], linewidth=0, marker=markers[ms_selection], ms=10, color=colors[j])
#plt.legend(loc='lower left')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Terminus change [km]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_ylim(-100, 1)
##plt.axes().set_yticks([-75, -50, -25, 0])
#plt.title('Terminus retreat of {} Greenland outlet glaciers 2006-2014 ERA-I, Tice=-10 C'.format(len(glaciers_to_plot)), fontsize=20)
#plt.show()

####SINGLE NETWORK - termini vs obs
for j, gid in enumerate(glaciers_to_plot):
    sim_termini = full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::]
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    plt.figure('Main line terminus change, GID{}'.format(gid))
    plt.plot(plot_years, -0.001*np.array(sim_termini), linewidth=2, color='k', linestyle='-', label='Modelled')
    plt.errorbar(obs_years, -1*obs_term_centr, yerr = e, fmt='D')
    plt.axes().tick_params(axis='both', length=5, width=2, labelsize=30)
    plt.axes().set_xlim(2006, 2014.5)
    plt.axes().set_xticks([2006, 2008, 2010, 2012, 2014])
    plt.axes().set_ylim(-20, 2)
    plt.axes().set_yticks([-20, -15, -10, -5, 0])
    if gid==175:
        plt.axes().set_xticklabels(['2006', '2008', '2010', '2012', '2014'])
        plt.axes().set_xlabel('Year', size=35)
    else:
        plt.axes().set_xticklabels([])
    if gid==3:
        plt.axes().set_yticklabels(['-20', '', '-10', '', '0'])
        plt.axes().set_ylabel('Terminus change [km]', size=35)
    else:
        plt.axes().set_yticklabels([])
    plt.axes().set_aspect(0.3)
    plt.show()


##
#####SINGLE NETWORK - splitting termini
#single_network_output = full_output_dicts['persistence']['GID10']
#plt.figure('Single network terminus change')
#for k in range(len(single_network_output)): #for each branch j
#    #colork = alt_colors[mod(k, len(colors))]
#    branch_termini = single_network_output[k]['Termini'] #for some reason enumerate doesn't work with loaded-in output, so we're stuck with this
#    #markerk = (k+2, mod(k+1, 3), 0)
#    plt.plot(testyears, -0.001*np.array(branch_termini[0:-1:]), linewidth=4, label='Branch {}'.format(k))
#    #plt.plot(testyears[::10], -0.001*np.array(branch_termini[0:-1:])[::10], linewidth=0, color=colork, marker=markerk, ms=10)
#plt.legend(loc='lower left')
#plt.axes().set_xlabel('Year of simulation', size=30)
#plt.axes().set_ylabel('Terminus change [km]', size=30)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_ylim(-50, 1)
##plt.axes().set_yticks([-50, -25, 0])
#plt.show()
#
#####Flux
#for j, gid in enumerate(glaciers_to_plot):
#    plt.figure('GID{}'.format(gid))
#    ms_selection = mod(gid, len(styles))
#    plt.plot(testyears, 1E-12*np.array(pernetwork_cumul_fx[j]), linewidth=4, label=str(gid), color=colors[ms_selection], linestyle = styles[ms_selection])
#    plt.plot(testyears[::20], 1E-12*np.array(pernetwork_cumul_fx[j][::20]), linewidth=0, ms=10, marker=markers[ms_selection])
#    plt.fill_between(testyears, y1=1E-12*np.array(pernetwork_cumul_fx[j]), y2=0, alpha=0.5)
#    plt.axes().set_xlabel('Year of simulation', size=20)
#    plt.axes().set_ylabel('Cumulative ice flux [Gt]', size=20)
#    plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#    plt.axes().set_xlim(0,50)
#    plt.axes().set_xticks([0,10,20, 30, 40, 50])
#    #plt.axes().set_ylim(0, 401)
#    #plt.axes().set_yticks([0, 50, 100, 150, 200, 250, 300, 350, 400])
#    plt.show()
    ##
#######Sea level equivalent
#plt.figure(figsize=(12,8))
#for j, gid in enumerate(glaciers_to_plot):
#    #if gid!=10:
#    ms_selection = mod(gid, len(styles))
#    plt.plot(testyears[::], scenario_sle[j], linewidth=4, color=colors[ms_selection], label=gid)
#    plt.plot(testyears[::5], scenario_sle[j][::5], linewidth=0, marker=markers[ms_selection], ms=10, color=colors[ms_selection])
#    if j==0:
#        plt.fill_between(testyears[::], y1=scenario_sle[j], y2=0, color=colors[ms_selection], alpha=0.7)  
#    else:
#        plt.fill_between(testyears[::], y1=scenario_sle[j], y2=scenario_sle[j-1], color=colors[ms_selection], alpha=0.7)     
##plt.plot([0, 20, 40, 60], [0, 14, 28, 42], color='k', linewidth=1, ls='-', alpha=0.8) #GRACE linear trend
##rect = mpatches.Rectangle((98,8.5), width=2, height=4.6, color='k', alpha=0.7) # Nick et al total projection for 2100, A1B
##rect2 = mpatches.Rectangle((98,11.3), width=2, height=6.2, color='k', alpha=0.7) # Nick et al total projection for 2100, RCP8.5 
##plt.axes().add_patch(rect)
##plt.axes().add_patch(rect2)
#plt.legend(loc='upper left')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_xlim(0, 100)
##plt.axes().set_xticks([0, 25, 50, 75, 100])
##plt.axes().set_ylim(0, 12)
##plt.axes().set_yticks([0, 2, 4, 6, 8, 10, 12])
#plt.show()