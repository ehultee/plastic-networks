## Scatter plot of year-by-year dL/dt observed versus simulated
## 26 Nov 2019  EHU

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy import interpolate
from scipy.ndimage import gaussian_filter
#from matplotlib.colors import LogNorm
from matplotlib import cm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
## Special import for SERMeQ modules
import sys
sys.path.insert(0, '/Users/lizz/Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *

### Topography needed to remove floating points from saved coords
###
print('Reading in surface topography')
gl_bed_path ='/Users/lizz/Documents/GitHub/Data_unsynced/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
b_raw = fh.variables['bed'][:].copy() # bed topo
e_raw = fh.variables['errbed'][:].copy() # error in bed topo or ice thickness
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
E = e_raw[::2, ::2]
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
errors = (5, 17, 18, 19, 29, 51, 71, 92, 95, 97, 100, 101, 102, 106, 107, 108, 109, 110, 113, 115, 117, 120, 121, 134, 168, 171) #glacier IDs that crashed in hindcasting 12 Mar 2019 *or* showed network problems 21 May 2019
## which ones get special treatment
added_jan19 = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)
seaward_projected = (61, 64, 82, 83, 99, 130, 132, 139, 140, 141, 156, 157, 158, 161, 167, 170, 178, 179, 180, 184) 
special_treatment = np.concatenate((added_jan19, seaward_projected))
rmv = np.concatenate((not_present, errors))
for n in rmv:
    try:
        glacier_ids.remove(n)
    except ValueError:
        pass
glaciers_to_plot=np.copy(glacier_ids).tolist()
for m in special_treatment:
    try:
        glaciers_to_plot.remove(m)
    except ValueError:
        pass

testyears = np.arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence',)
tempmarker = 'min10Cice' #and temperature of ice
timestepmarker = '8a_dt025a' #and total time and timestep

full_output_dicts = {}

for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for gid in glacier_ids:
        if gid in seaward_projected:
            fn = glob.glob('/Users/lizz/Documents/GitHub/Data_unsynced/Hindcasted_networks/advance_test/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] 
        else:
            fn = glob.glob('/Users/lizz/Documents/GitHub/Data_unsynced/Hindcasted_networks/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] #using glob * to select files of different run dates
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

gl_termpos_fldr = '/Users/lizz/Documents/GitHub/Data_unsynced/MEaSUREs-termini'
basefiles = ['/termini_0607_v01_2', '/termini_0708_v01_2', '/termini_0809_v01_2', '/termini_1213_v01_2', '/termini_1415_v01_2']
obs_years = [2006.75, 2007.75, 2009.0, 2013.0, 2014.75] #compare with years of hindcast, 2006-2014

termini = {}
for i,b in enumerate(basefiles):
    yr = obs_years[i]
    fn = gl_termpos_fldr+b
    termini[yr] = read_termini(fn, yr) #creating dictionary for each year
    print len(termini[yr])

nw_base_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Auto_selected-networks/Gld-autonetwork-GID'
seaward_coords_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Auto_selected-networks/Seaward_coords/Gld-advnetwork-GID' 
projected_termini = {gid: [] for gid in glacier_ids}
termpos_corrections = {gid: 0 for gid in glacier_ids}

for gid in glacier_ids: 
    print 'Reading in glacier ID: '+str(gid)
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
plot_years = 2006+np.array(testyears)
ids = [i for i in range(len(plot_years)) if plot_years[i] in obs_years] #which terminus positions will we compare
yr_colors = cm.get_cmap('Blues')(linspace(0.2, 0.9, num=len(ids)))


def make_error_boxes(ax, xdata, ydata, xerror, yerror, colorscheme_indices,
                     cmap='viridis', edgecolor='None', barcolor='None', alpha=0.5):
    """Make a collection of even-width lines (not rectangles) and add it to axes.
    Keyword args same as for mpl.patches.Rectangle, except 'colorscheme_indices'
    colorscheme_indices: array of same length as xdata, ydata setting assignment of facecolors"""
    
    # Plot horizontal lines
    for x, y, xe, ye, c in zip(xdata, ydata, xerror.T, yerror.T, colorscheme_indices):
        ax.hlines(y=y, xmin=x-xe[0], xmax=x+xe[1], lw=3.0, color=cm.get_cmap(cmap)((c-min(colorscheme_indices))/(max(colorscheme_indices)-min(colorscheme_indices))), alpha=alpha)
    artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                          fmt='None', ecolor=barcolor)
    return artists



# Create figure of 'un-futzed' glaciers
fig, ax = plt.subplots(1)
# Call function to create error boxes
for gid in glaciers_to_plot:
    tc = -1000*termpos_corrections[gid]
    sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    # sim_termini = np.zeros_like(obs_term_centr) # test constant 0 terminus
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    _ = make_error_boxes(ax, -1*obs_term_centr, -0.001*(tc + np.array(sim_termini)), xerror=e, yerror=0.1*np.ones(shape(e)), colorscheme_indices=obs_years)
ax.plot(range(-20,20), range(-20,20), c='k', linestyle='-.')
ax.axhline(y=0, linestyle='--', color='k')
ax.set_aspect(1)
ax.set_ylabel('Simulated $x_{term}$ [km]', fontsize=14)
ax.set_xlabel('Observed $x_{term}$ [km]', fontsize=14)
ax.set_title('Glaciers with usual processing')
plt.show()

## Same figure, 'special treatment' glaciers
fig2, ax2 = plt.subplots(1)
# Call function to create error boxes
for gid in np.setdiff1d(seaward_projected, added_jan19):
    tc = -1000*termpos_corrections[gid]
    sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    # sim_termini = np.zeros_like(obs_term_centr) # test constant 0 terminus
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    _ = make_error_boxes(ax2, -1*obs_term_centr, -0.001*(tc + np.array(sim_termini)), xerror=e, yerror=0.1*np.ones(shape(e)), colorscheme_indices=obs_years, alpha=0.8, cmap='Purples')
for gid in np.setdiff1d(added_jan19, seaward_projected): #sketchy flowline selection, but not projected forward
    if gid not in rmv:
        tc = -1000*termpos_corrections[gid]
        sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
        obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
        obs_term_centr = obs_termini[:,1]
        # sim_termini = np.zeros_like(obs_term_centr) # test constant 0 terminus
        e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
        _ = make_error_boxes(ax3, -1*obs_term_centr, -0.001*(tc + np.array(sim_termini)), xerror=e, yerror=0.1*np.ones(shape(e)), colorscheme_indices=obs_years, alpha=0.8, cmap='Reds')
ax2.plot(range(-20,20), range(-20,20), c='k', linestyle='-.')
ax2.set_aspect(1)
ax2.set_ylabel('Simulated $x_{term}$ [km]', fontsize=14)
ax2.set_xlabel('Observed $x_{term}$ [km]', fontsize=14)
ax2.set_title('Glaciers requiring special handling')
plt.show()

## All glaciers together, with futzed glaciers separated by colour
fig3, ax3 = plt.subplots(1)
x1 = range(-20, 20)
y1 = x1
ax3.plot(x1, y1, c='k', linestyle='-.')
ax3.fill_between(x1, y1=np.concatenate((range(-20,0), np.zeros(20))), y2=ax3.get_ylim()[0], color='DarkSlateGrey', alpha=0.5)
# Call function to create error boxes
for gid in glaciers_to_plot:
    tc = -1000*termpos_corrections[gid]
    sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    # sim_termini = np.zeros_like(obs_term_centr) # test constant 0 terminus
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    _ = make_error_boxes(ax3, -1*obs_term_centr, -0.001*(tc + np.array(sim_termini)), xerror=e, yerror=0.1*np.ones(shape(e)), colorscheme_indices=obs_years, alpha=0.8, cmap='Blues')
for gid in np.setdiff1d(seaward_projected, added_jan19): #projected forward, but otherwise trustworthy
    tc = -1000*termpos_corrections[gid]
    sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    # sim_termini = np.zeros_like(obs_term_centr) # test constant 0 terminus
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    _ = make_error_boxes(ax3, -1*obs_term_centr, -0.001*(tc + np.array(sim_termini)), xerror=e, yerror=0.1*np.ones(shape(e)), colorscheme_indices=obs_years, alpha=0.8, cmap='Purples')
for gid in np.setdiff1d(added_jan19, seaward_projected): #sketchy flowline selection, but not projected forward
    if gid not in rmv:
        tc = -1000*termpos_corrections[gid]
        sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
        obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
        obs_term_centr = obs_termini[:,1]
        # sim_termini = np.zeros_like(obs_term_centr) # test constant 0 terminus
        e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
        _ = make_error_boxes(ax3, -1*obs_term_centr, -0.001*(tc + np.array(sim_termini)), xerror=e, yerror=0.1*np.ones(shape(e)), colorscheme_indices=obs_years, alpha=0.8, cmap='Blues')
ax3.axhline(y=0, linestyle='--', color='Grey')
ax3.set_yscale('symlog')
ax3.set_xscale('symlog')
# ax3.set_aspect(1)
ax3.set_ylabel('Simulated $x_{term}$ [km]', fontsize=14)
ax3.set_xlabel('Observed $x_{term}$ [km]', fontsize=14)
ax3.set_xlim((-10, 15))
ax3.set_ylim((-20, 5))
ax3.set_title('All glaciers simulated')
plt.show()



## 1:1 sim-obs comparison with percent error and unit circle
def abs_percent_error(obs_list, sim_list):
    paired = zip(obs_list, sim_list)
    corrected = [(o, s) for (o, s) in paired if (abs(o)>0.1 and not np.isnan(o) and not np.isnan(s))] # avoid divide-by-zero errors
    perc_err = [100*abs((o-s)/(o)) for (o,s) in corrected]
    return perc_err

def unit_circle_compare(sim_pt, obs_pt):
    norm = np.sqrt(sim_pt**2 + obs_pt**2)
    return (sim_pt/norm, obs_pt/norm)

def range_normalized_centroid_diff(projected_obs, sim_termpos):
    """Calculate the difference between simulated terminus position and projected centroid of observed 2D terminus.
    Normalize by range of projected 2D terminus.
    Input:
        projected_obs: an array of shape (len(obs_years), 3) with an entry (advanced_extreme, centroid, retreated_extreme) for each year
        sim_termpos: an array of shape (len(obs_years)), with a simulated terminus position for each year"""
    joined = zip(projected_obs[:,0], projected_obs[:,2], projected_obs[:,1], sim_termpos)
    corrected = [(a, r, c, s) for (a,r,c,s) in joined if (abs(r-a)>0.1 and sum(np.isnan((a,r,c, s)))==0)]
    normalized_diff = [abs(c-s)/abs(r-a) for (a, r, c, s) in corrected]
    return normalized_diff        


# annual_pe_by_glacier = {gid: [] for gid in glaciers_to_plot}
# avg_pe_by_glacier = {gid:[] for gid in glaciers_to_plot}
# annual_uc_comp_by_glacier = {gid:[] for gid in glaciers_to_plot}
# avg_uc_comp_by_glacier = {gid:[] for gid in glaciers_to_plot}
# annual_rne_by_glacier = {gid: [] for gid in glaciers_to_plot}
# for gid in glaciers_to_plot:
#     sim_term_raw = [float(t) for t in full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini']]
#     sim_termini = -0.001*np.take(sim_term_raw, indices=ids)
#     obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
#     obs_term_centr = obs_termini[:,1]
#     sim_rates = diff(sim_termini)/diff(obs_years)
#     obs_rates = diff(obs_term_centr)/diff(obs_years)
#     annual_pe_by_glacier[gid] = abs_percent_error(obs_rates, sim_rates)
#     if len(annual_pe_by_glacier[gid])==0:
#         avg_pe_by_glacier[gid] = np.nan
#     else:
#         avg_pe_by_glacier[gid] = np.nanmean(annual_pe_by_glacier[gid])
#     annual_uc_comp_by_glacier[gid] = [unit_circle_compare(sim_rates[i], obs_rates[i]) for i in range(len(obs_rates))]
#     avg_uc_comp_by_glacier[gid] = unit_circle_compare(mean(sim_rates), mean(obs_rates))
#     annual_rne_by_glacier[gid] = range_normalized_centroid_diff(obs_termini, sim_termini)

# ## unit circle plot
# fig4 = plt.figure()
# ax4 = plt.axes([0,0,1,1])
# for gid in glaciers_to_plot:
#     uc_comp = np.asarray(annual_uc_comp_by_glacier[gid])
#     for j in range(len(uc_comp)):
#         ax4.plot((0, uc_comp[j,0]), (0, uc_comp[j,1]), color='LightGrey', lw=2.0)
# for gid in glaciers_to_plot: #plot average over top
#     avg_uc_comp = avg_uc_comp_by_glacier[gid]
#     ax4.plot((0, avg_uc_comp[0]), (0, avg_uc_comp[1]), color='DarkSlateGrey', lw=3.0, alpha=0.8)
# circ = plt.Circle((0, 0), radius=1., edgecolor='k', facecolor='None', lw=3.0, zorder=3)
# ax4.add_patch(circ)
# ax4.set_aspect(1)
# ax4.axhline(y=0, linestyle='--', color='k')
# ax4.axvline(x=0, linestyle='--', color='k')
# plt.axis('off')
# plt.margins(x=0.1, y=0.1)
# plt.show()

# ## absolute percent error plot
# pe_all = [annual_pe_by_glacier[gid] for gid in glaciers_to_plot]
# pe_arr = np.concatenate(pe_all)
# pe_weights = np.ones_like(pe_arr)/float(len(pe_arr))
# avg_pe = np.ravel([avg_pe_by_glacier[gid] for gid in glaciers_to_plot if not np.isnan(avg_pe_by_glacier[gid])])
# avg_pe_weights = np.ones_like(avg_pe)/float(len(avg_pe))
# pe_bins = np.linspace(0, 1000, num=20) 
# plt.figure('Absolute percent error observed - simulated annual retreat, Greenland outlets 2006-2014')
# plt.hist([pe_arr, avg_pe], bins=pe_bins, weights=[pe_weights, avg_pe_weights], color=['LightGrey', 'DarkSlateGrey'], label=['all', 'mean'])
# #plt.hist(avg_pe, bins=pe_bins, weights=avg_pe_weights, color='DarkViolet', alpha=0.5) # plot period-averaged percent diff
# plt.axes().tick_params(axis='both', length=5, width=2, labelsize=16)
# plt.xlabel('Percent difference $dL/dt$', fontsize=18)
# plt.ylabel('Normalized frequency', fontsize=18)
# plt.legend(loc='best')
# plt.axes().set_yticks([0, 0.1, 0.2])
# plt.axes().set_xlim((0,1000))
# plt.show()

# plt.figure('Cumulative distribution observed-simulated annual retreat, Greenland outlets 2006-2014')
# plt.hist([pe_arr, avg_pe], bins=pe_bins, weights=[pe_weights, avg_pe_weights], color=['LightGrey', 'DarkSlateGrey'], label=['all', 'mean'], cumulative=True)
# plt.axes().tick_params(axis='both', length=5, width=2, labelsize=16)
# plt.xlabel('Percent difference $dL/dt$ [m/a]', fontsize=18)
# plt.ylabel('Cumulative normalized frequency', fontsize=18)
# plt.legend(loc='best')
# plt.axes().set_xlim((0,1000))
# plt.show()


annual_rne_by_glacier = {gid: [] for gid in glaciers_to_plot}
for gid in glaciers_to_plot:
    sim_term_raw = [float(t) for t in full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini']]
    sim_termini = -0.001*np.take(sim_term_raw, indices=ids)    
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    annual_rne_by_glacier[gid] = range_normalized_centroid_diff(obs_termini, sim_termini)


## range-normalized error plot
rne_all = [annual_rne_by_glacier[gid] for gid in glaciers_to_plot]
rne_arr = np.concatenate(rne_all)
rne_weights = np.ones_like(rne_arr)/float(len(rne_arr))
rne_bins = np.linspace(0, 20, num=20) 
fig, ax = plt.subplots(1)
# plt.figure('Range normalized difference observed - simulated annual terminus position, Greenland outlets 2006-2014')
#plt.fill_between(x=(0, 1), y1=0.25, y2=0, color='DarkGreen', alpha=0.5)
#plt.fill_between(x=(1, 20), y1=0.25, y2=0, color='LightGreen', alpha=0.5)
ax.hist(rne_arr, bins=rne_bins, weights=rne_weights, color='DarkSlateGrey')
#plt.hist(avg_pe, bins=pe_bins, weights=avg_pe_weights, color='DarkViolet', alpha=0.5) # plot period-averaged percent diff
ax.tick_params(axis='both', length=5, width=2, labelsize=16)
ax.set_xlabel(r'Range-normalized difference $\frac{c_{obs}-c_{sim}}{\text{max}_{obs}-\text{min}_{obs}}$', fontsize=18)
ax.set_ylabel('Normalized frequency', fontsize=18)
ax.legend(loc='best')
ax.set(yticks=[0, 0.1, 0.2, 0.3, 0.4], xlim=(0,20))
plt.show()

## range-normalized retreat bound plot
mxrt_all = [annual_mxrt_by_glacier[gid] for gid in glaciers_to_plot]
mxrt_arr_cnct = np.concatenate(mxrt_all)
mxrt_arr = [m for m in mxrt_arr_cnct if np.isfinite(m)]
mxrt_weights = np.ones_like(mxrt_arr)/float(len(mxrt_arr))
mxrt_bins = np.linspace(0, 20, num=40) 
plt.figure('Normalized histo, range normalized difference annual minimum terminus position, Greenland outlets 2006-2014')
plt.hist(mxrt_arr, bins=mxrt_bins, weights=mxrt_weights, color='DarkSlateGrey')
#plt.hist(avg_pe, bins=pe_bins, weights=avg_pe_weights, color='DarkViolet', alpha=0.5) # plot period-averaged percent diff
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=16)
plt.xlabel('Range-normalized difference $(min_{obs}-c_{sim})/(max_{obs}-min_{obs})$', fontsize=18)
plt.ylabel('Normalized frequency', fontsize=18)
plt.legend(loc='best')
#plt.axes().set_yticks([0, 0.1, 0.2])
plt.axes().set_xlim((0,20))
plt.show()
