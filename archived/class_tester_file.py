## Testing Flowline class
## Loading Jakobshavn mainline and testing CV optimisation of tau_y
## 5 Feb 2018  EHU

## 7 Mar 2018 revision: testing Network functions with Helheim branches

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

#-------------------
## READING IN BED
## COMMENT OUT IF DATA IS ALREADY READ IN TO YOUR SESSION
#-------------------

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

#JAK_0 = {}
#JAK_1 = {}
#JAK_2 = {}
#JAK_0['line'], JAK_1['line'], JAK_2['line'] = Flowline_CSV('Documents/1. Research/2. Flowline networks/Model/Greenland_tests/Flowline_sets/jakobshavn-thead_lines.csv', 3)
#
#jakobshavn = Flowline(JAK_0['line'], name='Jakobshavn', index=0)

#KB['line'] = Flowline_CSV('Documents/1. Research/2. Flowline networks/Model/Greenland_tests/Flowline_sets/koge_bugt-thead_lines.csv', 1)[0]

helcoords_0, helcoords_1, helcoords_2 = Flowline_CSV('Documents/1. Research/2. Flowline networks/Model/Greenland_tests/Flowline_sets/helheim-branch_lines.csv', 3)
hel_0 = Branch(coords=helcoords_0, index=0, order=0)
hel_1 = Branch(coords=helcoords_1, index=1, order=1, flows_to=0)
hel_2 = Branch(coords=helcoords_2, index=2, order=1, flows_to=0)
hel_branches = (hel_0, hel_1, hel_2)

helheim = PlasticNetwork(name='Helheim', init_type='Branch', branches=hel_branches, main_terminus=helcoords_0[0])
#hel_lines = TopToTerm(hel_branches) #OLD FORM
helheim.make_full_lines() #NEW FORM

#hel_master = [{'line': hel_lines[linenum]} for linenum in hel_lines]  #OLD FORM
#
#dictlist = (hel_master)
#keylist = ('bed', 'surface', 'thickness')
#fieldlist = (B_interp, S_interp, H_interp)
#for d in dictlist:
#    for j,k in enumerate(keylist):
#        d[k] = FlowProcess(d['line'], fieldlist[j])
#bestfit_hel = [(330000.0, 270000.0), (150000.0, 65000.0), (265000.0, 205000.0)]
helheim.process_full_lines(bed_field=B_interp, surface_field=S_interp, thickness_field=H_interp) #NEW FORM

#helheim.optimize_network_yield(testrange=arange(50e3, 500e3, 5e3), check_all=False)
helheim.flowlines[0].yield_type = 'constant'
helheim.flowlines[0].optimal_tau = 330e3
helheim.network_yield_type = helheim.flowlines[0].yield_type
helheim.network_tau = helheim.flowlines[0].optimal_tau

helheim.network_ref_profiles()


### Test model Helheim with optimal value
#for j,d in enumerate(dictlist): 
#    tau_0 = bestfit_hel[j][1]
#    arcmax = ArcArray(d['line'])[-1]
#    modelprof = PlasticProfile(d['bed'], tau_0, B_var, 0, d['surface'](0)/H0, arcmax, 10000, d['surface'])
#    modelint = interpolate.interp1d(modelprof[0], modelprof[1], kind='linear', copy=True)
#    d['Modelled'] = modelprof
#    d['Ref-profile-func'] = modelint
#    #d['Best tau_y'] = tau_y
#    d['Best tau_0'] = tau_0
    
#for j,d in enumerate(dictlist): 
#    tau_0 = bestfit_kb[j][1]
#    arcmax = ArcArray(d['line'])[-1]
#    modelprof = PlasticProfile(d['bed'], tau_0, B_var, 0, d['surface'](0)/H0, arcmax, 10000, d['surface'])
#    modelint = interpolate.interp1d(modelprof[0], modelprof[1], kind='linear', copy=True)
#    d['Modelled'] = modelprof
#    d['Ref-profile-func'] = modelint
#    #d['Best tau_y'] = tau_y
#    d['Best tau_0'] = tau_0

### Func to model forward from GL_model_tools
#testyears = arange(100)
#helheim_modeldicts = PlasticEvol(dictlist, testyears) 
#KB_modeldict = PlasticEvol(dictlist, testyears)[0]

### Plot results
#for j,d in enumerate(dictlist):
#    profile = d['Modelled']
#    taubest = d['Best tau_0']
#    sarr = profile[0]
#    amax = sarr[-1]
#    
#    plt.figure('Line '+str(j))
#    plt.plot(10*sarr, 1000*np.array(profile[2]), color='#9c4e16') #, label='Bed - Main branch'
#    plt.plot(10*sarr, 1000*np.array(profile[3]), color='Black', alpha=0.8, lw=2.5, label='GIMP DEM') #For reference
#    #plt.plot((10*ref07[0][0], 10*ref07[0][0]), (bedf(ref07[0][0]),1000*(ref07[1][0]-0.002)), color='Black', lw=2.5) #line for ice front
#    #plt.plot(10*sarr, 1000*np.array(profile[1]), '--', color='Grey', lw=3.0, label='Plastic, Ty = ' +str(tau_yield/1000)+' kPa')
#    plt.plot(10*sarr, 1000*np.array(profile[1]), '-.', color='DimGrey', lw=3.0, label = 'Plastic, Ty = ' +str(taubest/1000)+' + mu*N kPa')
#    plt.fill_between(10*sarr, 1000*np.array(profile[2]), (plt.axes().get_ylim()[0]), where=None, interpolate=True, color='#9c4e16', hatch='//', alpha=0.6) #flling bedrock...darker version of 'Chocolate' colour previously used
#    plt.fill_between(10*sarr, 1000*np.array(profile[1]), 1000*np.array(profile[2]), where=None, interpolate=True, color='Gainsboro') #filling reference profile
#    #plt.fill_between(10*xarr, 1000*np.array(jak_model_1[2]), zeros(len(xarr)), where=waterwhere, interpolate=True, color='MediumBlue', alpha=0.3) #filling water
#    #plt.plot((0, 10*term07), (0,0), color='MediumBlue', lw=2.0, label='Sea level') #Water/sea level indication
#    #plt.plot((15,15),(sef(1.5)+100, sef(1.5)+100), marker='v', ms=12.0, color='RoyalBlue') #marking thinning ref point
#    plt.xlabel('Distance along centerline (km)', fontsize=22)
#    plt.ylabel('Elevation (m.a.s.l.)', fontsize=22)
#    #plt.title('Jakobshavn Isbrae 2006-2014: Plastic glacier profiles vs. composite surface observations', fontsize=26)
#    #plt.title('Jakobshavn Isbrae 2006-2014: unsmoothed bed topography', fontsize=26)
#    plt.legend(loc='upper right', fontsize=20)
#    plt.axes().set_aspect(0.005)
#    plt.axes().set_xlim(right=0, left=10*amax)
#    #loc, labels = plt.yticks()
#    #plt.yticks(loc, np.array(['-1000', '', '0', '', '', '1500']))
#    plt.axes().tick_params(axis='y', direction='in', length=5, width=2, labelsize=24)
#    plt.axes().tick_params(axis='x', direction='in', length=5, width=2, labelsize=24)
#    plt.show()