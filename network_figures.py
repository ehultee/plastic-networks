#Flowline network figures
#7 May 2018  EHU

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
Jakobshavn_main.load_network(filename='JakobshavnIsbrae-main_south.pickle', load_mainline_output=False)

jakcoords_sec = Flowline_CSV('Jakobshavn_secondary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_1 = Flowline(coords=jakcoords_sec, index=1, name='Jak secondary branch', has_width=True)
Jakobshavn_sec = PlasticNetwork(name='Jakobshavn Isbrae [secondary/central]', init_type='Flowline', branches=(jak_1), main_terminus=jakcoords_sec[0])
Jakobshavn_sec.load_network(filename='Jakobshavn_sec.pickle', load_mainline_output=False)

jakcoords_tert = Flowline_CSV('Jakobshavn_tertiary-flowline-w_width.csv', 1, has_width=True, flip_order=False)[0]
jak_2 = Flowline(coords=jakcoords_tert, index=2, name='Jak tertiary branch', has_width=True)
Jakobshavn_tert = PlasticNetwork(name='Jakobshavn Isbrae [tertiary/north]', init_type='Flowline', branches=(jak_2), main_terminus=jakcoords_tert[0])
Jakobshavn_tert.load_network(filename='Jakobshavn_tert.pickle', load_mainline_output=False)

kbcoords = Flowline_CSV('KogeBugt-mainline-w_width.csv', 1, has_width=True, flip_order=True)[0]
kb_line = Flowline(coords=kbcoords, index=0, name='Koge Bugt mainline', has_width=True)
KogeBugt = PlasticNetwork(name='Koge Bugt', init_type='Flowline', branches=(kb_line), main_terminus=kbcoords[0])
KogeBugt.load_network(filename='KogeBugt.pickle', load_mainline_output=False)

#
#### INTERSECTING LINES
helcoords_0, helcoords_1, helcoords_2 = Flowline_CSV('Helheim-network-w_width.csv', 3, has_width=True, flip_order=False)
hel_0 = Branch(coords=helcoords_0, index=0, order=0)
hel_1 = Branch(coords=helcoords_1, index=1, order=1, flows_to=0)
hel_2 = Branch(coords=helcoords_2, index=2, order=1, flows_to=0)
hel_branches = (hel_0, hel_1, hel_2)
Helheim = PlasticNetwork(name='Helheim', init_type='Branch', branches=hel_branches, main_terminus=helcoords_0[0])
Helheim.make_full_lines()
Helheim.load_network(filename='Helheim.pickle', load_mainline_output=False)

kangercoords_0, kangercoords_1, kangercoords_2, kangercoords_3, kangercoords_4 = Flowline_CSV('Documents/GitHub/plastic-networks/kangerlussuaq-network-w_width.csv', 5, has_width=True, flip_order=False)
kanger_0 = Branch(coords=kangercoords_0, index=0, order=0)
kanger_1 = Branch(coords=kangercoords_1, index=1, order=1, flows_to=0, intersect=174)
kanger_2 = Branch(coords=kangercoords_2, index=2, order=1, flows_to=0, intersect=191) #DIFFERENT FROM PREVIOUS BRANCH 2.  NEW FLOWLINE SET AS OF 31 MAR 2018
kanger_3 = Branch(coords=kangercoords_3, index=3, order=1, flows_to=0, intersect=146)
kanger_4 = Branch(coords=kangercoords_4, index=4, order=1, flows_to=0, intersect=61)
kanger_branches = (kanger_0, kanger_1, kanger_3, kanger_4)
Kanger = PlasticNetwork(name='Kangerlussuaq', init_type='Branch', branches=kanger_branches, main_terminus=kangercoords_0[0])
Kanger.make_full_lines()
Kanger.load_network(filename='Kangerlussuaq.pickle', load_mainline_output=False)


##-------------------
### PLOTTING FUNCTIONS
##-------------------

glacier_networks = (Jakobshavn_main, Jakobshavn_sec, Jakobshavn_tert, KogeBugt, Helheim, Kanger) #list which glaciers we're handling

for gl in glacier_networks:
    xcoords = np.array([fl.coords[:,0] for fl in gl.flowlines])
    ycoords = np.array([fl.coords[:,1] for fl in gl.flowlines])
    
    plt.figure(gl.name)
    plt.contourf(X, Y, S, cmap='Purples_r')
    for fl in gl.flowlines:
        plt.plot(fl.coords[:,0], fl.coords[:,1], color='k', linewidth=5, marker='.')
    fl0 = gl.flowlines[0]
    plt.plot(fl0.coords[::25,0], fl0.coords[::25,1], color='w', marker='o')
    plt.axes().set_aspect(1)
    plt.axes().set_xlim(left = min(flatten(xcoords))-10000, right = max(flatten(xcoords))+10000)
    plt.axes().set_ylim(bottom = min(flatten(ycoords))-10000, top = max(flatten(ycoords))+10000)
    plt.show()

#Jakobshavn plot
Jak_lines = (Jakobshavn_main.flowlines[0], Jakobshavn_sec.flowlines[0], Jakobshavn_tert.flowlines[0])
Jak_x = np.array([fl.coords[:,0] for fl in Jak_lines])
Jak_y = np.array([fl.coords[:,1] for fl in Jak_lines])
plt.figure('Jakobshavn combined')
plt.contourf(X, Y, S, cmap='Purples_r')
for line in Jak_lines:
    plt.plot(line.coords[:,0], line.coords[:,1], color='k', linewidth=5, marker='.')
fl0=Jak_lines[0]
plt.plot(fl0.coords[::25,0], fl0.coords[::25,1], color='w', marker='o')
plt.axes().set_aspect(1)
plt.axes().set_xlim(left = min(flatten(Jak_x))-10000, right = max(flatten(Jak_x))+10000)
plt.axes().set_ylim(bottom = min(flatten(Jak_y))-10000, top = max(flatten(Jak_y))+10000)
plt.show()