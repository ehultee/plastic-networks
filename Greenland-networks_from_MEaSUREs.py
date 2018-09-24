## Finding flowline networks for all ~200 MEaSUREs outlet glaciers
## 24 Sept 2018  EHU
import numpy as np
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter #use Savitzky-Golay filter to smooth catchments in Trace_wWidth


## Read in MEaSUREs velocity composite
print 'Reading MEaSUREs velocities'
x_comp, y_comp, v_comp_raw = read_velocities('Documents/GitHub/gld-velocity-composite.tif')
vx_comp_raw = read_velocities('Documents/GitHub/gld-x_velocity-composite.tif', return_grid=False)
vy_comp_raw = read_velocities('Documents/GitHub/gld-y_velocity-composite.tif', return_grid=False)
v_comp = np.ma.masked_invalid(v_comp_raw)
vx_comp = np.ma.masked_invalid(vx_comp_raw)
vy_comp = np.ma.masked_invalid(vy_comp_raw)
#v_excludemasked = np.ma.filled(v_comp, fill_value=np.nanmean(v_comp))
#vx_excludemasked = np.ma.filled(vx_comp, fill_value=np.nanmean(vx_comp))
#vy_excludemasked = np.ma.filled(vy_comp, fill_value=np.nanmean(vy_comp))


##Make 2D-interpolated function of velocity field for tracing
print 'Interpolating MEaSUREs velocity composites for tracing'
x_1d = x_comp[0,:]
y_1d = y_comp[:,0]
#v_interp = interpolate.RectBivariateSpline(x_1d, y_1d[::-1], v_excludemasked.T[::, ::-1])
func_vxcomp = interpolate.interp2d(x_1d, y_1d, vx_comp[::,::]) #check whether these need to be flipped along y-axis
func_vycomp = interpolate.interp2d(x_1d, y_1d, vy_comp[::,::])
func_vcomp = interpolate.interp2d(x_1d, y_1d[::-1], v_comp[::, ::])

xtest = np.linspace(min(x_1d), max(x_1d), 200)
ytest = np.linspace(min(y_1d), max(y_1d), 200)
plt.figure()
try:
    plt.contour(xtest, ytest, func_vcomp(xtest, ytest))
except ValueError: #catching nonsense warning
    pass
plt.show()

## Using read_termini from MEaSUREs-validation.py
gl_termpos_fldr = 'Documents/GitHub/plastic-networks/Data/MEaSUREs-termini'
terminus_basefile = '/termini_0607_v01_2'
init_year = 2006
fn = gl_termpos_fldr + terminus_basefile #filename to read in for termini that will be traced
termini_init = read_termini(fn, init_year)

##iterate over keys in termini_init to make dictionary of lines for each GlacierID
#ids_to_trace = termini_init.keys() #trace all points of all glaciers
ids_to_trace = (3, 153, 175) # IDs for only Jakobshavn, Kangerlussuaq, Helheim

all_lines = {}
for gid in ids_to_trace:
    lines = {}
    termcoords = termini_init[gid] #points spanning terminus for this glacier
    for j in range(len(termcoords)):
        p = termcoords[j]
        line_coords, width = Trace_wWidth(p[0], p[1], trace_up=True, xarr=x_comp, yarr=y_comp, Vx = func_vxcomp, Vy = func_vycomp, V = func_vcomp) #Uses Trace_wWidth and FilterMainTributaries from network_selection.py
        xyw = [(line_coords[n][0], line_coords[n][1], width[n]) for n in range(len(line_coords))]
        lines[j] = (xyw)
    filtered_tribs = FilterMainTributaries(lines, Vx = func_vxcomp, Vy = func_vycomp)
    all_lines[gid] = lines

        