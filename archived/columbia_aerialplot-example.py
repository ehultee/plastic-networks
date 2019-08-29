# Columbia retreat 1957-2007, viewed from above, contour plot by colour
# Jun/Nov 2016  EHU

#from numpy import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from plastic_utilities import *
from StringIO import StringIO

# Define Dimensional and Dimensionless parameters for model use
H0=1e3 #characteristic height for nondimensionalisation 
L0=10e3 #characteristic length (10km)
tau_yield = 150e3 
tau_0 = 130e3 #for effective pressure format
g = 9.8 
rho_ice = 920.0 #ice density kg/m^3
rho_sea=1020.0 #seawater density kg/m^3

#Bingham stress functions-----------------
def B_var(elev, thick, pos, time): #variation by altitude and ice thickness (effective pressure at base)...pos, time arguments required by plasticmodel
    if elev<0:
        D = -elev #Water depth D the nondim bed topography value when Z<0
    else:
        D = 0
    N = rho_ice*g*thick*H0 - rho_sea*g*D*H0
    mu = 0.01 #Coefficient between 0 and 1
    tau_y = tau_0 + mu*N
    return tau_y/(rho_ice*g*H0**2/L0) 
    
def B_const(elev, thick, pos, time):  #functional form of B if using a constant yield strength
    return tau_yield/(rho_ice*g*H0**2/L0)
#-----------------------------------------

#Reading in original Columbia dataset
columbiadata ='Data/McNabb_11J249_S1.nc'
fh = Dataset(columbiadata, mode='r')
X=fh.variables['X'][:].copy() # Get grid of Easting coordinates of points
Y=fh.variables['Y'][:].copy() # Same for Northing
H57=fh.variables['H57'][:].copy() #1957 ice thickness
H07=fh.variables['H07'][:].copy()
Z=fh.variables['Z'][:].copy() # Bed topography
xx = X[:,0] #X has same value in each row - value only changes by column, so flattening by taking first dimension of X
yy = Y[0,:]  #Y flattening in same way

unsmoothZ = Z
Z = gaussian_filter(Z,3)

#Set up 2d interpolation--this is not affected by choice of flowline
xx = X[:,0] #X has same value in each row - value only changes by column, so flattening by taking first dimension of X
yy = Y[0,:]  #Y flattening in same way
bedinterp = interpolate.interp2d(xx,yy,Z.T, kind='linear', copy=True, bounds_error=True) #need to interpolate transpose of Z because of how data is stored...unclear why this works this way
thickinterp = interpolate.interp2d(xx,yy,H57.T, kind='linear', copy=True, bounds_error=True)
seinterp = interpolate.interp2d(xx,yy,(H57+unsmoothZ).T, kind='linear', copy=True, bounds_error=True)
thick07interp = interpolate.interp2d(xx, yy, H07.T, kind='linear', copy=True, bounds_error=True)
se07interp = interpolate.interp2d(xx,yy,(H07+unsmoothZ).T, kind='linear', copy=True, bounds_error=True)


mainline = flowlinearr('Data/mainline.txt')
mainarr = evenspace(mainline, 250)
bedf = flowprocess(mainarr, bedinterp)
thickf = flowprocess(mainarr, thickinterp)
sef = flowprocess(mainarr, seinterp)
thick07f = flowprocess(mainarr, thick07interp)
se07f = flowprocess(mainarr, se07interp)

westline = flowlinearr('Data/westtribline.txt')
westarr = evenspace(westline, 250)
bedf_2 = flowprocess(westarr, bedinterp)
thickf_2 = flowprocess(westarr, thickinterp)
sef_2 = flowprocess(westarr, seinterp)
thick07f_2 = flowprocess(westarr, thick07interp)
se07f_2 = flowprocess(westarr, se07interp)

eastline = flowlinearr('Data/easttribline.txt')
eastarr = evenspace(eastline, 250)
bedf_3 = flowprocess(eastarr, bedinterp)
thickf_3 = flowprocess(eastarr, thickinterp)
sef_3 = flowprocess(eastarr, seinterp)
thick07f_3 = flowprocess(eastarr, thick07interp)
se07f_3 = flowprocess(eastarr, se07interp)


#arcmax1 = 6.453912514944193 #nondim max arclength of main flowline, retrieved from "arcarray(mainarr)[-1]"
arcmax1 = arcarray(mainarr)[-1]
arcmax2 = arcarray(westarr)[-1]
arcmax3 = arcarray(eastarr)[-1]

term07 = 1.61620408461 #2007 position where thickness makes water balance (tested with Ty=150kpa)

#print 'Run: 1957, main branch, Ty =' + str(0.001*tau_yield) + 'kPa'
#model57_1 = plasticmodel(bedf, B_const, 0, balancethick((bedf(0)/H0), B_const(bedf(0), thickf(0), 0, 0))+(bedf(0)/H0), arcmax1, 25000, sef) 
print 'Run: 1957, main branch, Ty=' + str(0.001*tau_0) + 'kPa +mu*N'
model57var_1 = plasticmodel(bedf, B_var, 0, balancethick((bedf(0)/H0), B_var(bedf(0)/H0, thickf(0)/H0, 0, 0))+(bedf(0)/H0), arcmax1, 25000, sef)
modelint_57_var_1 = interpolate.interp1d(model57var_1[0], model57var_1[1], kind='linear', copy=True)
#print '\n Run: 2007, main branch, Ty =' + str(0.001*tau_yield) + 'kPa'
#model07_1 = plasticmodel(bedf, B_const, term07, balancethick((bedf(term07)/H0), B_const(bedf(term07)/H0, thick07f(term07)/H0, 0, 0))+(bedf(term07)/H0), arcmax1, 25000, se07f)
print '\n Run: 2007, main branch, Ty=' + str(0.001*tau_0) + 'kPa +mu*N'
model07var_1 = plasticmodel(bedf, B_var, term07, balancethick((bedf(term07)/H0), B_var(bedf(term07)/H0, thick07f(term07)/H0, 0, 0))+(bedf(term07)/H0), arcmax1, 25000, se07f)
modelint_07_var_1 = interpolate.interp1d(model07var_1[0], model07var_1[1], kind='linear', copy=True)

#print '\n Run: 1957, west tributary, Ty =' + str(0.001*tau_yield) + 'kPa'
#model57_2 = plasticmodel(bedf_2, B_const, 0, balancethick((bedf_2(0)/H0), B_const(bedf_2(0)/H0, thickf_2(0)/H0, 0, 0))+(bedf_2(0)/H0), arcmax2, 25000, sef_2) 
print '\n Run: 1957, west tributary, Ty=' + str(0.001*tau_0) + 'kPa +mu*N'
model57var_2 = plasticmodel(bedf_2, B_var, 0, balancethick((bedf_2(0)/H0), B_var(bedf_2(0)/H0, thickf_2(0)/H0, 0, 0))+(bedf_2(0)/H0), arcmax2, 25000, sef_2)
modelint_57_var_2 = interpolate.interp1d(model57var_2[0], model57var_2[1], kind='linear', copy=True)
#print '\n Run: 2007, west tributary, Ty =' + str(0.001*tau_yield) + 'kPa'
#model07_2 = plasticmodel(bedf_2, B_const, term07, balancethick((bedf_2(term07)/H0), B_const(bedf_2(term07)/H0, thick07f_2(term07)/H0, 0, 0))+(bedf_2(term07)/H0), arcmax2, 25000, se07f_2)
print '\n Run: 2007, west tributary, Ty=' + str(0.001*tau_0) + 'kPa +mu*N'
model07var_2 = plasticmodel(bedf_2, B_var, term07, balancethick((bedf_2(term07)/H0), B_var(bedf_2(term07)/H0, thick07f_2(term07)/H0, 0, 0))+(bedf_2(term07)/H0), arcmax2, 25000, se07f_2)
modelint_07_var_2 = interpolate.interp1d(model07var_2[0], model07var_2[1], kind='linear', copy=True)

#print '\n Run: 1957, east tributary, Ty =' + str(0.001*tau_yield) + 'kPa'
#model57_3 = plasticmodel(bedf_3, B_const, 0, balancethick((bedf_3(0)/H0), B_const(bedf_3(0)/H0, thickf_3(0)/H0, 0, 0))+(bedf_3(0)/H0), arcmax3, 25000, sef_3) 
print '\n Run: 1957, east tributary, Ty=' + str(0.001*tau_0) + 'kPa +mu*N'
model57var_3 = plasticmodel(bedf_3, B_var, 0, balancethick((bedf_3(0)/H0), B_var(bedf_3(0)/H0, thickf_3(0)/H0, 0, 0))+(bedf_3(0)/H0), arcmax3, 25000, sef_3)
modelint_57_var_3 = interpolate.interp1d(model57var_3[0], model57var_3[1], kind='linear', copy=True)
#print '\n Run: 2007, east tributary, Ty =' + str(0.001*tau_yield) + 'kPa'
#model07_3 = plasticmodel(bedf_3, B_const, term07, balancethick((bedf_3(term07)/H0), B_const(bedf_3(term07)/H0, thick07f_3(term07)/H0, 0, 0))+(bedf_3(term07)/H0), arcmax3, 25000, se07f_3)
print '\n Run: 2007, east tributary, Ty=' + str(0.001*tau_0) + 'kPa +mu*N'
model07var_3 = plasticmodel(bedf_3, B_var, term07, balancethick((bedf_3(term07)/H0), B_var(bedf_3(term07)/H0, thick07f_3(term07)/H0, 0, 0))+(bedf_3(term07)/H0), arcmax3, 25000, se07f_3)
modelint_07_var_3 = interpolate.interp1d(model07var_3[0], model07var_3[1], kind='linear', copy=True)


## Plotting thickness along flowline-------------------------

all_lines_arr = np.concatenate((mainarr, westarr, eastarr))

#FOR 1957
modelsepts_main = modelint_57_var_1(arcarray(mainarr)[:-3]) #have to stop before end to avoid interpolation error
modelthickpts_main = np.array(modelint_57_var_1(arcarray(mainarr)[:-3]))-0.001*np.array(bedf(arcarray(mainarr)[:-3]))
modelsepts_west = modelint_57_var_2(arcarray(westarr))
modelthickpts_west = np.array(modelint_57_var_2(arcarray(westarr)))-0.001*np.array(bedf_2(arcarray(westarr)))
modelsepts_east = modelint_57_var_3(arcarray(eastarr))
modelthickpts_east = np.array(modelint_57_var_3(arcarray(eastarr)))-0.001*np.array(bedf_3(arcarray(eastarr)))

##Comparison with obs values
#obssepts_main = sef(arcarray(mainarr)[:-3]) #have to stop before end to avoid interpolation error
#obsthickpts_main = thickf(arcarray(mainarr)[:-3])
#obssepts_west = sef_2(arcarray(westarr)) #have to stop before end to avoid interpolation error
#obsthickpts_west = thickf_2(arcarray(westarr))
#obssepts_east = sef_3(arcarray(eastarr)) #have to stop before end to avoid interpolation error
#obsthickpts_east = thickf_3(arcarray(eastarr))


##arrays for plotting 2007
mainarr07 = []
for j, val in enumerate(arcarray(mainarr)):
    if val > term07:
        mainarr07.append(mainarr[j])
mainarr07 = np.array((mainarr07))
westarr07 = []
for j, val in enumerate(arcarray(westarr)):
    if val > term07:
        westarr07.append(westarr[j])
westarr07 = np.array(westarr07)
eastarr07 = []
for j, val in enumerate(arcarray(eastarr)):
    if val > term07:
        eastarr07.append(eastarr[j])
eastarr07 = np.array(eastarr07)

mainstart = next(j for j, val in enumerate(arcarray(mainarr)) if val>term07)
weststart = next(j for j, val in enumerate(arcarray(westarr)) if val>term07)
eaststart = next(j for j, val in enumerate(arcarray(eastarr)) if val>term07)

#2007 PLOT
#modelsepts_main_07 = modelint_07_var_1(arcarray(mainarr)[mainstart:])
modelthickpts_main_07 = np.array(modelint_07_var_1(arcarray(mainarr)[mainstart:]))-0.001*np.array(bedf(arcarray(mainarr)[mainstart:]))
#modelsepts_west_07 = modelint_07_var_2(arcarray(westarr)[weststart:])
modelthickpts_west_07 = np.array(modelint_07_var_2(arcarray(westarr)[weststart:]))-0.001*np.array(bedf_2(arcarray(westarr)[weststart:]))
#modelsepts_east_07 = modelint_07_var_3(arcarray(eastarr)[eaststart:])
modelthickpts_east_07 = np.array(modelint_07_var_3(arcarray(eastarr)[eaststart:]))-0.001*np.array(bedf_3(arcarray(eastarr)[eaststart:]))

#For plotting - comparisons
all_lines_thick = H0*np.array(np.concatenate((modelthickpts_main, modelthickpts_west, modelthickpts_east)))
all_lines_thick_07 = H0*np.array(np.concatenate((modelthickpts_main_07, modelthickpts_west_07, modelthickpts_east_07)))
#Lost ice arrays
lostpts_main = zeros(len(modelthickpts_main)-len(modelthickpts_main_07))
main_comparray_07 = np.concatenate((lostpts_main, modelthickpts_main_07))
main_diffarray = np.array(modelthickpts_main)-np.array(main_comparray_07)
lostpts_west = zeros(len(modelthickpts_west)-len(modelthickpts_west_07))
west_comparray_07 = np.concatenate((lostpts_west, modelthickpts_west_07))
west_diffarray = np.array(modelthickpts_west)-np.array(west_comparray_07)
lostpts_east = zeros(len(modelthickpts_east)-len(modelthickpts_east_07))
east_comparray_07 = np.concatenate((lostpts_east, modelthickpts_east_07))
east_diffarray = np.array(modelthickpts_east)-np.array(east_comparray_07)
#
all_lines_diff = -H0*np.array(np.concatenate((main_diffarray, west_diffarray, east_diffarray)))


#Composite figure
plt.figure()
ax1 = plt.subplot(131) #Model 1957
ax2 = plt.subplot(132) #Model 2007
ax3 = plt.subplot(133) #Difference plot


class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
#Making colour scheme based on 1957 to display changes
mainthick_57 = np.array(modelint_57_var_1(arcarray(mainarr)[:-3]))-0.001*np.array(bedf(arcarray(mainarr)[:-3]))
westthick_57 = np.array(modelint_57_var_2(arcarray(westarr)))-0.001*np.array(bedf_2(arcarray(westarr)))
eastthick_57 = np.array(modelint_57_var_3(arcarray(eastarr)))-0.001*np.array(bedf_3(arcarray(eastarr)))
thickmax = H0*max(mainthick_57.max(), westthick_57.max(), eastthick_57.max())

##UNCOMMENT BELOW FOR MAIN BRANCH PLOT
#sc=ax1.scatter(0.001*np.array(mainarr[:-3,0]), 0.001*np.array(mainarr[:-3,1]), c=all_lines_thick[:len(mainarr)-3], vmin=0, vmax=H0*mainthick_57.max(), s=80, cmap = cm.Blues, edgecolor='none')
##UNCOMMENT BELOW FOR FULL THREE BRANCH PLOT
sc=ax1.scatter(0.001*np.array(mainarr[:-3,0]), 0.001*np.array(mainarr[:-3,1]), c=all_lines_thick[:len(mainarr)-3], vmin=0, vmax=thickmax, s=80, cmap = cm.Blues, edgecolor='none')
ax1.scatter(0.001*np.array(westarr[:, 0]), 0.001*np.array(westarr[:,1]), c=all_lines_thick[len(mainarr)-3:len(mainarr)+len(westarr)-3], vmin=0, vmax=thickmax, cmap=cm.Blues, s=80, edgecolor='none')
ax1.scatter(0.001*np.array(eastarr[:, 0]), 0.001*np.array(eastarr[:,1]), c=all_lines_thick[len(mainarr)+len(westarr)-3:], vmin=0, vmax=thickmax, cmap=cm.Blues, s=80, edgecolor='none')
#ax1.title('Columbia Glacier network: ice thickness', fontsize=20)
#ax1.set_xlabel('Easting (km, UTM 6V)', fontsize = 20)
#ax1.ylabel('Northing (km, UTM 6V)', fontsize = 20)
ax1.set_aspect(1)
ax1.tick_params(axis='both', direction='in', length=5, width=2, labelsize=18)
ax1.set_xticks([480, 490, 500, 510, 520, 525])
ax1.set_xticklabels(['', '10 km'])
ax1.set_yticks([6765, 6775, 6785, 6795, 6805, 6815])
ax1.set_yticklabels([])
ax1.set_ylim(6760, 6815)
ax1.add_patch(Rectangle((480, 6760), 10, 2, facecolor='Black', edgecolor='Black', lw=3.0))
##append axes for colorbar
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", "5%", pad="3%")
##UNCOMMENT BELOW FOR MAIN BRANCH PLOT
#cbar1 = plt.colorbar(sc, ticks=[0, 300, 600, 900], cax=cax1)
#cbar1.ax.set_yticklabels(['0', '300', '600', '900'])
##UNCOMMENT BELOW FOR FULL THREE BRANCH PLOT
cbar1 = plt.colorbar(sc, ticks=[0, 250, 500, 750, 1000], cax=cax1)
cbar1.ax.set_yticklabels(['0', '250', '500', '750', '1000'])
cbar1.ax.tick_params(axis='y', direction='in', length=5, width=2, labelsize=18)


##UNCOMMENT BELOW FOR MAIN BRANCH PLOT
#sc2=ax2.scatter(0.001*np.array(mainarr07[:-3,0]), 0.001*np.array(mainarr07[:-3,1]), c=all_lines_thick_07[:len(mainarr07)-3], vmin=0, vmax=H0*mainthick_57.max(), s=80, cmap = cm.Blues, edgecolor='none')
##UNCOMMENT BELOW FOR FULL THREE BRANCH PLOT
sc2=ax2.scatter(0.001*np.array(mainarr07[:-3,0]), 0.001*np.array(mainarr07[:-3,1]), c=all_lines_thick_07[:len(mainarr07)-3], vmin=0, vmax=thickmax, s=80, cmap = cm.Blues, edgecolor='none')
ax2.scatter(0.001*np.array(westarr07[:, 0]), 0.001*np.array(westarr07[:,1]), c=all_lines_thick_07[len(mainarr07)-3:len(mainarr07)+len(westarr07)-3], vmin=0, vmax=thickmax, cmap=cm.Blues, s=80, edgecolor='none')
ax2.scatter(0.001*np.array(eastarr07[:, 0]), 0.001*np.array(eastarr07[:,1]), c=all_lines_thick_07[len(mainarr07)+len(westarr07)-3:-3], vmin=0, vmax=thickmax, cmap=cm.Blues, s=80, edgecolor='none')
#ax2.title('Columbia Glacier network: 2007 ice thickness', fontsize=20)
#ax2.set_xlabel('Easting (km, UTM 6V)', fontsize = 20)
#ax2.ylabel('Northing (km, UTM 6V)', fontsize = 20)
ax2.set_aspect(1)
ax2.tick_params(axis='both', direction='in', length=5, width=2, labelsize=18)
ax2.set_xticks([480, 490, 500, 510, 520, 525])
ax2.set_xticklabels(['', '10 km'])
ax2.set_yticks([6765, 6775, 6785, 6795, 6805, 6815])
ax2.set_yticklabels([])
ax2.set_ylim(6760, 6815)
ax2.add_patch(Rectangle((480, 6760), 10, 2, facecolor='Black', edgecolor='Black', lw=3.0))
#append axes for colorbar
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", "5%", pad="3%")
##UNCOMMENT BELOW FOR MAIN BRANCH PLOT
#cbar2 = plt.colorbar(sc2, ticks=[0, 300, 600, 900], cax=cax2)
#cbar2.ax.set_yticklabels(['0', '300', '600', '900'])
##UNCOMMENT BELOW FOR FULL THREE BRANCH PLOT
cbar2 = plt.colorbar(sc2, ticks=[0, 250, 500, 750, 1000], cax=cax2)
cbar2.ax.set_yticklabels(['0', '250', '500', '750', '1000'])
cbar2.ax.tick_params(axis='y', direction='in', length=5, width=2, labelsize=18)


sc3=ax3.scatter(0.001*np.array(mainarr[:-3,0]), 0.001*np.array(mainarr[:-3,1]), c=all_lines_diff[:len(mainarr)-3], norm=MidpointNormalize(vmin=min(all_lines_diff), vmax=max(all_lines_diff), midpoint=0.), s=80, cmap = cm.RdBu, edgecolor='none')
##UNCOMMENT BELOW FOR FULL THREE BRANCH PLOT
ax3.scatter(0.001*np.array(westarr[:, 0]), 0.001*np.array(westarr[:,1]), c=all_lines_diff[len(mainarr)-3:len(mainarr)+len(westarr)-3], norm=MidpointNormalize(vmin=min(all_lines_diff), vmax=max(all_lines_diff), midpoint=0.), cmap = cm.RdBu, s=80, edgecolor='none')
ax3.scatter(0.001*np.array(eastarr[:, 0]), 0.001*np.array(eastarr[:,1]), c=all_lines_diff[len(mainarr)+len(westarr)-3:], norm=MidpointNormalize(vmin=min(all_lines_diff), vmax=max(all_lines_diff), midpoint=0.), cmap = cm.RdBu, s=80, edgecolor='none')
#ax3.title('Columbia Glacier network: change in ice thickness 1957-2007', fontsize=20)
#ax3.set_xlabel('Easting (km, UTM 6V)', fontsize = 20)
#ax3.ylabel('Northing (km, UTM 6V)', fontsize = 20)
ax3.set_aspect(1)
ax3.tick_params(axis='both', direction='in', length=5, width=2, labelsize=18)
ax3.set_xticks([480, 490, 500, 510, 520, 525])
ax3.set_xticklabels(['', '10 km'])
ax3.set_yticks([6765, 6775, 6785, 6795, 6805, 6815])
ax3.set_yticklabels([])
ax3.set_ylim(6760, 6815)
ax3.add_patch(Rectangle((480, 6760), 10, 2, facecolor='Black', edgecolor='Black', lw=3.0))
#append axes for colorbar
divider = make_axes_locatable(ax3)
cax3 = divider.append_axes("right", "5%", pad="3%")
#cbar3 = plt.colorbar(sc3, cax=cax3)
cbar3 = plt.colorbar(sc3, ticks=[-900, -600, -300, 0, 150], cax=cax3)
cbar3.ax.set_yticklabels(['-900', '-600', '-300', '0', '+150']) #change in ice thickness 1957->2007, in m
cbar3.ax.tick_params(axis='y', direction='in', length=5, width=2, labelsize=18)

plt.show()