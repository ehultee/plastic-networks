## Creating maps of flowline networks projected on satellite imagery
## 12 Jun 2018  EHU

from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np
import matplotlib.pyplot as plt

 
#selection of which map appearance you want, projection settings
services = ['World_Physical_Map', 'World_Shaded_Relief', 'World_Topo_Map', 'NatGeo_World_Map', 'ESRI_Imagery_World_2D', 'World_Street_Map', 'World_Imagery', 'Ocean_Basemap']
wgs84 = pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
psn_gl = pyproj.Proj("+init=epsg:3413") # Polar Stereographic North used by BedMachine 



def Greenland_map(service='ESRI_Imagery_World_2D', epsg=3413, xpixels=5000):
    """Function using Basemap to plot map for all of Greenland
    """
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=300, llcrnrlat=57, urcrnrlon=20, urcrnrlat=80, resolution='h')
    
    plt.figure()
    m.arcgisimage(service=service, xpixels=xpixels)
    plt.show()
    return m

def SermeqKujalleq_map(service='ESRI_Imagery_World_2D', epsg=3413, xpixels=2000):
    """Function using Basemap to plot map for only Sermeq Kujalleq
    """
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=309.6, llcrnrlat=68.9, urcrnrlon=314.6, urcrnrlat=69.6, resolution='h')
    
    plt.figure()
    m.arcgisimage(service=service, xpixels=xpixels)
    plt.show()
    return m

def Helheim_map(service='ESRI_Imagery_World_2D', epsg=3413, xpixels=2000):
    """Function using Basemap to plot map for only Helheim Glacier
    """
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=320.3, llcrnrlat=66.1, urcrnrlon=322.4, urcrnrlat=67.6, resolution='h')
    
    plt.figure()
    m.arcgisimage(service=service, xpixels=xpixels)
    plt.show()
    return m

def Kanger_map(service='ESRI_Imagery_World_2D', epsg=3413, xpixels=2000):
    """Function using Basemap to plot map for only Kangerlussuaq
    """
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=325.4, llcrnrlat=68.4, urcrnrlon=328.2, urcrnrlat=69.2, resolution='h')
    
    plt.figure()
    m.arcgisimage(service=service, xpixels=xpixels)
    plt.show()
    return m

def KogeBugt_map(service='ESRI_Imagery_World_2D', epsg=3413, xpixels=2000):
    """Function using Basemap to plot map for only Kangerlussuaq
    """
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=318.2, llcrnrlat=64.9, urcrnrlon=319.4, urcrnrlat=65.4, resolution='h')
    
    plt.figure()
    m.arcgisimage(service=service, xpixels=xpixels)
    plt.show()
    return m    
    
##Convert coords into lat/lon so that Basemap can convert them back (don't know why this is necessary, but it works)
def flowline_latlon(coords, fromproj=pyproj.Proj("+init=epsg:3413"), toproj=pyproj.Proj("+init=EPSG:4326")):
    """Convert coords into lat/lon so that Basemap can convert them back for plotting (don't know why this is necessary, but it works)
    Defaults:
        fromproj = NSIDC Polar Stereographic North
        toproj = WGS84 lat-lon
    """
    xs = coords[:,0]
    ys = coords[:,1]
    x_lon, y_lat = pyproj.transform(fromproj, toproj, xs, ys)
    latlon_coords = np.asarray(zip(x_lon, y_lat))
    return latlon_coords

glacier_networks = (Jakobshavn_main, Jakobshavn_sec, Jakobshavn_tert, KogeBugt, Helheim, Kanger) #all glaciers defined in Greenland-summary_plotting
all_coords_latlon = []
for gl in glacier_networks:
    for fl in gl.flowlines:
        latlon_coords = flowline_latlon(fl.coords)
        all_coords_latlon.append(latlon_coords)
all_coords_latlon = np.asarray(all_coords_latlon)        

### Reading in SENTINEL velocity map
print 'Now reading in (vector) velocity map'
v_path = 'Documents/1. Research/2. Flowline networks/Model/Data/ESA-Greenland/greenland_iv_500m_s1_20161223_20170227_v1_0.nc'
fh2 = Dataset(v_path, mode='r')
xv = fh2.variables['x'][:].copy()
yv = fh2.variables['y'][:].copy()
#yv = yv_flipped[::-1]
v_raw = fh2.variables['land_ice_surface_velocity_magnitude'][:].copy() #this is v(y, x)
vx_raw = fh2.variables['land_ice_surface_easting_velocity'][:].copy()
vy_raw =fh2.variables['land_ice_surface_northing_velocity'][:].copy()
v_upper = np.ma.masked_greater(v_raw, 10000)
vx_upper = np.ma.masked_greater(vx_raw, 10000)
vy_upper = np.ma.masked_greater(vy_raw, 10000)
fh2.close()

yv_flipped = yv[::-1]
xvm, yvm = np.meshgrid(xv, yv_flipped)
#x_lon, y_lat = pyproj.transform(psn_gl, wgs84, xvm, yvm)
#regridded_v = interpolate.griddata((xvm.ravel(), yvm.ravel()), v_upper.ravel(), (x_lon, y_lat), method='nearest') #do we need to regrid, or just call contour?

#Plot flowlines (or just termini) on map of Greenland
grnld = Greenland_map()
x_lon, y_lat = grnld(xvm, yvm, inverse=True)
grnld.contourf(x_lon, y_lat, v_upper, latlon=True, cmap='viridis', norm=matplotlib.colors.LogNorm(vmin=0.01, vmax=20.0), levels=(0.01, 0.5, 1.0, 5.0, 10.0, 50.0), zorder=3, alpha=0.5) #overlay Sentinel velocity
#for j in range(len(all_coords_latlon)):
#    grnld.scatter(all_coords_latlon[j][:,0], all_coords_latlon[j][:,1], s=5, color='DarkSlateBlue', latlon=True) #flowlines
#    grnld.scatter(all_coords_latlon[j][0,0], all_coords_latlon[j][0,1], s=100, marker='*', color='Cyan', edgecolors='DarkSlateBlue', linewidths=1, latlon=True) #termini
#grnld.plot(Jak_main_lon, Jak_main_lat, latlon=True)
plt.show()

plt.figure()
sk = SermeqKujalleq_map()
for j in (0, 1, 2):
    sk.scatter(all_coords_latlon[j][:,0], all_coords_latlon[j][:,1], color='DarkSlateBlue', latlon=True) #flowlines
plt.show()

plt.figure()
kb = KogeBugt_map()
kb.scatter(all_coords_latlon[3][:,0], all_coords_latlon[3][:,1], color='DarkSlateBlue', latlon=True)
kb.scatter(all_coords_latlon[3][0,0], all_coords_latlon[3][0,1], s=150, marker='*', color='Cyan', edgecolors='DarkSlateBlue', linewidths=1, latlon=True)
plt.show()

plt.figure()
hh = Helheim_map()
for j in (4, 5, 6):
    hh.scatter(all_coords_latlon[j][:,0], all_coords_latlon[j][:,1], color='DarkSlateBlue', latlon=True)
hh.scatter(all_coords_latlon[4][0,0], all_coords_latlon[4][0,1], s=150, marker='*', color='Cyan', edgecolors='DarkSlateBlue', linewidths=1, latlon=True)
plt.show()

plt.figure()
k = Kanger_map()
for j in (7,8,9,10):
    k.scatter(all_coords_latlon[j][:,0], all_coords_latlon[j][:,1], color='DarkSlateBlue', latlon=True)
k.scatter(all_coords_latlon[7][0,0], all_coords_latlon[7][0,1], s=150, marker='*', color='Cyan', edgecolors='DarkSlateBlue', linewidths=1, latlon=True)
plt.show()



## Contour sentinel velocities separately for eventual overlay
plt.figure(frameon=False)
vs = plt.contourf(xvm, yvm, v_upper, latlon=False, cmap='viridis', norm=matplotlib.colors.LogNorm(), levels=(0.001, 0.01, 0.1, 1.0, 10.0, 50.0, 100.0))
cbar=plt.colorbar(vs)
cbar.set_ticks([0.01, 0.1, 1.0, 10.0, 50.0, 100.0])
cbar.set_ticklabels(['0.01', '0.1', '1.0', '10.0', '50.0', '100.0'])
plt.clim(0.01, 50.0)
plt.axes().set_yticks([])
plt.axes().set_xticks([])
plt.axes().set_aspect(1)
plt.show()
