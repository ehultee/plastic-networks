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
    m = Basemap(projection='npstere', boundinglat=70, lon_0=315, epsg=epsg, llcrnrlon=320.3, llcrnrlat=66.1, urcrnrlon=322.4, urcrnrlat=67.5, resolution='h')
    
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



#Plot flowlines (or just termini) on map of Greenland
grnld = Greenland_map()
for j in range(len(all_coords_latlon)):
    grnld.scatter(all_coords_latlon[j][0,0], all_coords_latlon[j][0,1], marker='*', color='Cyan', linewidths=8, latlon=True) #termini
    #grnld.scatter(all_coords_latlon[j][:,0], all_coords_latlon[j][:,1], color='k', latlon=True) #flowlines
#grnld.plot(Jak_main_lon, Jak_main_lat, latlon=True)
plt.show()
