import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0,zorder=map_order+6)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False)
    axs.add_collection(lc)

def plot_map(ax):
    # Map specs
    m = Basemap(llcrnrlat=lat[0,0]-m_off,urcrnrlat=lat[-1,0]+m_off,\
           llcrnrlon=lon[0,0]-m_off,urcrnrlon=lon[0,-1]+m_off, resolution='i')

    m.drawmeridians([-160,-140], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
    m.drawparallels([55,60], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
    m.pcolormesh(pc_lon,pc_lat,pc_sst,cmap='seismic',zorder=map_order)
    polygon_patch(m,ax)
    ax.text(-162.5,61, str_val[n], fontsize=14)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# FIGURE

m_off = 0.1
map_order=1
str_val = ['a','b','c','d','e','f']

nplot = 6

fig = plt.figure(figsize=(11,4))
gs = gridspec.GridSpec(3,2)
gs.update(wspace=0.05, hspace=0.05)

# Iterate over different plots

for n in range(nplot):
    fid = nc.Dataset('AKG_test.nc')
    if n==0:

       lat = fid.variables['TLAT'][:]
       lon = fid.variables['TLONG'][:]
       lon[lon>180]=lon[lon>180]-360

       # Bilinear interpolation of lat/lon for pcolor corners
       pc_lat = (lat[:-1,:-1] + lat[1:,:-1] + lat[1:,1:] + lat[:-1,1:])/4
       pc_lon = (lon[:-1,:-1] + lon[1:,:-1] + lon[1:,1:] + lon[:-1,1:])/4

    sst = fid.variables['SST'][0,:].squeeze()

    ###########
    # MASKING #
    ###########

    for j in range(2,10):
        for i in range(j):
            sst[j,i-1] = np.ma.masked

    pc_sst = sst[1:-1,1:-1]

    ax = plt.subplot(gs[n])
    plot_map(ax)

plt.show()

