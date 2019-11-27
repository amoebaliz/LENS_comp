import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.dates as pltd
import datetime as dt

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0,zorder=map_order+6)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False)
    axs.add_collection(lc)

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_map(ax):
    # PLOT: Map object
    m = Basemap(llcrnrlat=lat[0,0]-m_off,urcrnrlat=lat[-1,0]+m_off,\
           llcrnrlon=lon[0,0]-m_off,urcrnrlon=lon[0,-1]+m_off, resolution='i')

    # PLOT: SST pcolormesh
    m.pcolormesh(pc_lon,pc_lat,pc_sst,vmin=3.0,vmax=5.5,cmap=new_cmap,zorder=map_order)
    polygon_patch(m,ax)

    # PLOT: Subplot ID
    ax.text(-162.5,tx_lat, ids[n], fontsize=14)
    # PLOT: Endmember ID
    enm_str = 'EM ' + dat_fils[n][:3] 
    ax.text(-137, tx_lat, enm_str, fontsize=12) 

    # PLOT: Latitudes (left plots)
    m.drawparallels([55,60], labels=[int((n+1)%2),0,0,0], fmt='%d', fontsize=12,zorder=map_order+5)
    # PLOT: Longitudes (bottom plots)
    m.drawmeridians([-160,-140], labels=[0,0,0,int(n/4)], fmt='%d', fontsize=12,zorder=map_order+5)

    # PLOT: Titles (top plots)
    if n<2:
       ax.set_title(titles[n], fontsize=14)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# FIGURE

m_off = 0.1
map_order=1
ids = ['a','b','c','d','e','f']
titles = ['10 years', '30 years']
nplot = 6
tx_lat = 61

fig = plt.figure(figsize=(11,6))
gs = gridspec.GridSpec(3,2)
gs.update(wspace=0.05, hspace=0.05)

# SETTING UP COLORMAP
cmap = plt.get_cmap('hot_r')
new_cmap = truncate_colormap(cmap, 0.2, 1)

# LENS VALUES TO PLOT
dat_fils = ['013_average_10yr_dif', '013_average_30yr_dif',\
            'ALL_average_10yr_dif', 'ALL_average_30yr_dif',\
            '031_average_10yr_dif', '031_average_30yr_dif' ] 

fid = nc.Dataset('AKG_test.nc')
lat = fid.variables['TLAT'][:]
lon = fid.variables['TLONG'][:]
lon[lon>180]=lon[lon>180]-360

# Bilinear interpolation of lat/lon for pcolor corners
pc_lat = (lat[:-1,:-1] + lat[1:,:-1] + lat[1:,1:] + lat[:-1,1:])/4
pc_lon = (lon[:-1,:-1] + lon[1:,:-1] + lon[1:,1:] + lon[:-1,1:])/4

for n in range(nplot):
    # OPEN: SST file
    sst = np.load(dat_fils[n]).squeeze()

    # MASK: Bering Sea
    for j in range(2,11):
        for i in range(j-1):
            sst[j,i] = np.ma.masked

    # SUBSET: For pcolor plot
    pc_sst = sst[1:-1,1:-1]

    print ''
    print 'PLOT', ids[n]
    print np.min(pc_sst),np.max(pc_sst)

    # PLOT: Maps
    ax = plt.subplot(gs[n])
    plot_map(ax)

# ADJUST: plot positions for colorbar and timeline
plt.subplots_adjust(right=.88, bottom=0.26)
fig.canvas.draw_idle()

# PLOT: Colorbar
cb_x0 = fig.axes[nplot-1].get_position().x1 + \
       (fig.axes[nplot-1].get_position().x0-fig.axes[0].get_position().x1)
cb_y0 = fig.axes[nplot-1].get_position().y0
cb_yn = fig.axes[0].get_position().y1-fig.axes[nplot-1].get_position().y0
cax = fig.add_axes([cb_x0, cb_y0, 0.02, cb_yn])

plt.colorbar(ax.collections[0],cax=cax,extend='both')  
cax.get_yaxis().labelpad = 25
cax.set_ylabel(r'$\Delta$ SST ($\!^\circ\!$C)', size=14,rotation=270)
cax.tick_params(labelsize=11) 

# PLOT: Timeline
tm_x0 = fig.axes[0].get_position().x0

tm_xn = fig.axes[nplot-1].get_position().x1 - \
        fig.axes[0].get_position().x0

ax2 = plt.axes([tm_x0,.1,tm_xn,.05], facecolor=(1,1,1,0)) 
#ax2.set_title('Differential (Future-Historical) Timeframes',fontsize=14)
ax2.set_title('Timeperiods Averaged for Differencing (Future-Historical)',fontsize=12)

tim_lin = []
for ny in np.arange(1970,2111,10):
    tim_lin.append(pltd.date2num(dt.datetime(ny,1,1)))

# PLOT: Full timeline
ax2.plot([tim_lin[0],tim_lin[-1]],[0,0],'k')

# PLOT: 30-yr timelines
ax2.plot([pltd.date2num(dt.datetime(1971,1,1)),\
          pltd.date2num(dt.datetime(2000,12,31))],\
          [.2,.2],'-b',linewidth=4)

ax2.plot([pltd.date2num(dt.datetime(2071,1,1)),\
          pltd.date2num(dt.datetime(2100,12,31))],\
          [.2,.2],'-b',linewidth=4)

# PLOT: 10-yr timelines
ax2.plot([pltd.date2num(dt.datetime(1981,1,1)),\
          pltd.date2num(dt.datetime(1990,12,31))],\
          [.4,.4],'-g',linewidth=4)

ax2.plot([pltd.date2num(dt.datetime(2081,1,1)),\
          pltd.date2num(dt.datetime(2090,12,31))],\
          [.4,.4],'-g',linewidth=4)

ax2.set_xticks(tim_lin)
ax2.tick_params(axis="x", labelsize=11, labelrotation=40)
ax2.xaxis_date()
ax2.set_ylim(0,.7)
ax2.set_xlim(tim_lin[0],tim_lin[-1])
ax2.get_yaxis().set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

# SAVE: Figure
plt.savefig('GoAK_aliasing_example')
plt.show()

