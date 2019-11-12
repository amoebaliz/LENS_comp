import numpy as np
import netCDF4 as nc
import matplotlib
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

def plot_annual_maps():
    # Map specs
    m_lat = [20,50]; m_lon = [-140,-110]; m_off=3
    str_val = ['a','b']
    # Customize colormap
    disc_cmap = cmap_discretize(matplotlib.cm.nipy_spectral, end_mem)
    map_order=30

    # Rendering
    nplot=2
    fig = plt.figure(figsize=(11,5))
    gs = gridspec.GridSpec(1,nplot)
    gs.update(wspace=0.025)

    for n in range(nplot):
        ax = plt.subplot(gs[n])
        m = Basemap(llcrnrlat=m_lat[0]-m_off,urcrnrlat=m_lat[1]+m_off,\
            llcrnrlon=m_lon[0]-m_off,urcrnrlon=m_lon[1]+m_off, resolution='i', ax=ax)
        P = m.pcolormesh(lon,lat,dif_all_ind[-1*n,:].squeeze(),vmin=.5,vmax=end_mem+.5,cmap=disc_cmap)
        polygon_patch(m,ax)
        m.drawmeridians(m_lon, labels=[0,0,1,0], fmt='%d', fontsize=14)
        m.drawparallels(m_lat, labels=[0,0,0,0], fmt='%d')
        ax.text(m_lon[0]-2, m_lat[-1]+.8, str_val[n], fontsize=14)
        if n:
           print 'MEEP' 
           # plot colorbar axis
           bbox = ax.get_position()
           cax = fig.add_axes([bbox.xmax*1.0+.02, bbox.ymin, bbox.width*0.08, bbox.height-.085])
           cb = plt.colorbar(P,cax=cax,ticks=np.linspace(5,35,7))
           cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=12)
           cb.ax.set_title('LENS\nEndmember',horizontalalignment='center',\
                           multialignment='center')
           plt_str = str(nm).zfill(2) + '_MIN_DIF_SST_ENDMEMBER'
           plt.savefig(plt_str)
        else:
           m.drawparallels(m_lat, labels=[1,0,0,0], fmt='%d',fontsize=14)
        
def plot_clim_difs():
     # NOTE #5: Line plot
     print avg_dif_sst.shape
     # Plot each endmember against month
     #for nl in range(avg_dif_sst.shape[1]): 
     #    plt.plot(np.arange(1,12+1),avg_dif_sst[:,nl].squeeze())
     #    plt.show()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# NCAR directory for Large Ensemble files
dir = '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/'


# File name strings
his_fil_bas = 'b.e11.B20TRC5CNBDRD.f09_g16.'
fut_fil_bas = 'b.e11.BRCP85C5CNBDRD.f09_g16.'

fil_bas2 = '.pop.h.SST.'

# GULF OF AK DOMAIN
a = 314
b = 328
c = 206
d = 237

# WHETHER CALCULATING CLIMATOLOGY OR ANNUAL MEAN
# 0 = annual mean; 1 = monthly climatology
CLIM=1

if CLIM:
   nmons = 12
else:
   nmons = 1
print nmons

# Number of Years in climatology
nclim_yr = 10
nclim_yr = 30

# Number of NCAR LE endmembers
end_mem = 35

# Customize colormap
disc_cmap = cmap_discretize(matplotlib.cm.nipy_spectral, end_mem)
map_order=30

# Initializing difference matrix for line plot
avg_dif_sst=np.zeros((nmons,end_mem))

# EACH MONTH (if climatology)
for nm in range(nmons):
    print nm

    # STORAGE ARRAYS FOR MONTHLY CLIM MEAN/STD
    sst_all_stor = np.empty((0,b-a,d-c))
    dif_all_stor = np.empty((0,b-a,d-c))

    for ne in range(1,end_mem+1):
        ############################
        ## HISTORICAL CLIMATOLOGY ##
        ############################
        if ne == 1: 
           his_str_yr = 1850
        else:
           his_str_yr = 1920

        # File name
        his_ncfil  = dir + his_fil_bas + str(ne).zfill(3) + fil_bas2 + str(his_str_yr) + '01-' + str(200512)+ '.nc'
        fid_his = nc.Dataset(his_ncfil)

        # Montly indicies for historical fields
        nyr_h = 2005-his_str_yr+1
        Ih = np.arange(nm + 12*(nyr_h-nclim_yr), 12*nyr_h, nmons) 

        # Historical climatology
        SST_H = np.mean(fid_his.variables['SST'][Ih,:,a:b,c:d].squeeze(), axis=0)
        
        ######################## 
        ## FUTURE CLIMATOLOGY ##
        ########################
        fut_str_yr = [2006,2081]
        fut_end_yr = [2080,2100]

        if ne < 34:
           # Two files: 2006-2080, 2081-2100
           tot_yr_f = [nclim_yr-np.diff(fut_end_yr),np.diff(fut_end_yr)]
        else:
           # One file: 2006-2100
           tot_yr_f = [nclim_yr]

        # Opening and appending variable number of files
        SST_F = np.empty((0,b-a,d-c))
        nfils = len(tot_yr_f)
        for nf in range(nfils):
            fut_ncfil = dir + fut_fil_bas + str(ne).zfill(3) + fil_bas2 + str(fut_str_yr[nf]) + '01-' + str(fut_end_yr[nf-nfils])+ '12.nc'
            fid_fut = nc.Dataset(fut_ncfil)
            # Number of years contained in future file 
            nyr_f = fut_end_yr[nf]-fut_str_yr[nf-nfils]+1
            # (Montly) indicies for future fields
            If = np.arange(nm + int(12*(nyr_f-tot_yr_f[nf])), 12*nyr_f, nmons)
            # Monthly climatology/Annual Mean
            SST_F = np.ma.append(SST_F, fid_fut.variables['SST'][If,:,a:b,c:d].squeeze(), axis=0)
        # Climatology
        SST_F = np.mean(SST_F,axis=0)

        # STORE VALUES FOR ANALYSIS 
        sst_all_stor = np.ma.append(sst_all_stor, np.ma.expand_dims(SST_F,axis=0), axis=0)
        dif_all_stor = np.ma.append(dif_all_stor, np.ma.expand_dims(SST_F-SST_H,axis=0), axis=0)

    avg_dif_sst[nm,:] = np.ma.mean(np.ma.mean(dif_all_stor,axis=2),axis=1).squeeze()
    #print avg_dif_sst
    ##########
    # OUTPUT #
    ##########
    print "COLDEST"
    #print 'ABS', np.argsort(avg_fut_sst)[:5]+1
    print 'DIF', np.argsort(avg_dif_sst[nm,:].squeeze())[:5]+1

    print "WARMEST"
    #print 'ABS', np.argsort(avg_fut_sst)[-5:]+1
    print 'DIF', np.argsort(avg_dif_sst[nm,:].squeeze())[-5:]+1

    #print "MIDDLE"
    #print np.argsort(sst_stor)[(33/2)-2:(33/2)+3]+1   
    #print np.argsort(dif_stor)[(33/2)-2:(33/2)+3]+1

    ###########
    # FIGURES #
    ###########




#ax.plot(np.arange(1,33+1),sst_stor,'o')
#ax.errorbar(np.arange(1,33+1),sst_stor,yerr=sst_std,fmt='o')   
#ax.plot([1,33],sst_stor[6-1]*np.ones(2),'-b')
#ax.plot([1,33],sst_stor[16-1]*np.ones(2),'-r')
#ax.plot([1,33],sst_stor[33-1]*np.ones(2),'-k')
#plt.plot(SST)
    
#ax.errorbar(np.arange(1,33+1),sst_stor,yerr=sst_std,fmt='o')
#ax.plot(np.arange(1,33+1),sst_stor,'o')
plt.show()
