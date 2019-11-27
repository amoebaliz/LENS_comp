import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# NCAR directory for Large Ensemble files
dir = '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/SST/'

# File name strings
his_fil_bas = 'b.e11.B20TRC5CNBDRD.f09_g16.'
fut_fil_bas = 'b.e11.BRCP85C5CNBDRD.f09_g16.'

fil_bas2 = '.pop.h.SST.'

# GULF OF AK DOMAIN
# for ncks -d nlat,313,329 -d nlon,206,238
# in python: same start value, +1 last value
a = 313
b = 330
c = 206
d = 239
nmons=1

# Number of Years in climatology
nclim_yr = 30

# Number of NCAR LE endmembers
end_mem = 35

# DIFS ARE FOR 30-yr CLIM DIFS
#ne_vals = [31] # COLDEST
#ne_vals = [13] # WARMEST
#ne_vals = np.arange(35)+1 # ALL

ne_vals = [[31],[13],np.arange(35)+1]

#dif_all_30_stor = np.empty((0,b-a,d-c))
#dif_all_10_stor = np.empty((0,b-a,d-c))

for ne_val in ne_vals:
    dif_all_30_stor = np.empty((0,b-a,d-c))
    dif_all_10_stor = np.empty((0,b-a,d-c))
    for ne in ne_val:

        print 'Endmember #', ne
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
        nyr_h = 2000-his_str_yr+1
        Ih = np.arange(12*(nyr_h-nclim_yr), 12*nyr_h, nmons) 

        # Historical climatology
        sst = fid_his.variables['SST'][Ih,:,a:b,c:d].squeeze()
        lats = fid_his.variables['TLAT'][a:b,c:d].squeeze()
        #dates = fid_his.variables['time'][Ih].squeeze()
        #print dates[0]/365,dates[0]%365,dates[-1]/365,dates[-1]%365
        SST_H_30 = np.mean(sst, axis=0) 
        SST_H_10 = np.mean(sst[10*12:-10*12,:], axis=0)

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
            If = np.arange(int(12*(nyr_f-tot_yr_f[nf])), 12*nyr_f, nmons)
            # Monthly climatology/Annual Mean
            SST_F = np.ma.append(SST_F, fid_fut.variables['SST'][If,:,a:b,c:d].squeeze(), axis=0)
            #dates = fid_fut.variables['time'][If].squeeze() 
            #print dates[0]/365,dates[0]%365,dates[-1]/365,dates[-1]%365
            # Bering Sea Masking #
            for j in range(1,9):
                for i in range(j):
                    SST_F[:,j,i-1] = np.ma.masked
        # Climatology
        SST_F_30 = np.mean(SST_F,axis=0)
        SST_F_10 = np.mean(SST_F[10*12:-10*12,:],axis=0)

        # STORE VALUES FOR ANALYSIS 
        dif_all_30_stor = np.ma.append(dif_all_30_stor, np.ma.expand_dims(SST_F_30-SST_H_30,axis=0), axis=0)
        dif_all_10_stor = np.ma.append(dif_all_10_stor, np.ma.expand_dims(SST_F_10-SST_H_10,axis=0), axis=0)

    # SAVE 2D FIELD FOR FIGURE GENERATION
    if len(ne_val)==1:
       enm_str = str(ne_val[0]).zfill(3)
    else:
       enm_str = 'ALL'

    filstr_30 = enm_str + '_average_' + str(30) + 'yr_dif'
    filstr_10 = enm_str + '_average_' + str(10) + 'yr_dif'

    np.mean(dif_all_30_stor,axis=0).dump(filstr_30)
    np.mean(dif_all_10_stor,axis=0).dump(filstr_10)
