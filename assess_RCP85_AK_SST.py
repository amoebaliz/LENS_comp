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
a = 313
b = 330
c = 206
d = 239
#nlat,313,329 -d nlon,206,238

nmons=1

# Number of Years in climatology
nclim_yr = 10
nclim_yr = 30

# Number of NCAR LE endmembers
end_mem = 35

# EACH MONTH (if climatology)
for nm in range(nmons):

    # STORAGE ARRAYS FOR MONTHLY CLIM MEAN/STD
    dif_30_all_stor = np.empty((0,b-a,d-c))
    dif_10_all_stor = np.empty((0,b-a,d-c))

    for ne in range(1,end_mem+1):
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
        Ih = np.arange(nm + 12*(nyr_h-nclim_yr), 12*nyr_h, nmons) 

        # Historical climatology
        sst = fid_his.variables['SST'][Ih,:,a:b,c:d].squeeze()
        print sst.shape
        dates = fid_his.variables['time'][Ih].squeeze()
        print dates[0]/365,dates[0]%365,dates[-1]/365,dates[-1]%365
        for j in range(1,9):
            for i in range(j):
                sst[:,j,i-1] = np.ma.masked

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
            If = np.arange(nm + int(12*(nyr_f-tot_yr_f[nf])), 12*nyr_f, nmons)
            # Monthly climatology/Annual Mean
            SST_F = np.ma.append(SST_F, fid_fut.variables['SST'][If,:,a:b,c:d].squeeze(), axis=0)
            dates = fid_fut.variables['time'][If].squeeze() 
            print dates[0]/365,dates[0]%365,dates[-1]/365,dates[-1]%365
            # Bering Sea Masking #
            for j in range(2,11):
                for i in range(j-1):
                    SST_F[:,j,i] = np.ma.masked
        # Climatology
        SST_F_30 = np.mean(SST_F,axis=0)
        SST_F_10 = np.mean(SST_F[10*12:-10*12,:],axis=0)
        # STORE VALUES FOR ANALYSIS 
        dif_30_all_stor = np.ma.append(dif_30_all_stor, np.ma.expand_dims(SST_F_30-SST_H_30,axis=0), axis=0)
        dif_10_all_stor = np.ma.append(dif_10_all_stor, np.ma.expand_dims(SST_F_10-SST_H_10,axis=0), axis=0)

    avg_dif_30_sst = np.ma.mean(np.ma.mean(dif_30_all_stor,axis=2),axis=1).squeeze()
    avg_dif_10_sst = np.ma.mean(np.ma.mean(dif_10_all_stor,axis=2),axis=1).squeeze()
    #print avg_dif_sst
    ##########
    # OUTPUT #
    ##########
    print "COLDEST 30yr"
    #print 'ABS', np.argsort(avg_fut_sst)[:5]+1
    print 'DIF', np.argsort(avg_dif_30_sst.squeeze())[:5]+1

    print "COLDEST 10yr"
    #print 'ABS', np.argsort(avg_fut_sst)[:5]+1
    print 'DIF', np.argsort(avg_dif_10_sst.squeeze())[:5]+1

    print "WARMEST 30yr"
    #print 'ABS', np.argsort(avg_fut_sst)[-5:]+1
    print 'DIF', np.argsort(avg_dif_30_sst.squeeze())[-10:]+1

    print "WARMEST 10yr"
    #print 'ABS', np.argsort(avg_fut_sst)[-5:]+1
    print 'DIF', np.argsort(avg_dif_10_sst.squeeze())[-10:]+1


