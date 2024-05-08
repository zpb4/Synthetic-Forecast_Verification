import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from util import water_day
import cvxpy as cvx
from numba import njit
#from scipy.optimize import differential_evolution as DE, minimize

# using cvxpy for the inner loop prevents numba compiling

kcfs_to_tafd = 2.29568411*10**-5 * 86400
K = 317 # TAF
#Rmax = 12.5 * kcfs_to_tafd # estimate - from MBK, this is correct
Rmax = 4 * kcfs_to_tafd

#ramp rate from WCM is <= 2000 cfs per 2 hr period
#ramping_rate = 12 * kcfs_to_tafd
ramping_rate = 2 * kcfs_to_tafd

def extract_obs(sd,ed,Rsyn_path,loc,site):
    
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 

    dowy = np.array([water_day(d) for d in df.index.dayofyear])
    
    return Q,dowy
    
def extract(sd,ed,forecast_type,syn_sample,Rsyn_path,forecast_param,loc,site):

    
 
    if forecast_type=='syn':
        path = '%s/out/Qf-%s.nc' % (Rsyn_path,forecast_type)
    else:
        path = '%s/out/Qf-%s.nc' % (Rsyn_path,forecast_type)
        
    da = xr.open_dataset(path)[forecast_type]
    df = pd.read_csv('%s/data/observed_flows.csv' %(Rsyn_path), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df['ADOC1'].values 

    dowy = np.array([water_day(d) for d in df.index.dayofyear])
    tocs = get_tocs(dowy)
    #dowy -=1
    # if the start date is in a leap year adjust by 1
    
    # (ensemble: 4, site: 2, date: 15326, trace: 42, lead: 15)
    if forecast_type == 'hefs':
        if gen_path == 'py-gen':
            Qf = da.sel(ensemble='HEFS', site='ADOC1', date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
        if gen_path == 'r-gen':
            Qf = da.sel(ensemble=0, site=0, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    if forecast_type == 'syn':
        if gen_path == 'py-gen':
            Qf = da.sel(ensemble=syn_sample, site='ADOC1', date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
        if gen_path == 'r-gen':
            Qf = da.sel(ensemble=int(syn_sample[5:])-1, site=0, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    
    #recommend not presorting ensembles because it will mix and match ensemble members
    #Qf.sort(axis = 1)
    #Qf_MSG.sort(axis = 1)
    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx

def extract86(sd,ed,Rsyn_path):
    path = '%s/out/YRS/Qf-hefs86.nc' % (Rsyn_path)
    da = xr.open_dataset(path)['hefs']
    df = pd.read_csv('%s/data/YRS/observed_flows.csv' %(Rsyn_path), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q_MRY = df['MRYC1L'].values 
    Q_NBB = df['NBBC1'].values
    Q_ORD = df['ORDC1'].values

    dowy = np.array([water_day(d) for d in df.index.dayofyear])
    tocs = get_tocs(dowy)
    
    Qf_MRY = da.sel(ensemble=0, site=0, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    Qf_NBB = da.sel(ensemble=0, site=1, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    Qf_ORD = da.sel(ensemble=0, site=2, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)

    df_idx = df.index
    
    return Q_MRY,Q_NBB,Q_ORD,Qf_MRY,Qf_NBB,Qf_ORD,dowy,tocs,df_idx


def get_tocs(d): 
  tp = [0, 60, 170, 252, 366]
  sp = [K, 0.5*K, 0.5*K, K, K] # Note: actual NHG tocs min is 152 TAF, about 0.48 * K
  return np.interp(d, tp, sp)





