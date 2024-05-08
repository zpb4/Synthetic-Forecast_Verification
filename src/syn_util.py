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

    if forecast_type=='hefs':
        path = '%s/out/%s/Qf-%s.nc' % (Rsyn_path,loc,forecast_type)
    else:
        path = '%s/out/%s/Qf-%s.nc' % (Rsyn_path,loc,forecast_type+forecast_param)
        
    da = xr.open_dataset(path)[forecast_type]
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 

    dowy = np.array([water_day(d) for d in df.index.dayofyear])
    tocs = get_tocs(dowy)
    
    ado = {'ADOC1':0}
    nhg = {'MSGC1L':0,'NHGC1':1}
    lam = {'HOPC1L':0,'LAMC1':1,'UKAC1':2}
    yrs = {'MRYC1L':0,'NBBC1':1,'ORDC1':2}
    locs = {'ADO':ado,'NHG':nhg,'LAM':lam,'YRS':yrs}
    
    site_id = locs[loc][site]
    
    # (ensemble: 4, site: 2, date: 15326, trace: 42, lead: 15)
    if forecast_type == 'hefs':
        Qf = da.sel(ensemble=0, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    if forecast_type == 'syn':
        Qf = da.sel(ensemble=int(syn_sample[5:])-1, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    
    #recommend not presorting ensembles because it will mix and match ensemble members
    #Qf.sort(axis = 1)
    #Qf_MSG.sort(axis = 1)
    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx

def extract86(sd,ed,Rsyn_path,loc,site):
    path = '%s/out/%s/Qf-hefs86.nc' % (Rsyn_path,loc)
    da = xr.open_dataset(path)['hefs']
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 
    
    ado = {'ADOC1':0}
    nhg = {'MSGC1L':0,'NHGC1':1}
    lam = {'HOPC1L':0,'LAMC1':1,'UKAC1':2}
    yrs = {'MRYC1L':0,'NBBC1':1,'ORDC1':2}
    locs = {'ADO':ado,'NHG':nhg,'LAM':lam,'YRS':yrs}
    
    site_id = locs[loc][site]

    dowy = np.array([water_day(d) for d in df.index.dayofyear])
    tocs = get_tocs(dowy)
    
    Qf = da.sel(ensemble=0, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)

    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx


def get_tocs(d): 
  tp = [0, 60, 170, 252, 366]
  sp = [K, 0.5*K, 0.5*K, K, K] # Note: actual NHG tocs min is 152 TAF, about 0.48 * K
  return np.interp(d, tp, sp)





