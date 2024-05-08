# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:12:25 2024

@author: zpb4
"""
import sys
import os
sys.path.insert(0, os.path.abspath('./src'))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import syn_util
import ensemble_verification_functions as verify
from time import localtime, strftime
from datetime import datetime

now=datetime.now()
print('ecrps-c start',now.strftime("%H:%M:%S"))

# to simulate a policy, replace 'firo_pool' and 'risk_thresholds' with those found in the train.py file
kcfs_to_tafd = 2.29568411*10**-5 * 86400
K = 317 # TAF
Rmax = 12.5 * kcfs_to_tafd # estimate - from MBK

loc='ADO'
site='ADOC1'
val_samps=6
upr_pcnt = (0.9,1)
lwr_pcnt = (0,0.9)

sd_syn = '1985-10-15' 
ed_syn = '2019-08-15'

idx_syn = pd.date_range(start = sd_syn, end = ed_syn )

sd = '1990-10-01' 
ed = '2019-08-15'
sl_idx = idx_syn.slice_indexer(sd,ed)
save_figs = True  # T to show plots, F to save .png files to ./plot repo
syn_vers1 = 'v1'    # synthetic forecast version; 'v1' or 'v2'
syn_vers1_param = 'a' 
#syn_path1 = 'z:/Synthetic-Forecast-%s-FIRO-DISES/' %(syn_vers1) # path to R synthetic forecast repo for 'r-gen' setting below
syn_path1 = '../Synthetic-Forecast-%s-FIRO-DISES/' %(syn_vers1) # path to R synthetic forecast repo for 'r-gen' setting below
syn_vers2 = 'v2'    # synthetic forecast version; 'v1' or 'v2'
syn_vers2_param = 'l'
#syn_path2 = 'z:/Synthetic-Forecast-%s-FIRO-DISES/' %(syn_vers2) # path to R synthetic forecast repo for 'r-gen' setting below
syn_path2 = '../Synthetic-Forecast-%s-FIRO-DISES/' %(syn_vers2) # path to R synthetic forecast repo for 'r-gen' setting below
gen_path = 'r-gen'
nsamps = 10

vals = np.loadtxt("%s/data/%s/opt_val_years_samp=%s.csv" %(syn_path1,loc,val_samps),skiprows=1) 
val_yrs = np.array(vals)
#val_yrs = (1991,1994,1996,1997,2018,2020)

Q_hefs_inp,Qf_hefs_inp,dowy_hefs_trn,tocs,df_idx_hefs = syn_util.extract(sd,ed,forecast_type='hefs',syn_sample='',Rsyn_path=syn_path1,forecast_param=syn_vers1_param,loc=loc,site=site)
Qf_hefs_trn = verify.onesamp_forecast_rearrange_cumul(Qf_hefs_inp)
Q_hefs_trn = verify.tgt_cumul(Q_hefs_inp, nl=np.shape(Qf_hefs_inp)[2])

wy_vec = df_idx_hefs.year.values
wy_vec[np.isin(df_idx_hefs.month,[10,11,12])] = wy_vec[np.isin(df_idx_hefs.month,[10,11,12])]+1

val_idx = np.arange(len(df_idx_hefs))[np.isin(wy_vec,val_yrs)]
df_idx_val = df_idx_hefs[val_idx]

Qf_v1_inp = np.load('data/%s-%s_Qf_syn-forecast%s_nsamp=%s.npz' %(loc,site,syn_vers1+syn_vers1_param,nsamps))['arr']
Qf_v1_trn = verify.multisamp_forecast_rearrange_cumul(Qf_v1_inp[:,sl_idx,:,:])

Qf_v2_inp = np.load('data/%s-%s_Qf_syn-forecast%s_nsamp=%s.npz' %(loc,site,syn_vers2+syn_vers2_param,nsamps))['arr']
Qf_v2_trn = verify.multisamp_forecast_rearrange_cumul(Qf_v2_inp[:,sl_idx,:,:])

#4. eCRPS diagrams
Qf_hefs = Qf_hefs_trn[val_idx,:,:]
Q_hefs = Q_hefs_trn[val_idx]
Qf_v1 = Qf_v1_trn[:,val_idx,:,:]
Qf_v2 = Qf_v2_trn[:,val_idx,:,:]
dowy_hefs = dowy_hefs_trn[val_idx]

ref_st = '1985-10-15'
ref_end = '2019-08-15'
lds = np.arange(np.shape(Qf_hefs)[2])
sset_idx = np.where((dowy_hefs>60) & (dowy_hefs<170))[0]

Q_ref_inp,dowy_ref = syn_util.extract_obs(ref_st, ref_end, Rsyn_path=syn_path1,loc=loc,site=site)
Qf_ref_ens_inp = verify.create_obs_ref_ens(Q_ref_inp, dowy_ref, sd, ed)

Q_ref = verify.tgt_cumul(tgt = Q_ref_inp, nl = np.shape(Qf_hefs)[2])
Qf_ref_ens_temp = verify.ref_forecast_rearrange_cumul(Qf_ref_ens_inp,nl=np.shape(Qf_hefs)[2])
Qf_ref_ens = Qf_ref_ens_temp[val_idx,:]

ecrps_ss_hefs = np.empty((2,len(lds)))
ecrps_ss_v1 = np.empty((2,len(lds),nsamps))
ecrps_ss_v2 = np.empty((2,len(lds),nsamps))

for i in range(len(lds)):
    Qf_v1_ecrps = Qf_v1[:,:,:,lds[i]]
    Qf_v2_ecrps = Qf_v2[:,:,:,lds[i]]

    ecrps_hefs = verify.onesamp_ecrps(ensemble = Qf_hefs[sset_idx,:,lds[i]], tgt = Q_hefs[sset_idx,lds[i]], pcntile = upr_pcnt)
    ecrps_ref = verify.onesamp_ecrps(ensemble = Qf_ref_ens[sset_idx,:,lds[i]], tgt = Q_hefs[sset_idx,lds[i]], pcntile = upr_pcnt)

    ecrps_v1 = verify.multisamp_ecrps(ensemble = Qf_v1_ecrps[:,sset_idx,:], tgt = Q_hefs[sset_idx,lds[i]], pcntile = upr_pcnt)
    ecrps_v2 = verify.multisamp_ecrps(ensemble = Qf_v2_ecrps[:,sset_idx,:], tgt = Q_hefs[sset_idx,lds[i]], pcntile = upr_pcnt)

    ecrps_ss_hefs[0,i] = 1 - (ecrps_hefs[1] / ecrps_ref[1])
    ecrps_ss_v1[0,i,:] = verify.ecrps_ss(ens_ecrps=ecrps_v1, ref_ecrps=ecrps_ref[1])
    ecrps_ss_v2[0,i,:] = verify.ecrps_ss(ens_ecrps=ecrps_v2, ref_ecrps=ecrps_ref[1])

    ecrps_hefs = verify.onesamp_ecrps(ensemble = Qf_hefs[sset_idx,:,lds[i]], tgt = Q_hefs[sset_idx,lds[i]], pcntile = lwr_pcnt)
    ecrps_ref = verify.onesamp_ecrps(ensemble = Qf_ref_ens[sset_idx,:,lds[i]], tgt = Q_hefs[sset_idx,lds[i]], pcntile = lwr_pcnt)

    ecrps_v1 = verify.multisamp_ecrps(ensemble = Qf_v1_ecrps[:,sset_idx,:], tgt = Q_hefs[sset_idx,lds[i]], pcntile = lwr_pcnt)
    ecrps_v2 = verify.multisamp_ecrps(ensemble = Qf_v2_ecrps[:,sset_idx,:], tgt = Q_hefs[sset_idx,lds[i]], pcntile = lwr_pcnt)

    ecrps_ss_hefs[1,i] = 1 - (ecrps_hefs[1] / ecrps_ref[1])
    ecrps_ss_v1[1,i,:] = verify.ecrps_ss(ens_ecrps=ecrps_v1, ref_ecrps=ecrps_ref[1])
    ecrps_ss_v2[1,i,:] = verify.ecrps_ss(ens_ecrps=ecrps_v2, ref_ecrps=ecrps_ref[1])
    now=datetime.now()
    print(i,now.strftime("%H:%M:%S"))

np.savez_compressed('data/%s-%s_ecrps-ss-cumul_syn-forecast%s_nsamp=%s.npz' %(loc,site,syn_vers1+syn_vers1_param,nsamps), arr=ecrps_ss_v1)
np.savez_compressed('data/%s-%s_ecrps-ss-cumul_syn-forecast%s_nsamp=%s.npz' %(loc,site,syn_vers2+syn_vers2_param,nsamps), arr=ecrps_ss_v2)

np.savez_compressed('data/%s-%s_ecrps-ss-cumul_hefs.npz', arr=ecrps_ss_hefs)

now=datetime.now()
print('ecrps-c end',now.strftime("%H:%M:%S"))