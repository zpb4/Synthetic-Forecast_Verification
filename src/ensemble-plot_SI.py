# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 11:22:54 2024

@author: zpb4
"""
import sys
import os
sys.path.insert(0, os.path.abspath('./src'))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import syn_util
from time import localtime, strftime
from datetime import datetime
import matplotlib.dates as mdates
import seaborn as sns

col_cb = sns.color_palette('colorblind')
#sns.palplot(col_cb)  #plot coloblind palette for comparison
colv1 = sns.color_palette('PuRd',10)
colv2 = sns.color_palette('YlOrBr',10)

cv2 = colv1[4]  #version 1 is pink-purplish color
cv1 = colv2[6]  #version 2 is orangey-brown color
chefs = col_cb[0]  #hefs is colorblind blue

# to simulate a policy, replace 'firo_pool' and 'risk_thresholds' with those found in the train.py file
kcfs_to_tafd = 2.29568411*10**-5 * 86400
K = 317 # TAF
Rmax = 12.5 * kcfs_to_tafd # estimate - from MBK
sd_hefs = '1990-10-01' 
sd_syn = '1985-10-15' 
ed = '2019-08-15'

sd_86 = '1986-01-15' 
ed_86 = '1986-03-10'

loc = 'YRS'
site = 'ORDC1'
vers_out1 = 'test'
vers_out2 = 'test'
has_86 = False

syn_vers1 = 'v1'
syn_vers1_param = 'a'
syn_vers1_setup = '5fold'
syn_path1 = '../Synthetic-Forecast-%s-FIRO-DISES' %(syn_vers1) # path to R synthetic forecast repo for 'r-gen' setting below

syn_vers2 = 'v2'
syn_vers2_pct = 0.99
syn_vers2_setup = '5fold-test'
syn_path2 = '../Synthetic-Forecast-%s-FIRO-DISES' %(syn_vers2) # path to R synthetic forecast repo for 'r-gen' setting below

nsamps = 100

disp_samps = 6  #use 8 for 3 x 3 plot with 1 HEFS and 2:9 syn-hefs

Q_hefs,Qf_hefs,dowy_hefs,tocs,df_idx_hefs = syn_util.extract(sd_hefs,ed,forecast_type='hefs',syn_sample='',Rsyn_path='../Synthetic-Forecast-v1-FIRO-DISES',syn_vers='',forecast_param='',loc=loc,site=site,opt_pcnt=0.99,gen_setup='')
if has_86 == True:
    Q_hefs_86,Qf_hefs_86,dowy_hefs_86,tocs_86,df_idx_hefs_86 = syn_util.extract86(sd_86,ed_86,Rsyn_path='../Synthetic-Forecast-v1-FIRO-DISES',loc=loc,site=site)
Q,Qf,dowy,tocs,df_idx = syn_util.extract(sd_syn,ed,forecast_type='syn',syn_sample='Synth1',Rsyn_path='../Synthetic-Forecast-v1-FIRO-DISES',syn_vers='v1',forecast_param='a',loc=loc,site=site,opt_pcnt=0.99,gen_setup='cal')

Qf_syn_plot_arr1 = np.empty((disp_samps,np.shape(Qf)[0],np.shape(Qf)[1],np.shape(Qf)[2]))

for i in range(disp_samps):
    syn_samp = 'Synth'+str(i+1)
    Q,Qf,dowy,tocs,df_idx = syn_util.extract(sd_syn,ed,forecast_type='syn',syn_sample=syn_samp,Rsyn_path=syn_path1,syn_vers=syn_vers1,forecast_param=syn_vers1_param,loc=loc,site=site,opt_pcnt='',gen_setup=syn_vers1_setup)
    Qf_syn_plot_arr1[i,:,:,:] = Qf
    
Qf_syn_plot_arr2 = np.empty(np.shape(Qf_syn_plot_arr1))

for i in range(disp_samps):
    syn_samp = 'Synth'+str(i+1)
    Q,Qf,dowy,tocs,df_idx = syn_util.extract(sd_syn,ed,forecast_type='syn',syn_sample=syn_samp,Rsyn_path=syn_path2,syn_vers=syn_vers2,forecast_param='',loc=loc,site=site,opt_pcnt=syn_vers2_pct,gen_setup=syn_vers2_setup)
    Qf_syn_plot_arr2[i,:,:,:] = Qf


def fig_title(fig: matplotlib.figure.Figure, txt: str, loc=(0.5,0.98), fontdict=None, **kwargs):
    """Alternative to fig.suptitle that behaves like ax.set_title.
    DO NOT use with suptitle.

    See also:
    https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_title.html
    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.suptitle.html
    https://stackoverflow.com/a/77063164/8954109
    """
    if fontdict is not None:
        kwargs = {**fontdict, **kwargs}
    if "fontsize" not in kwargs and "size" not in kwargs:
        kwargs["fontsize"] = plt.rcParams["axes.titlesize"]

    if "fontweight" not in kwargs and "weight" not in kwargs:
        kwargs["fontweight"] = plt.rcParams["figure.titleweight"]

    if "verticalalignment" not in kwargs:
        kwargs["verticalalignment"] = "top"
    if "horizontalalignment" not in kwargs:
        kwargs["horizontalalignment"] = "center"

    # Tell the layout engine that our text is using space at the top of the figure
    # so that tight_layout does not break.
    # Is there a more direct way to do this?
    fig.suptitle(" ")
    text = fig.text(loc[0], loc[1], txt, transform=fig.transFigure, in_layout=True, **kwargs)

    return text


#---------------------------------------------------
#4x4 forecast plot 
forecast_date = '1996-12-24'
event_date = '1997-01-02'
forecast_idx=df_idx.get_loc(forecast_date)
event_idx=df_idx.get_loc(event_date)
ne,leads = np.shape(Qf)[1:]
ld=event_idx - forecast_idx

st_idx=forecast_idx
end_idx=forecast_idx+leads
dt_idx=df_idx.slice_indexer(df_idx[st_idx],df_idx[end_idx])

forecast_idx_hefs=df_idx_hefs.get_loc(forecast_date)
event_idx_hefs=df_idx_hefs.get_loc(event_date)

st_idx_hefs=forecast_idx_hefs
end_idx_hefs=forecast_idx_hefs+leads
dt_idx_hefs=df_idx_hefs.slice_indexer(df_idx_hefs[st_idx_hefs],df_idx_hefs[end_idx_hefs])

#------------------------------
#1. Ensemble plots
#plotting params
sns.set_theme()
sns.set_style('ticks')
sns.set_context('paper')

if has_86==True:
    fig = plt.figure(layout='constrained',figsize=(10,8))
    gs0 = fig.add_gridspec(2,1)
    gs1 = gs0[0].subgridspec(2,4)
    gs2 = gs0[1].subgridspec(2,4)
    
if has_86==False:
    fig = plt.figure(layout='constrained',figsize=(10,4))
    gs0 = fig.add_gridspec(1,1)
    gs1 = gs0[0].subgridspec(2,4)
    
dt_format=mdates.DateFormatter('%m-%d')
f_arr=np.empty((leads+1))
f_arr[0]=Q[st_idx]


ax1 = fig.add_subplot(gs1[4])
#plot HEFS ensemble
ax1.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
for i in range(ne):
    f_arr[1:]=Qf_hefs[st_idx_hefs,i,:]
    ax1.plot(df_idx_hefs[dt_idx_hefs], f_arr, c=chefs,alpha=0.1)
ax1.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
ax1.xaxis.set_major_formatter(dt_format)
#plt.gca().xaxis.set_ticklabels([])
ax1.tick_params(axis='x',labelrotation=45, labelsize='small')
ax1.set_ylabel('Inflow (TAF/d)')
ax1.set_ylim([0, np.max(Q[dt_idx])*1.5])
ax1.legend(['Obs','HEFS'])
ax1.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
ax1.text(df_idx[event_idx-ld+3],np.max(Q[dt_idx])*1.5*0.9,event_date,fontsize='medium')
ax1.text(df_idx[event_idx-ld],np.max(Q[dt_idx])*1.3,'a)',fontsize='xx-large',fontweight='bold')
#plt.show()

labs = ('b)','c)','d)')
#plot syn-HEFS ensemble
for k in range(3):
    ax2 = fig.add_subplot(gs1[k+1])
    ax2.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr2[k,st_idx,i,:]
        ax2.plot(df_idx[dt_idx], f_arr, c=cv2,alpha=0.2)
    ax2.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    ax2.xaxis.set_ticklabels([])
    if np.isin(k,np.array([1,2])) == True:
        ax2.yaxis.set_ticklabels([])
    if np.isin(k,np.array([0])) == True:
        ax2.legend(['Obs','sHEFS-V1'])
        ax2.set_ylabel('Inflow (TAF/d)')
    ax2.set_ylim([0, np.max(Q[dt_idx])*1.5])
    ax2.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
    ax2.text(df_idx[event_idx-ld],np.max(Q[dt_idx])*1.3,labs[k],fontsize='xx-large',fontweight='bold')

labs2 = ('e)','f)','g)')
#plot syn-HEFS ensemble
for k in range(3):
    ax3 = fig.add_subplot(gs1[k+5])
    ax3.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr2[k+3,st_idx,i,:]
        ax3.plot(df_idx[dt_idx], f_arr, c=cv2,alpha=0.2)
    ax3.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    ax3.yaxis.set_ticklabels([])
    ax3.xaxis.set_major_formatter(dt_format)
    ax3.tick_params(axis='x',labelrotation=45, labelsize='small')
    if np.isin(k,np.array([0])) == True:
        ax3.legend(['Obs','sHEFS-V2'])
    if np.isin(k,np.array([1,2])) == True:
        ax3.yaxis.set_ticklabels([])
    ax3.set_ylim([0, np.max(Q[dt_idx])*1.5])
    ax3.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
    ax3.text(df_idx[event_idx-ld],np.max(Q[dt_idx])*1.3,labs2[k],fontsize='xx-large',fontweight='bold')

if has_86==False:
    fig_title(fig,'%s' %(site),loc=(0.5,1),fontsize='xx-large',fontweight='bold',rotation=0,ha='center')
    fig_title(fig,'1997 Event (hindcast)',loc=(-0.025,0.5),fontsize='xx-large',fontweight='bold',rotation=90,va='center')

    plt.savefig('./figs/comparison/%s/%s/hefs-syn_2x4-ens-plot_SI_synv1%s-synv2%s_forc1=%s_evt1=%s_no86.png' %(loc,site,vers_out1,vers_out2,forecast_date,event_date),dpi=300,bbox_inches='tight')
#-----------------------


sys.modules[__name__].__dict__.clear()
#-----------------------------------------------------end-----------------------------------------------------------------------------