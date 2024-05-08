# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 11:22:54 2024

@author: zpb4
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import model
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
sd_hefs = '1989-10-01' 
sd_syn = '1979-10-02' 
ed = '2019-09-15'

sd_86 = '1986-01-15' 
ed_86 = '1986-03-11'

save_figs = True  # T to show plots, F to save .png files to ./plot repo
syn_vers1 = 'v1'    # synthetic forecast version; 'v1' or 'v2'
syn_vers1_param = 'a'  
syn_path1 = 'z:/Synthetic-Forecast-%s-FIRO-DISES/' %(syn_vers1) # path to R synthetic forecast repo for 'r-gen' setting below
syn_vers2 = 'v2'    # synthetic forecast version; 'v1' or 'v2'
syn_vers2_param = 'a'  
syn_path2 = 'z:/Synthetic-Forecast-%s-FIRO-DISES/' %(syn_vers2) # path to R synthetic forecast repo for 'r-gen' setting below
syn_color = 'pink' # color of synthetic forecast traces to distinguish between the 2 models (magenta for v1, orange for v2)
gen_path = 'r-gen'
nsamps = 100

disp_samps = 3  #use 8 for 3 x 3 plot with 1 HEFS and 2:9 syn-hefs

Q_hefs,Q_MSG_hefs,Qf_hefs,Qf_MSG_hefs,dowy_hefs,tocs,df_idx_hefs = model.extract(sd_hefs,ed,forecast_type='hefs',syn_sample='Synth1',gen_path=gen_path,Rsyn_path=syn_path1)
Q_hefs_86,Q_MSG_hefs_86,Qf_hefs_86,Qf_MSG_hefs_86,dowy_hefs_86,tocs_86,df_idx_hefs_86 = model.extract86(sd_86,ed_86,Rsyn_path=syn_path1)
Q,Q_MSG,Qf,Qf_MSG,dowy,tocs,df_idx = model.extract(sd_syn,ed,forecast_type='syn',syn_sample='Synth1',gen_path=gen_path,Rsyn_path=syn_path1,forecast_param=syn_vers1_param)

Qf_syn_plot_arr1 = np.empty((disp_samps,np.shape(Qf)[0],np.shape(Qf)[1],np.shape(Qf)[2]))
Qf_MSG_syn_plot_arr1 = np.empty(np.shape(Qf_syn_plot_arr1))

for i in range(disp_samps):
    syn_samp = 'Synth'+str(i+1)
    Q,Q_MSG,Qf,Qf_MSG,dowy,tocs,df_idx = model.extract(sd_syn,ed,forecast_type='syn',syn_sample=syn_samp,gen_path=gen_path,Rsyn_path=syn_path1,forecast_param=syn_vers1_param)
    Qf_syn_plot_arr1[i,:,:,:] = Qf
    Qf_MSG_syn_plot_arr1[i,:,:,:] = Qf_MSG
    
Qf_syn_plot_arr2 = np.empty(np.shape(Qf_syn_plot_arr1))
Qf_MSG_syn_plot_arr2 = np.empty(np.shape(Qf_syn_plot_arr1))

for i in range(disp_samps):
    syn_samp = 'Synth'+str(i+1)
    Q,Q_MSG,Qf,Qf_MSG,dowy,tocs,df_idx = model.extract(sd_syn,ed,forecast_type='syn',syn_sample=syn_samp,gen_path=gen_path,Rsyn_path=syn_path2,forecast_param=syn_vers2_param)
    Qf_syn_plot_arr2[i,:,:,:] = Qf
    Qf_MSG_syn_plot_arr2[i,:,:,:] = Qf_MSG


def fig_title(
    fig: matplotlib.figure.Figure, txt: str, loc=(0.5,0.98), fontdict=None, **kwargs
):
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
forecast_date = '1996-12-31'
event_date = '1997-01-03'
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

fig = plt.figure(layout='constrained',figsize=(10,8))
gs0 = fig.add_gridspec(2,1)
gs1 = gs0[0].subgridspec(2,4)
gs2 = gs0[1].subgridspec(2,4)
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
ax1.text(df_idx[event_idx-ld+3],32,event_date,fontsize='medium')
ax1.text(df_idx[event_idx-ld],np.max(Q[dt_idx])*1.3,'a)',fontsize='xx-large',fontweight='bold')
#plt.show()

labs = ('b)','c)','d)')
#plot syn-HEFS ensemble
for k in range(disp_samps):
    ax2 = fig.add_subplot(gs1[k+1])
    ax2.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr1[k,st_idx,i,:]
        ax2.plot(df_idx[dt_idx], f_arr, c=cv1,alpha=0.2)
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
for k in range(disp_samps):
    ax3 = fig.add_subplot(gs1[k+5])
    ax3.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr2[k,st_idx,i,:]
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
    
#-----------------------
#1986 event
forecast_date2 = '1986-02-15'
event_date2 = '1986-02-18'
forecast_idx=df_idx.get_loc(forecast_date2)
event_idx=df_idx.get_loc(event_date2)
ne,leads = np.shape(Qf)[1:]
ld=event_idx - forecast_idx

st_idx=forecast_idx
end_idx=forecast_idx+leads
dt_idx=df_idx.slice_indexer(df_idx[st_idx],df_idx[end_idx])

forecast_idx_hefs=df_idx_hefs_86.get_loc(forecast_date2)
event_idx_hefs=df_idx_hefs_86.get_loc(event_date2)

st_idx_hefs=forecast_idx_hefs
end_idx_hefs=forecast_idx_hefs+leads
dt_idx_hefs=df_idx_hefs.slice_indexer(df_idx_hefs[st_idx_hefs],df_idx_hefs[end_idx_hefs])

#------------------------------
#86 event

ax1 = fig.add_subplot(gs2[4])
#plot HEFS ensemble
ax1.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
for i in range(ne):
    f_arr[1:]=Qf_hefs_86[st_idx_hefs,i,:]
    ax1.plot(df_idx_hefs_86[dt_idx_hefs], f_arr, c=chefs,alpha=0.1)
ax1.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
ax1.xaxis.set_major_formatter(dt_format)
#plt.gca().xaxis.set_ticklabels([])
ax1.tick_params(axis='x',labelrotation=45, labelsize='small')
ax1.set_ylabel('Inflow (TAF/d)')
ax1.set_ylim([0, np.max(Q[dt_idx])*1.5])
ax1.legend(['Obs','HEFS'])
ax1.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
ax1.text(df_idx[event_idx-ld+3],32,event_date2,fontsize='medium')
ax1.text(df_idx[event_idx-ld],np.max(Q[dt_idx])*1.3,'h)',fontsize='xx-large',fontweight='bold')
#plt.show()

labs = ('i)','j)','k)')
#plot syn-HEFS ensemble
for k in range(disp_samps):
    ax2 = fig.add_subplot(gs2[k+1])
    ax2.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr1[k,st_idx,i,:]
        ax2.plot(df_idx[dt_idx], f_arr, c=cv1,alpha=0.2)
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

labs2 = ('l)','m)','n)')
#plot syn-HEFS ensemble
for k in range(disp_samps):
    ax3 = fig.add_subplot(gs2[k+5])
    ax3.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr2[k,st_idx,i,:]
        ax3.plot(df_idx[dt_idx], f_arr, c=cv2,alpha=0.2)
    ax3.plot(df_idx[dt_idx], Q[dt_idx], c = 'black')
    ax3.yaxis.set_ticklabels([])
    ax3.xaxis.set_major_formatter(dt_format)
    ax3.tick_params(axis='x',labelrotation=45, labelsize='small')
    if np.isin(k,np.array([0])) == True:
        ax3.legend(['Obs','sHEFS-V2'])
    plt.ylim([0, np.max(Q[dt_idx])*1.5])
    plt.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
    plt.text(df_idx[event_idx-ld],np.max(Q[dt_idx])*1.3,labs2[k],fontsize='xx-large',fontweight='bold')

fig_title(fig,'1997 Event (In-sample)',loc=(-0.025,0.75),fontsize='xx-large',fontweight='bold',rotation=90,va='center')
fig_title(fig,'1986 Event (Out-of-sample)',loc=(-0.025,0.25),fontsize='xx-large',fontweight='bold',rotation=90,va='center')


plt.savefig('./plot/ensemble/hefs-syn_4x4-ens-plot_forc1=%s_evt1=%s_forc2=%s_evt2=%s.png' %(forecast_date,event_date,forecast_date2,event_date2),dpi=300,bbox_inches='tight')
plt.show()


#---------------------------------------------------
#4x4 cumulativeforecast plot 
forecast_date = '1996-12-31'
event_date = '1997-01-03'
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

def cumul_fun(forecast):
    out = np.zeros_like(forecast)
    for i in range(len(forecast)):
        out[i] = np.sum(forecast[0:i])
    return out

#------------------------------
#1. Ensemble plots
#plotting params
sns.set_theme()
sns.set_style('ticks')
sns.set_context('paper')

fig = plt.figure(layout='constrained',figsize=(10,8))
gs0 = fig.add_gridspec(2,1)
gs1 = gs0[0].subgridspec(2,4)
gs2 = gs0[1].subgridspec(2,4)
dt_format=mdates.DateFormatter('%m-%d')
f_arr=np.empty((leads+1))
f_arr[0]=Q[st_idx]


ax1 = fig.add_subplot(gs1[4])
#plot HEFS ensemble
ax1.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
for i in range(ne):
    f_arr[1:]=Qf_hefs[st_idx_hefs,i,:]
    ax1.plot(df_idx_hefs[dt_idx_hefs], cumul_fun(f_arr), c=chefs,alpha=0.1)
ax1.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
ax1.xaxis.set_major_formatter(dt_format)
#plt.gca().xaxis.set_ticklabels([])
ax1.tick_params(axis='x',labelrotation=45, labelsize='small')
ax1.set_ylabel('Inflow (TAF/d)')
ax1.set_ylim([0, 125])
ax1.legend(['Obs','HEFS'],loc='upper right')
ax1.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
ax1.text(df_idx[event_idx-ld+3],5,event_date,fontsize='medium')
ax1.text(df_idx[event_idx-ld],105,'a)',fontsize='xx-large',fontweight='bold')
#plt.show()

labs = ('b)','c)','d)')
#plot syn-HEFS ensemble
for k in range(disp_samps):
    ax2 = fig.add_subplot(gs1[k+1])
    ax2.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr1[k,st_idx,i,:]
        ax2.plot(df_idx[dt_idx], cumul_fun(f_arr), c=cv1,alpha=0.2)
    ax2.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    ax2.xaxis.set_ticklabels([])
    if np.isin(k,np.array([1,2])) == True:
        ax2.yaxis.set_ticklabels([])
    if np.isin(k,np.array([0])) == True:
        ax2.legend(['Obs','sHEFS-V1'],loc='upper right')
        ax2.set_ylabel('Inflow (TAF/d)')
    ax2.set_ylim([0, 125])
    ax2.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
    ax2.text(df_idx[event_idx-ld],105,labs[k],fontsize='xx-large',fontweight='bold')

labs2 = ('e)','f)','g)')
#plot syn-HEFS ensemble
for k in range(disp_samps):
    ax3 = fig.add_subplot(gs1[k+5])
    ax3.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr2[k,st_idx,i,:]
        ax3.plot(df_idx[dt_idx], cumul_fun(f_arr), c=cv2,alpha=0.2)
    ax3.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    ax3.yaxis.set_ticklabels([])
    ax3.xaxis.set_major_formatter(dt_format)
    ax3.tick_params(axis='x',labelrotation=45, labelsize='small')
    if np.isin(k,np.array([0])) == True:
        ax3.legend(['Obs','sHEFS-V2'])
    if np.isin(k,np.array([1,2])) == True:
        ax3.yaxis.set_ticklabels([])
    ax3.set_ylim([0, 125])
    ax3.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
    ax3.text(df_idx[event_idx-ld],105,labs2[k],fontsize='xx-large',fontweight='bold')
    
#-----------------------
#1986 event
forecast_date2 = '1986-02-15'
event_date2 = '1986-02-18'
forecast_idx=df_idx.get_loc(forecast_date2)
event_idx=df_idx.get_loc(event_date2)
ne,leads = np.shape(Qf)[1:]
ld=event_idx - forecast_idx

st_idx=forecast_idx
end_idx=forecast_idx+leads
dt_idx=df_idx.slice_indexer(df_idx[st_idx],df_idx[end_idx])

forecast_idx_hefs=df_idx_hefs_86.get_loc(forecast_date2)
event_idx_hefs=df_idx_hefs_86.get_loc(event_date2)

st_idx_hefs=forecast_idx_hefs
end_idx_hefs=forecast_idx_hefs+leads
dt_idx_hefs=df_idx_hefs.slice_indexer(df_idx_hefs[st_idx_hefs],df_idx_hefs[end_idx_hefs])

#------------------------------
#86 event

ax1 = fig.add_subplot(gs2[4])
#plot HEFS ensemble
ax1.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
for i in range(ne):
    f_arr[1:]=Qf_hefs_86[st_idx_hefs,i,:]
    ax1.plot(df_idx_hefs_86[dt_idx_hefs], cumul_fun(f_arr), c=chefs,alpha=0.1)
ax1.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
ax1.xaxis.set_major_formatter(dt_format)
#plt.gca().xaxis.set_ticklabels([])
ax1.tick_params(axis='x',labelrotation=45, labelsize='small')
ax1.set_ylabel('Inflow (TAF/d)')
ax1.set_ylim([0, 125])
ax1.legend(['Obs','HEFS'])
ax1.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
ax1.text(df_idx[event_idx-ld+3],5,event_date2,fontsize='medium')
ax1.text(df_idx[event_idx-ld],105,'h)',fontsize='xx-large',fontweight='bold')
#plt.show()

labs = ('i)','j)','k)')
#plot syn-HEFS ensemble
for k in range(disp_samps):
    ax2 = fig.add_subplot(gs2[k+1])
    ax2.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr1[k,st_idx,i,:]
        ax2.plot(df_idx[dt_idx], cumul_fun(f_arr), c=cv1,alpha=0.2)
    ax2.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    ax2.xaxis.set_ticklabels([])
    if np.isin(k,np.array([1,2])) == True:
        ax2.yaxis.set_ticklabels([])
    if np.isin(k,np.array([0])) == True:
        ax2.legend(['Obs','sHEFS-V1'])
        ax2.set_ylabel('Inflow (TAF/d)')
    ax2.set_ylim([0, 125])
    ax2.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
    ax2.text(df_idx[event_idx-ld],105,labs[k],fontsize='xx-large',fontweight='bold')

labs2 = ('l)','m)','n)')
#plot syn-HEFS ensemble
for k in range(disp_samps):
    ax3 = fig.add_subplot(gs2[k+5])
    ax3.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    for i in range(ne):
        f_arr[1:]=Qf_syn_plot_arr2[k,st_idx,i,:]
        ax3.plot(df_idx[dt_idx], cumul_fun(f_arr), c=cv2,alpha=0.2)
    ax3.plot(df_idx[dt_idx], cumul_fun(Q[dt_idx]), c = 'black')
    ax3.yaxis.set_ticklabels([])
    ax3.xaxis.set_major_formatter(dt_format)
    ax3.tick_params(axis='x',labelrotation=45, labelsize='small')
    if np.isin(k,np.array([0])) == True:
        ax3.legend(['Obs','sHEFS-V2'])
    plt.ylim([0, 125])
    plt.axvline(df_idx[event_idx],linewidth=1,linestyle='--')
    plt.text(df_idx[event_idx-ld],105,labs2[k],fontsize='xx-large',fontweight='bold')

fig_title(fig,'1997 Event (In-sample)',loc=(-0.025,0.75),fontsize='xx-large',fontweight='bold',rotation=90,va='center')
fig_title(fig,'1986 Event (Out-of-sample)',loc=(-0.025,0.25),fontsize='xx-large',fontweight='bold',rotation=90,va='center')


plt.savefig('./plot/ensemble/hefs-syn_4x4-ens-plot-cumul_forc1=%s_evt1=%s_forc2=%s_evt2=%s.png' %(forecast_date,event_date,forecast_date2,event_date2),dpi=300,bbox_inches='tight')
plt.show()

#-----------------------------------------------------end-----------------------------------------------------------------------------