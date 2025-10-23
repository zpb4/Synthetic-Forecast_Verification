

rm(list=ls());gc()
library(fields)
library(scales)
library(zoo)
library(abind)
#----------------------------------------
setwd('z:/Synthetic-Forecast_Verification/')

syn_vers = 2
loc = 'SOD'
opt_site = 'SRWC1'
disp_site = 'SRWC1'
opt_pcnt = 0.99
cal_val_setup = '5fold' # 'cal' '5fold' '5fold-test'
opt_strat = 'ecrps-dts'
obj_pwr = 0
has_86 = F

if (!dir.exists(paste('e:/Projects/FIRO/firo_syn-forecast_production/figs/',disp_site,'/',cal_val_setup,opt_pcnt,opt_strat,obj_pwr,
'/',sep=''))) {
  dir.create(paste('e:/Projects/FIRO/firo_syn-forecast_production/figs/',disp_site,'/',cal_val_setup,opt_pcnt,opt_strat,obj_pwr,
                   '/',sep=''),recursive=T)
}

path_out = paste('e:/Projects/FIRO/firo_syn-forecast_production/figs/',disp_site,'/',cal_val_setup,opt_pcnt,opt_strat,obj_pwr,
      '/',sep='')
#plot setup
disp_pcnt <- 0.99
n_samp <- 3
obs_rank <- 3  #pick which obs event to plot (1 largest, 2 second largest, etc)
#LAMC1: 1 = 2005 evt, 2 = 1995 evt, 3 = 2019 evt
seed<-1
path = paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

if(has_86==T){
  load(paste(path,'out/',loc,'/data_prep_rdata86.RData',sep=''))
  cur_site <- which(site_names==disp_site)
  idx_site <- which(site_names==opt_site)
  ixx_hefs86 <- ixx_hefs
  hefs_fwd_86 <- hefs_forward[cur_site,,,]}

load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))
cur_site <- which(site_names==disp_site)
idx_site <- which(site_names==opt_site)

syn_hefs_forward <- readRDS(paste(path,'out/',loc,'/syn_hefs_forward_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',opt_site,'_',cal_val_setup,'_plot-ens.rds',sep=''))
hefs_scale <- readRDS(file=paste(path,'out/',loc,'/hefs-scale_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',opt_site,'_',cal_val_setup,'.rds',sep=''))
hefs_resamps <- readRDS(file=paste(path,'out/',loc,'/hefs-resamps_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',opt_site,'_',cal_val_setup,'.rds',sep=''))
shefs_fwd <- syn_hefs_forward[,cur_site,,,]
hefs_sc <- hefs_scale[,cur_site,,]
hefs_fwd_sset <- hefs_forward[cur_site,,,]

rm(hefs_scale,syn_hefs_forward,hefs_forward,hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,obs_forward_all_leads_hind)
gc()

hefs_fwd <- shefs_fwd[1,,,]
hefs_idx <- ixx_obs_forward%in%ixx_hefs
hefs_fwd[,hefs_idx,] <- hefs_fwd_sset

if(has_86==T){
  hefs86_idx <- ixx_obs_forward%in%ixx_hefs86
  hefs_fwd[,hefs86_idx,] <- hefs_fwd_86
  ixx_keep <- ixx_obs_forward%in%c(ixx_hefs86,ixx_hefs)}

if(has_86==F){
  ixx_keep <- ixx_obs_forward%in%c(ixx_hefs)}

source('./src/forecast_verification_functions.R')

ixx_obs_forward_wy <- wy_fun(ixx_obs_forward)

##if(loc=='SOD'){
  ##hefs_fwd[,ixx_obs_forward_wy$year%in%97:100,] <- hefs_fwd[,ixx_obs_forward_wy$year%in%93:96,]
##}

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
##############Ensemble Plots###############
n_evts = 10
sep = 15

obs_extract <- obs[,idx_site]
if(has_86==T){
  obs_extract[ixx_obs%in%c(ixx_hefs,ixx_hefs86)] <- 0}
if(has_86==F){
  obs_extract[ixx_obs%in%ixx_hefs] <- 0}
obs_evt_idx <- declust_evts_extract(obs_extract,n_evts,sep,max_lds = 15)
obs_gen <- obs[,cur_site]
obs_fwd_gen <- obs_forward_all_leads[cur_site,,]

obs_dates <- ixx_obs[obs_evt_idx]

#load packages
library(tidyverse)
library(data.table)
library(gridExtra)
library(ks)
library(viridis)
library(RColorBrewer)

#setup plot parameters
clrs<-palette.colors()

##########################sHEFS plot###########################################
lds = c(10,5,3,1)
hefs_plts = vector('list',length(lds))
ylm_scale = 1.5
plt_evts = 10

samps = 3
shefs_plts = vector('list',length(lds)*samps)
#ylm_scale = 1.5
#show_x = T
show_x = c(rep(F,(samps-1)*length(lds)),rep(T,length(lds)))
show_dtg = c(rep(T,length(lds)),rep(F,(samps-1)*length(lds)))

set.seed(seed)
samp_idx <- sample(1:10,samps);samp_vec <- rep(samp_idx,each=length(lds))
ldsv <- rep(lds,samps)

for(k in 1:plt_evts){
for(i in 1:(length(lds)*samps)){
  hefs_plt_dt <- data.table(x=0:leads,hefs=rbind(rep(obs_gen[obs_evt_idx[k]-ldsv[i]],dim(shefs_fwd)[2]),t(shefs_fwd[samp_vec[i],,obs_evt_idx[k]-ldsv[i],])),
                            obs=c(obs_gen[obs_evt_idx[k]-ldsv[i]],obs_fwd_gen[obs_evt_idx[k]-ldsv[i],]))
  
  ylm = ylm_scale * max(hefs_plt_dt$obs)
  #convert to long form
  df <- melt(hefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  df$col <- factor(ifelse(df$variable=="obs",1,2),levels=1:2,labels=c("obs","shefs"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=col,linewidth=col,alpha=col))+
    scale_color_manual(values=c('obs'=clrs[[1]],'shefs'=clrs[[7]]))+
    scale_linewidth_manual(values=c('obs'=1.5,'shefs'=.75))+
    scale_alpha_manual(values=c('obs'=1,'shefs'=.2))+
    scale_x_continuous(breaks=0:leads,labels=0:leads)+
    geom_vline(xintercept = ldsv[i],linetype='dotted',linewidth=0.5)+
    {if(show_dtg[i]==T & lds[i]>5 & i!=1)annotate('text',x=lds[i],y=0.95*ylm,label=obs_dates[k],size=4,hjust=1.15)}+
    {if(show_dtg[i]==T & lds[i]<=5 & i!=1)annotate('text',x=lds[i],y=0.95*ylm,label=obs_dates[k],size=4,hjust=-0.15)}+
    {if(i==1)annotate('text',x=leads,y=0.95*ylm,label=obs_dates[k],size=4,hjust=1.15)}+
    #annotate('text',x=1,y=850,label='b)',size=6)+
    labs(x='forecast lead (hrs)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,leads),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1),
          panel.grid.major.x = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,.9),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=1)theme(legend.position = 'none')}+
    {if(ldsv[i]!=lds[1])
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x[i]==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  shefs_plts[[i]] <- hefs_ens_plt}

layout_mat <- matrix(1:(samps*length(lds)),nrow=samps,ncol=length(lds),byrow=T)
shefs_gplot<-marrangeGrob(shefs_plts,nrow=samps,ncol=length(lds),top='',layout_matrix = layout_mat)
#shefs_gplot

ggsave(paste(path_out,loc,'_',disp_site,'_shefs-ens-plot-xhc_setup=',cal_val_setup,'_pct=',opt_pcnt,'_pwr=',obj_pwr,'_strat=',opt_strat,'_evt=',k,'_lds=',str_flatten(lds,collapse='-'),'_seed=',seed,'.png',sep=''),shefs_gplot,dpi=320,width=3*length(lds),height=2.5*(samps),unit='in')

#////////////////////////////////////////////////////////////////////////////

##########################sHEFS sample plot###########################################
if(loc=='SOD'){
  ixx_obs_wy <- wy_fun(ixx_obs)
  obs_gen[ixx_obs_wy$year%in%97:100] <- obs_gen[ixx_obs_wy$year%in%93:96]
  obs_fwd_gen[ixx_obs_forward_wy$year%in%97:100,] <- obs_fwd_gen[ixx_obs_forward_wy$year%in%93:96,]
}

samps = 3
shefs_plts = vector('list',length(lds)*samps)
#ylm_scale = 1.5
#show_x = T
show_x = c(rep(F,(samps-1)*length(lds)),rep(T,length(lds)))

set.seed(seed)
samp_idx <- sample(1:10,samps);samp_vec <- rep(samp_idx,each=length(lds))
ldsv <- rep(lds,samps)

for(i in 1:(length(lds)*samps)){
  hefs_samp_idx <- which(as.character(ixx_obs_forward)==hefs_resamps[samp_vec[i],obs_evt_idx[k]-ldsv[i]])
  hefs_plt_dt <- data.table(x=0:leads,hefs=rbind(rep(obs_gen[hefs_samp_idx],dim(shefs_fwd)[2]),t(hefs_fwd[,hefs_samp_idx,])),
                                      obs=c(obs_gen[hefs_samp_idx],obs_fwd_gen[hefs_samp_idx,]))
  
  ylm = ylm_scale * max(hefs_plt_dt$obs)
  #convert to long form
  df <- melt(hefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  df$col <- factor(ifelse(df$variable=="obs",1,2),levels=1:2,labels=c("obs","hefs"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=col,linewidth=col,alpha=col))+
    scale_color_manual(values=c('obs'=clrs[[1]],'hefs'=clrs[[9]]))+
    scale_linewidth_manual(values=c('obs'=1.5,'hefs'=.75))+
    scale_alpha_manual(values=c('obs'=1,'hefs'=.35))+
    scale_x_continuous(breaks=0:leads,labels=0:leads)+
    geom_vline(xintercept = ldsv[i],linetype='dotted',linewidth=0.5)+
    {if(ldsv[i]>5)annotate('text',x=ldsv[i],y=0.95*ylm,label=ixx_obs_forward[hefs_samp_idx+ldsv[i]],size=4,hjust=1.15)}+
    {if(ldsv[i]<=5)annotate('text',x=ldsv[i],y=0.95*ylm,label=ixx_obs_forward[hefs_samp_idx+ldsv[i]],size=4,hjust=-0.15)}+
    #annotate('text',x=1,y=850,label='b)',size=6)+
    labs(x='forecast lead (hrs)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,leads),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1),
          panel.grid.major.x = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.1,.9),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=1)theme(legend.position = 'none')}+
    #{if(ldsv[i]!=lds[1])
      #theme(axis.text.y=element_blank(),
            #axis.title.y=element_blank())}+
    {if(show_x[i]==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  shefs_plts[[i]] <- hefs_ens_plt}

layout_mat <- matrix(1:(samps*length(lds)),nrow=samps,ncol=length(lds),byrow=T)
shefs_gplot<-marrangeGrob(shefs_plts,nrow=samps,ncol=length(lds),top='',layout_matrix = layout_mat)
#shefs_gplot

ggsave(paste(path_out,loc,'_',disp_site,'_shefs-hefs-samp-plot-xhc_setup=',cal_val_setup,'_pct=',opt_pcnt,'_pwr=',obj_pwr,'_strat=',opt_strat,'_evt=',k,'_lds=',str_flatten(lds,collapse='-'),'_seed=',seed,'.png',sep=''),shefs_gplot,dpi=320,width=3*length(lds),height=2.5*(samps),unit='in')


##########################sHEFS scale plot###########################################
samps = 3
shefs_plts = vector('list',length(lds)*samps)
#ylm_scale = 1.5
#show_x = T
show_x = c(rep(F,(samps-1)*length(lds)),rep(T,length(lds)))

set.seed(seed)
samp_idx <- sample(1:10,samps);samp_vec <- rep(samp_idx,each=length(lds))
ldsv <- rep(lds,samps)

for(i in 1:(length(lds)*samps)){
  hefs_plt_dt <- data.table(x=1:leads,hefs=hefs_sc[samp_vec[i],obs_evt_idx[k]-ldsv[i],])
  
  ylm = 1.1 * max(hefs_plt_dt$hefs)
  #convert to long form
  df <- melt(hefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value),linewidth=2,color='black')+
    scale_x_continuous(breaks=1:leads,labels=1:leads)+
    geom_vline(xintercept = ldsv[i],linetype='dotted',linewidth=0.5)+
    geom_hline(yintercept = 1,linetype='dotted',linewidth=0.5,color='orange4')+
    labs(x='forecast lead (hrs)',y='scale factor')+
    {if(i%in%c(1:length(lds)))
      annotate('text',x=lds[i],y=0.95*ylm,label=obs_dates[k],size=4,hjust=-0.05)}+
    coord_cartesian(xlim=c(0,leads),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          panel.grid.major.x = element_line(linewidth=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,.9),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=1)theme(legend.position = 'none')}+
    #{if(ldsv[i]!=lds[1])
      #theme(axis.text.y=element_blank(),
            #axis.title.y=element_blank())}+
    {if(show_x[i]==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  shefs_plts[[i]] <- hefs_ens_plt}

layout_mat <- matrix(1:(samps*length(lds)),nrow=samps,ncol=length(lds),byrow=T)
shefs_gplot<-marrangeGrob(shefs_plts,nrow=samps,ncol=length(lds),top='',layout_matrix = layout_mat)
#shefs_gplot

ggsave(paste(path_out,loc,'_',disp_site,'_shefs-scale-plot-xhc_setup=',cal_val_setup,'_pct=',opt_pcnt,'_pwr=',obj_pwr,'_strat=',opt_strat,'_evt=',k,'_lds=',str_flatten(lds,collapse='-'),'_seed=',seed,'.png',sep=''),shefs_gplot,dpi=320,width=3*length(lds),height=2.5*(samps),unit='in')

}


##############################END###############################################################

