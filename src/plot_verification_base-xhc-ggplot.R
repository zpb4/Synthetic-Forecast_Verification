#Script to print out plots of verification metrics calculated in 'calc_ecrps-rankhist_base_xhindcast.R' script
#Specifications here much match those of the calculation script, or it will need to be rerun
#The x-hindcast verification statistics are compared to the hindcast verification statistics for reference

#Note: The current 'ylm' specifications for each plot are currently hard set and may need to be modified for 
#different sites; will update this to adjust to the data

#Load packages
rm(list=ls());gc()
library(fields)
library(scales)
library(zoo)
library(scales)
library(tidyverse)
library(data.table)
library(gridExtra)
library(ks)
library(viridis)
library(RColorBrewer)
library(latex2exp)
clrs<-palette.colors()

#root directory
setwd('z:/Synthetic-Forecast_Verification/')

#Primary modifiable input parameters
#////////////////////////////////////////////////////////////////////////////////
#location and site info
loc = 'YRS'             #overall location
opt_site = 'ORDC1'      #keysite for the synthetic forecasting run
disp_site = 'ORDC1'     #what site you want to display

#synthetic forecast setup specifics; should match generation setup you want to look at
syn_vers = 2            #which synthetic version to use (probably 2)
opt_pcnt = 0.99         #what percentile of the data was the synthetic forecast optimized to
cal_val_setup = 'cal'   #what was the optimization setup? options: 'cal' '5fold' '5fold-test'
opt_strat = 'ecrps-dts' #what was the loss function strateg? default: 'ecrps-dts'
obj_pwr = 0             #what was the objection function weighting across leads? default: 0 
has_86 = T              #does the HEFS training data include 1986 special run?

#plot setup (must match what was calculated in the 'calc_ecrps-rankhist_base.R' script)
disp_pcnt <- 0.99

#path to synthetic forecast repo
path = paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

#output path
path_out = paste('e:/Projects/FIRO/firo_syn-forecast_production/figs/',disp_site,'/',cal_val_setup,opt_pcnt,opt_strat,obj_pwr,
                 '/',sep='')

#//////////////////////////////////////////////////////////////////////////////////////////////////////

if (!dir.exists(path_out)) {
  dir.create(path_out,recursive=T)
}

if(has_86==T){
  load(paste(path,'out/',loc,'/data_prep_rdata86.RData',sep=''))
  cur_site <- which(site_names==disp_site)
  idx_site <- which(site_names==opt_site)
  ixx_hefs86 <- ixx_hefs}

load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward,hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,obs_forward_all_leads_hind)
gc()


hefs_ecrps_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-ecrps-vec.rds',sep=''))
shefs_ecrps_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-vec_xhc.rds',sep=''))
hefs_ecrps_ss<-readRDS(paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-ecrps-ss.rds',sep=''))
shefs_ecrps_ss<-readRDS(paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-ss_xhc.rds',sep=''))
hefs_rank_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-rank-vec.rds',sep=''))
shefs_rank_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-rank-vec_xhc.rds',sep=''))
hefs_mse_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-mse-vec.rds',sep=''))
shefs_mse_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-mse-vec_xhc.rds',sep=''))
hefs_sdev_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-sdev-vec.rds',sep=''))
shefs_sdev_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-sdev-vec_xhc.rds',sep=''))
hefs_pbias_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-pbias-vec.rds',sep=''))
shefs_pbias_vec<-readRDS(paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-pbias-vec_xhc.rds',sep=''))
hefs_ecrps_pk<-readRDS(paste('./data/',loc,'-',disp_site,'_hefs-ecrps-peak10.rds',sep=''))
shefs_ecrps_pk<-readRDS(paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-peak10_xhc.rds',sep=''))
hefs_ecrps_ss_pk<-readRDS(paste('./data/',loc,'-',disp_site,'_hefs-ecrps-ss-peak10.rds',sep=''))
shefs_ecrps_ss_pk<-readRDS(paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-ss-peak10_xhc.rds',sep=''))


source('./src/forecast_verification_functions.R')

obs_eval <- readRDS(paste('./data/',loc,'-',disp_site,'_obs-eval.rds',sep=''))
obs_key <- readRDS(paste('./data/',loc,'-',opt_site,'_obs-key.rds',sep=''))
ixx_eval <- readRDS(paste('./data/',loc,'-',disp_site,'_ixx-eval.rds',sep=''))

obs_eval_xhc <- readRDS(paste('./data/',loc,'-',disp_site,'_obs-eval_xhc.rds',sep=''))
obs_key_xhc <- readRDS(paste('./data/',loc,'-',opt_site,'_obs-key_xhc.rds',sep=''))
ixx_eval_xhc <- readRDS(paste('./data/',loc,'-',disp_site,'_ixx-eval_xhc.rds',sep=''))

#///////////////////////////////////////////////////////////////////////////////////////////////////////////

#########################################eCRPS boxplot##############################################
lds = rbind(1:5,6:10,c(11:14,14))
ecrps_plts = vector('list',dim(lds)[1])
show_x = T
ylm = c(35,65,75)
#ylm = c(3,3,3)

for(l in 1:dim(lds)[1]){
  shefs_ecrps_in = apply(shefs_ecrps_vec[,,lds[l,]],3,as.vector)
  hefs_ecrps_in = hefs_ecrps_vec[,lds[l,]]
  
  if(dim(shefs_ecrps_vec)[2]>dim(hefs_ecrps_in)[1]){
    hefs_ecrps_in <- rbind(hefs_ecrps_in,array(NA,c(dim(shefs_ecrps_vec)[2]-dim(hefs_ecrps_in)[1],dim(hefs_ecrps_in)[2])))
  }
  
  if(dim(shefs_ecrps_vec)[2]<dim(hefs_ecrps_in)[1]){
    shefs_sset = shefs_ecrps_vec[,,lds[l,]]
    shefs_in = abind(shefs_sset,array(NA,c(dim(shefs_sset)[1],dim(hefs_ecrps_in)[1] - dim(shefs_ecrps_vec)[2],dim(shefs_sset)[3])),along=2)
    shefs_ecrps_in = apply(shefs_in,3,as.vector)
  }
  
  ecrps_dt <- data.table(x=1,hefs=hefs_ecrps_in,syn=shefs_ecrps_in)

  df <- melt(ecrps_dt,id='x')

  df$x <- rep(rep(lds[l,],each=dim(shefs_ecrps_in)[1]),2)
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))
  
  ylm[l] = sort(as.vector(hefs_ecrps_vec[,lds[l,]]))[round(0.9*length(as.vector(hefs_ecrps_vec[,lds[l,]])))]

  xlbs = paste('ld',lds)

  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=clrs[[9]],'shefs'=clrs[[3]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    scale_x_discrete(breaks=lds,labels=xlbs,name='')+
    scale_y_continuous(name = 'CRPS')+
    coord_cartesian(ylim = c(0,ylm[l]))+
    #labs(title=paste('ld',lds[l]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,1.1),
            legend.direction = 'horizontal',
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none')}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    guides(fill = guide_legend(byrow = F))
      

  ecrps_plts[[l]] <- ecrps_bplt}
  
#mat<-matrix(1:length(lds),ncol=length(lds))
ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=dim(lds)[1],ncol=1,top='')
#ecrps_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_ecrps-bplots_pcntile=',disp_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),ecrps_gplot,dpi=320,width=4.5,height=6,unit='in')

#########################################eCRPS-SS boxplot##############################################
lds = rbind(1:5,6:10,c(11:14,14))
ecrps_plts = vector('list',dim(lds)[1])
show_x = T

for(l in 1:dim(lds)[1]){
  shefs_ecrps_in = apply(shefs_ecrps_ss[,,lds[l,]],3,as.vector)
  hefs_ecrps_in = hefs_ecrps_ss[,lds[l,]]
  
  if(dim(shefs_ecrps_ss)[2]>dim(hefs_ecrps_in)[1]){
    hefs_ecrps_in <- rbind(hefs_ecrps_in,array(NA,c(dim(shefs_ecrps_ss)[2]-dim(hefs_ecrps_in)[1],dim(hefs_ecrps_in)[2])))
  }
  
  if(dim(shefs_ecrps_ss)[2]<dim(hefs_ecrps_in)[1]){
    shefs_sset = shefs_ecrps_ss[,,lds[l,]]
    shefs_in = abind(shefs_sset,array(NA,c(dim(shefs_sset)[1],dim(hefs_ecrps_in)[1] - dim(shefs_ecrps_ss)[2],dim(shefs_sset)[3])),along=2)
    shefs_ecrps_in = apply(shefs_in,3,as.vector)
  }
  
  ecrps_dt <- data.table(x=1,hefs=hefs_ecrps_in,syn=shefs_ecrps_in)
  
  df <- melt(ecrps_dt,id='x')
  
  df$x <- rep(rep(lds[l,],each=dim(shefs_ecrps_in)[1]),2)
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))

  xlbs = paste('ld',lds)
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=clrs[[9]],'shefs'=clrs[[3]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    scale_x_discrete(breaks=lds,labels=xlbs,name='')+
    scale_y_continuous(name = 'CRPS-SS')+
    coord_cartesian(ylim = c(0,1),expand=T)+
    #labs(title=paste('ld',lds[l]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,1.1),
            legend.direction = 'horizontal',
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none')}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    guides(fill = guide_legend(byrow = F))
  
  
  ecrps_plts[[l]] <- ecrps_bplt}

#mat<-matrix(1:length(lds),ncol=length(lds))
ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=dim(lds)[1],ncol=1,top='')
#ecrps_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_ecrps-ss-bplots_pcntile=',disp_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),ecrps_gplot,dpi=320,width=4.5,height=6,unit='in')


#########################################mse boxplot##############################################
lds = rbind(1:5,6:10,c(11:14,14))
ecrps_plts = vector('list',dim(lds)[1])
show_x = T
ylm = c(1000,5000,6000)
#ylm = c(2,2,2)

for(l in 1:dim(lds)[1]){
  shefs_mse_in = apply(shefs_mse_vec[,,lds[l,]],3,as.vector)
  hefs_mse_in = hefs_mse_vec[,lds[l,]]
  
  if(dim(shefs_mse_vec)[2]>dim(hefs_mse_in)[1]){
    hefs_mse_in <- rbind(hefs_mse_in,array(NA,c(dim(shefs_mse_vec)[2]-dim(hefs_mse_in)[1],dim(hefs_mse_in)[2])))
  }
  
  if(dim(shefs_mse_vec)[2]<dim(hefs_mse_in)[1]){
    shefs_sset = shefs_mse_vec[,,lds[l,]]
    shefs_in = abind(shefs_sset,array(NA,c(dim(shefs_sset)[1],dim(hefs_mse_in)[1] - dim(shefs_mse_vec)[2],dim(shefs_sset)[3])),along=2)
    shefs_mse_in = apply(shefs_in,3,as.vector)
  }
  
  ecrps_dt <- data.table(x=1,hefs=hefs_mse_in,syn=shefs_mse_in)
  
  ylm[l] = sort(as.vector(hefs_mse_vec[,lds[l,]]))[round(0.9*length(as.vector(hefs_mse_vec[,lds[l,]])))]
  
  df <- melt(ecrps_dt,id='x')
  
  df$x <- rep(rep(lds[l,],each=dim(shefs_ecrps_in)[1]),2)
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))
  
  xlbs = paste('ld',lds)
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=clrs[[9]],'shefs'=clrs[[3]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    scale_x_discrete(breaks=lds,labels=xlbs,name='')+
    scale_y_continuous(name = 'MSE')+
    coord_cartesian(ylim = c(0,ylm[l]))+
    #labs(title=paste('ld',lds[l]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,1.1),
            legend.direction = 'horizontal',
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none')}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    guides(fill = guide_legend(byrow = F))
  
  
  ecrps_plts[[l]] <- ecrps_bplt}

#mat<-matrix(1:length(lds),ncol=length(lds))
ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=dim(lds)[1],ncol=1,top='')
#ecrps_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_mse-bplots_pcntile=',disp_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),ecrps_gplot,dpi=320,width=4.5,height=6,unit='in')


#########################################stdev boxplot##############################################
lds = rbind(1:5,6:10,c(11:14,14))
ecrps_plts = vector('list',dim(lds)[1])
show_x = T
ylm = c(35,65,75)
#ylm = c(0.5,0.5,0.5)

for(l in 1:dim(lds)[1]){
  shefs_sdev_in = apply(shefs_sdev_vec[,,lds[l,]],3,as.vector)
  hefs_sdev_in = hefs_sdev_vec[,lds[l,]]
  
  if(dim(shefs_sdev_vec)[2]>dim(hefs_sdev_in)[1]){
    hefs_sdev_in <- rbind(hefs_sdev_in,array(NA,c(dim(shefs_sdev_vec)[2]-dim(hefs_sdev_in)[1],dim(hefs_sdev_in)[2])))
  }
  
  if(dim(shefs_sdev_vec)[2]<dim(hefs_sdev_in)[1]){
    shefs_sset = shefs_sdev_vec[,,lds[l,]]
    shefs_in = abind(shefs_sset,array(NA,c(dim(shefs_sset)[1],dim(hefs_sdev_in)[1] - dim(shefs_sdev_vec)[2],dim(shefs_sset)[3])),along=2)
    shefs_sdev_in = apply(shefs_in,3,as.vector)
  }
  
  ecrps_dt <- data.table(x=1,hefs=hefs_sdev_in,syn=shefs_sdev_in)
  
  ylm[l] = sort(as.vector(hefs_sdev_vec[,lds[l,]]))[round(0.9*length(as.vector(hefs_sdev_vec[,lds[l,]])))]
  
  df <- melt(ecrps_dt,id='x')
  
  df$x <- rep(rep(lds[l,],each=dim(shefs_ecrps_in)[1]),2)
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))
  
  xlbs = paste('ld',lds)
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=clrs[[9]],'shefs'=clrs[[3]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    scale_x_discrete(breaks=lds,labels=xlbs,name='')+
    scale_y_continuous(name = 'StDev')+
    coord_cartesian(ylim = c(0,ylm[l]))+
    #labs(title=paste('ld',lds[l]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,1.1),
            legend.direction = 'horizontal',
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none')}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    guides(fill = guide_legend(byrow = F))
  
  
  ecrps_plts[[l]] <- ecrps_bplt}

#mat<-matrix(1:length(lds),ncol=length(lds))
ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=dim(lds)[1],ncol=1,top='')
#ecrps_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_sdev-bplots_pcntile=',disp_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),ecrps_gplot,dpi=320,width=4.5,height=6,unit='in')


#########################################PBias boxplot##############################################
lds = rbind(1:5,6:10,c(11:14,14))
ecrps_plts = vector('list',dim(lds)[1])
show_x = T
ylm = c(2,2,2)
#ylm = c(1.5,1.5,1.5)

for(l in 1:dim(lds)[1]){
  shefs_pbias_in = apply(shefs_pbias_vec[,,lds[l,]],3,as.vector)
  hefs_pbias_in = hefs_pbias_vec[,lds[l,]]
  
  if(dim(shefs_pbias_vec)[2]>dim(hefs_pbias_in)[1]){
    hefs_pbias_in <- rbind(hefs_pbias_in,array(NA,c(dim(shefs_pbias_vec)[2]-dim(hefs_pbias_in)[1],dim(hefs_pbias_in)[2])))
  }
  
  if(dim(shefs_pbias_vec)[2]<dim(hefs_pbias_in)[1]){
    shefs_sset = shefs_pbias_vec[,,lds[l,]]
    shefs_in = abind(shefs_sset,array(NA,c(dim(shefs_sset)[1],dim(hefs_pbias_in)[1] - dim(shefs_pbias_vec)[2],dim(shefs_sset)[3])),along=2)
    shefs_pbias_in = apply(shefs_in,3,as.vector)
  }
  
  ecrps_dt <- data.table(x=1,hefs=hefs_pbias_in,syn=shefs_pbias_in)
  
  ylm[l] = sort(as.vector(hefs_pbias_vec[,lds[l,]]))[round(0.9*length(as.vector(hefs_pbias_vec[,lds[l,]])))]
  
  df <- melt(ecrps_dt,id='x')
  
  df$x <- rep(rep(lds[l,],each=dim(shefs_ecrps_in)[1]),2)
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))
  
  xlbs = paste('ld',lds)
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=clrs[[9]],'shefs'=clrs[[3]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    scale_x_discrete(breaks=lds,labels=xlbs,name='')+
    scale_y_continuous(name = 'Pbias')+
    coord_cartesian(ylim = c(0,ylm[l]))+
    #labs(title=paste('ld',lds[l]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,1.1),
            legend.direction = 'horizontal',
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none')}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    guides(fill = guide_legend(byrow = F))
  
  
  ecrps_plts[[l]] <- ecrps_bplt}

#mat<-matrix(1:length(lds),ncol=length(lds))
ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=dim(lds)[1],ncol=1,top='')
#ecrps_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_pbias-bplots_pcntile=',disp_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),ecrps_gplot,dpi=320,width=4.5,height=6,unit='in')


#########################################CRPS plots for specific events##############################################
neval_evts = 10
eval_evts <- declust_evts_extract_xhc(obs_key_xhc,neval_evts,sep=15,max_lds = 14,ixx_hefs86,ixx_hefs,ixx_eval_xhc)

ecrps_plts = vector('list',neval_evts)

show_x = rep(c(F,T),each=5)
show_y = rep(c(T,F,F,F,F),2)

for(l in 1:neval_evts){
  set.seed(1)
  ecrps_dt <- data.table(x=1,syn=shefs_ecrps_pk[,l,])
  ecrps_d2t <- data.table(x=1:leads-.25,hefs=shefs_ecrps_pk[sample(1:dim(shefs_ecrps_pk)[1],1),l,])
  
  df <- melt(ecrps_dt,id='x')
  
  df$x <- rep(1:leads,each=dim(shefs_ecrps_pk)[1])
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))
  
  xlbs = paste(1:leads)
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.35,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('shefs'=clrs[[3]]), labels=c(TeX('$sHEFS$')))+
    scale_x_discrete(breaks=1:leads,labels=xlbs,name='leads')+
    {if(show_y[l]==T)
        scale_y_continuous(name = 'CRPS')}+
    {if(show_y[l]==F)
      scale_y_continuous(name = '')}+
    geom_point(aes(x=x,y=hefs,color='hefs'),data=ecrps_d2t,shape=18,size=2)+
    scale_color_manual(name='',breaks='hefs',labels=c(TeX('$sHEFS$ $samp$')),values=c(clrs[[1]]))+
    coord_cartesian(xlim=c(0,leads),ylim = c(0.95*min(df$value),1.05*max(df$value)))+
    annotate('text',x=12,y=1.1*min(df$value),label=ixx_eval_xhc[eval_evts[l]],size=4)+
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=8,angle=90),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.5,1.1),
            legend.box = 'horizontal',
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none')}+
    {if(show_x[l]==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    {if(show_y[l]==F)
      theme(axis.title.y=element_blank())}
  
  ecrps_plts[[l]] <- ecrps_bplt}

mat<-matrix(1:neval_evts,ncol=5,byrow=T)
ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=2,ncol=5,layout_matrix=mat,top='')
#ecrps_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_ecrps-bplots_top',neval_evts,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),ecrps_gplot,dpi=320,width=12,height=5,unit='in')

#########################################CRPS-SS plots for specific events##############################################
neval_evts = 10
eval_evts <- declust_evts_extract_xhc(obs_key_xhc,neval_evts,sep=15,max_lds = 14,ixx_hefs86,ixx_hefs,ixx_eval_xhc)

ecrps_plts = vector('list',neval_evts)

show_x = rep(c(F,T),each=5)
show_y = rep(c(T,F,F,F,F),2)

for(l in 1:neval_evts){
  set.seed(1)
  ecrps_dt <- data.table(x=1,syn=shefs_ecrps_ss_pk[,l,])
  ecrps_d2t <- data.table(x=1:leads-.25,hefs=shefs_ecrps_ss_pk[sample(1:dim(shefs_ecrps_ss_pk)[1],1),l,])
  
  df <- melt(ecrps_dt,id='x')
  
  df$x <- rep(1:leads,each=dim(shefs_ecrps_ss_pk)[1])
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))
  
  xlbs = paste(1:leads)
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.35,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('shefs'=clrs[[3]]), labels=c(TeX('$sHEFS$')))+
    scale_x_discrete(breaks=1:leads,labels=xlbs,name='leads')+
    {if(show_y[l]==T)
      scale_y_continuous(name = 'CRPS-SS')}+
    {if(show_y[l]==F)
      scale_y_continuous(name = '')}+
    geom_point(aes(x=x,y=hefs,color='hefs'),data=ecrps_d2t,shape=18,size=2)+
    scale_color_manual(name='',breaks='hefs',labels=c(TeX('$sHEFS$ $samp$')),values=c(clrs[[1]]))+
    coord_cartesian(xlim=c(0,leads+0.5),ylim = c(-0.1,1),expand=F)+
    annotate('text',x=12,y=0.9,label=ixx_eval_xhc[eval_evts[l]],size=4)+
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=8,angle=90,vjust=0.5),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.5,1.1),
            legend.box = 'horizontal',
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none')}+
    {if(show_x[l]==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    {if(show_y[l]==F)
      theme(axis.title.y=element_blank())}
  
  ecrps_plts[[l]] <- ecrps_bplt}

mat<-matrix(1:neval_evts,ncol=5,byrow=T)
ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=2,ncol=5,layout_matrix=mat,top='')
#ecrps_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_ecrps-ss-bplots_top',neval_evts,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),ecrps_gplot,dpi=320,width=12,height=5,unit='in')

#///////////////////////////////////////////////////////////////////////////////////////////////////////////

#########################################cumulative rank histogram plot##############################################
#calculations
lds = 1:leads
nsamp = dim(shefs_rank_vec)[1]
show_x = T

hefs_cumul_frac <- array(NA,c(n_ens+2,length(lds)))
shefs_cumul_frac <- array(NA,c(nsamp,n_ens+2,length(lds)))

for(i in 1:length(lds)){
  hefs_count <- hist(hefs_rank_vec[,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
  hefs_cumul_frac[,i] <- c(0,roll_sum(hefs_count$counts)/length(hefs_rank_vec[,i]))
  for(s in 1:nsamp){
    shefs_count<-hist(shefs_rank_vec[s,,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
    shefs_cumul_frac[s,,i]<-c(0,roll_sum(shefs_count$counts)/length(shefs_rank_vec[s,,i]))
  }
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Line plot
#plot
rank_plts <- vector('list',length(lds))
show_x = c(rep(F,10),rep(T,5))
show_y = rep(c(T,rep(F,4)),3)

for(i in 1:length(lds)){
  rank_dt <- data.table(x=0:(n_ens+1),syn=t(shefs_cumul_frac[,,i]),hefs=hefs_cumul_frac[,i])
  
  #convert to long form
  df <- melt(rank_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  df$grp <- factor(x=as.numeric(id_vec),levels=1:2,labels=c("hefs","shefs"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=grp,linewidth=grp,alpha=grp))+
    #scale_color_manual(values=c('hefs'='gray40','shefs'=clrs[[4]]))+
    scale_linewidth_manual(values=c('hefs'=1.5,'shefs'=.75), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    scale_alpha_manual(values=c('hefs'=1,'shefs'=.15), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    scale_color_manual(values=c('hefs'=clrs[[9]],'shefs'=clrs[[3]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS$')))+
    labs(x='ensemble rank',y='cumulative fraction')+
    coord_cartesian(xlim=c(0,n_ens+2),ylim=c(-0.01,1),expand=F)+
    #labs(title=paste(lds[i],'hr'),tag=paste(letters[i],')',sep=''))+
    labs(title=paste('ld',lds[i]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1),
          panel.grid.major.x = element_line(size=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(i==15)
      theme(legend.position = 'inside',
            legend.position.inside = c(.3,.8),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=14))}+
    {if(i!=15)
      theme(legend.position = 'none')}+
    {if(show_y[i]==F)
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x[i]==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  
  rank_plts[[i]] <- hefs_ens_plt}


mat<-matrix(1:15,ncol=5,byrow=T)
rank_gplot<-marrangeGrob(rank_plts,nrow=3,ncol=5,top='',layout_matrix = mat,widths=c(1.35,rep(1,4)))
#rank_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_rankhist-lineplot_pcntile=',disp_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),rank_gplot,dpi=320,width=9,height=7,unit='in')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#bounds plot
shefs_rank_upr<-apply(shefs_cumul_frac,c(2,3),max)
shefs_rank_lwr<-apply(shefs_cumul_frac,c(2,3),min)

shefs_plts <- vector('list',length(lds))

for(i in 1:length(lds)){
  rank_dt <- data.table(x=0:(n_ens+1),shefs_upr=shefs_rank_upr[,i],shefs_lwr=shefs_rank_lwr[,i],hefs=hefs_cumul_frac[,i])
  
  hefs_ens_plt<-ggplot(rank_dt)+theme_minimal()+
    geom_ribbon(mapping=aes(x=x,ymax=shefs_upr,ymin=shefs_lwr,fill='shefs'),alpha=0.25)+
    geom_line(mapping=aes(x=x,y=hefs,color='hefs'),linewidth=1)+
    scale_color_manual(values=c('hefs'=clrs[[9]]), labels=c(TeX('$HEFS$')))+
    scale_fill_manual(values=c('shefs'=clrs[[3]]), labels=c(TeX('$sHEFS$')))+
    labs(x='ensemble rank',y='cumulative fraction')+
    coord_cartesian(xlim=c(0,n_ens+2),ylim=c(-0.01,1),expand=F)+
    labs(title=paste('ld',lds[i]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1),
          panel.grid.major.x = element_line(size=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(i==15)
      theme(legend.position = 'inside',
            legend.position.inside = c(.3,.8),
            legend.box = 'vertical',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=14))}+
    {if(i!=15)
      theme(legend.position = 'none')}+
    {if(show_y[i]==F)
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x[i]==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    guides(shape = guide_legend(order=2),col=guide_legend(order=1))
  
  
  shefs_plts[[i]] <- hefs_ens_plt}

mat<-matrix(1:15,ncol=5,byrow=T)
shefs_gplot<-marrangeGrob(shefs_plts,nrow=3,ncol=5,layout_matrix = mat,top='',widths=c(1.35,rep(1,4)))
#shefs_gplot

ggsave(paste(path_out,loc,'-',disp_site,'_rankhist-shadeplot_pcntile=',disp_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_xhc.png',sep=''),shefs_gplot,dpi=320,width=9,height=7,unit='in')


#########################################################END########################################################
