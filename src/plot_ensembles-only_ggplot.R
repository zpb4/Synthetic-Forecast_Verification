

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
shefs_fwd <- syn_hefs_forward[,cur_site,,,]
hefs_fwd_sset <- hefs_forward[cur_site,,,]

rm(syn_hefs_forward,hefs_forward,hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,obs_forward_all_leads_hind)
gc()

hefs_fwd <- shefs_fwd[1,,,]
hefs_idx <- ixx_obs_forward%in%ixx_hefs
hefs_fwd[,hefs_idx,] <- hefs_fwd_sset[,ixx_hefs%in%ixx_obs_forward,]

if(has_86==T){
  hefs86_idx <- ixx_obs_forward%in%ixx_hefs86
  hefs_fwd[,hefs86_idx,] <- hefs_fwd_86}

source('./src/forecast_verification_functions.R')

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
##############Ensemble Plots###############
n_evts = 10
sep = 15

obs_extract <- obs[,idx_site]
obs_extract[!ixx_obs%in%ixx_hefs] <- 0
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


##########################HEFS plot###########################################
lds = c(10,5,3,1)
hefs_plts = vector('list',length(lds))
ylm_scale = 1.5
plt_evts = 5

for(k in 1:plt_evts){
  show_x = F
for(i in 1:length(lds)){
  hefs_plt_dt <- data.table(x=0:leads,hefs=rbind(rep(obs_gen[obs_evt_idx[k]-lds[i]],dim(hefs_fwd)[1]),t(hefs_fwd[,obs_evt_idx[k]-lds[i],])),
                            obs=c(obs_gen[obs_evt_idx[k]-lds[i]],obs_fwd_gen[obs_evt_idx[k]-lds[i],]))
  
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
    geom_vline(xintercept = lds[i],linetype='dotted',linewidth=0.5)+
    {if(lds[i]>5 & i!=1)annotate('text',x=lds[i],y=0.95*ylm,label=obs_dates[k],size=4,hjust=1.15)}+
    {if(lds[i]<=5 & i!=1)annotate('text',x=lds[i],y=0.95*ylm,label=obs_dates[k],size=4,hjust=-0.15)}+
    {if(i==1)annotate('text',x=leads,y=0.95*ylm,label=obs_dates[k],size=4,hjust=1.15)}+
    #labs(title=paste('Lead',lds[i]))+
    #{if(i==1)annotate('text',x=-0.1*leads,y=0.5*ylm,label='HEFS',size=8,hjust=0.5,angle=90,fontface=2)}+
    #annotate('text',x=1,y=850,label='b)',size=6)+
    labs(x='forecast lead (hrs)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,leads),ylim=c(0,ylm),expand=F,clip='off')+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          plot.title=element_text(hjust=0.5,size=20,face='bold'),
          panel.grid.major.y = element_line(size=1),
          panel.grid.major.x = element_line(size=1))+
          #plot.margin=unit(c(0.25,0.25,0.25,0.25),'cm'))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,.9),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
            #plot.margin=unit(c(0.25,0.25,0.25,1.5),'cm'))}+
    {if(i!=1)
      theme(legend.position = 'none',
            axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  hefs_plts[[i]] <- hefs_ens_plt}

mat<-matrix(1:length(lds),ncol=length(lds))
hefs_gplot<-marrangeGrob(hefs_plts,nrow=1,ncol=length(lds),top='')
#hefs_gplot

#ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_production/figs/ind_plots/',loc,'_',disp_site,'_hefs-ens-plot_evt=',obs_rank,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),hefs_gplot,dpi=320,width=9,height=3,unit='in')

##########################sHEFS plot###########################################
samps = 3
shefs_plts = vector('list',length(lds)*samps)
#ylm_scale = 1.5
#show_x = T
show_x = c(rep(F,(samps-1)*length(lds)),rep(T,length(lds)))

set.seed(seed)
samp_idx <- sample(1:10,samps);samp_vec <- rep(samp_idx,each=length(lds))
ldsv <- rep(lds,samps)

for(i in 1:(length(lds)*samps)){
  hefs_plt_dt <- data.table(x=0:leads,hefs=rbind(rep(obs_gen[obs_evt_idx[k]-ldsv[i]],dim(hefs_fwd)[1]),t(shefs_fwd[samp_vec[i],,obs_evt_idx[k]-ldsv[i],])),
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
    #{if(ldsv[i]!=1)annotate('text',x=lds[i]*24,y=0.95*ylm,label=obs_date,size=4,hjust=1.15)}+
    #{if(ldsv[i]==1)annotate('text',x=lds[i]*24,y=0.95*ylm,label=obs_date,size=4,hjust=-0.15)}+
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

#ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_production/figs/ind_plots/',loc,'_',disp_site,'_shefs-ens-plot_evt=',obs_rank,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),shefs_gplot,dpi=320,width=3*length(lds),height=3*samps,unit='in')


comb_mat <- matrix(1:((samps+1)*length(lds)),nrow=samps+1,ncol=length(lds),byrow=T)
comb_gplot<-marrangeGrob(c(hefs_plts,shefs_plts),nrow=samps+1,ncol=length(lds),top='',layout_matrix = comb_mat)
#comb_gplot

ggsave(paste(path_out,loc,'_',disp_site,'_hefs-shefs-ens-plot_setup=',cal_val_setup,'_pct=',opt_pcnt,'_pwr=',obj_pwr,'_strat=',opt_strat,'_evt=',k,'_lds=',str_flatten(lds,collapse='-'),'_seed=',seed,'.png',sep=''),comb_gplot,dpi=320,width=3*length(lds),height=2.5*(samps+1),unit='in')

}
#////////////////////////////////////////////////////////////////////////////


##############################END###############################################################

