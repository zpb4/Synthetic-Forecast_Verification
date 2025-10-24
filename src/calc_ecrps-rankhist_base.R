#Script to calculate verification statistics for the HINDCAST period

#NOTE: This script can be time-consuming to run and requires a lot of RAM; best to run on HPC if possible
print(paste('calc start',Sys.time()))

#set root directory
setwd('z:/Synthetic-Forecast_Verification/')

#Load packages
library(lubridate)

#Primary modifiable input parameters
#///////////////////////////////////////////////////////////////////////////////////////////////////////////
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

#calculation setup
disp_pcnt <- 0.99       #what percentile of the data to calculate statistics against? 

#directory for synthetic forecasts
path = paste('../Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

#path to output data; default is to output to a 'data' subrepo in the specified root directory in Line 7 above
path_out = './data'

#//////////////////////////////////////////////////////////////////////////////////

if (!dir.exists(path_out)) {
  dir.create(path_out,recursive=T)
}

#load data
if(has_86==T){
load(paste(path,'out/',loc,'/data_prep_rdata86.RData',sep=''))
cur_site <- which(site_names==disp_site)
idx_site <- which(site_names==opt_site)
ixx_hefs86 <- ixx_hefs
hefs_fwd_86 <- hefs_forward[cur_site,,,]}

load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))
cur_site <- which(site_names==disp_site)
idx_site <- which(site_names==opt_site)

syn_hefs_forward <- readRDS(paste(path,'out/',loc,'/syn_hefs_forward_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',opt_site,'_',cal_val_setup,'.rds',sep=''))
shefs_fwd <- syn_hefs_forward[,cur_site,,,]
hefs_fwd_sset <- hefs_forward[cur_site,,,]
obs_fwd <- obs_forward_all_leads[cur_site,,]

rm(syn_hefs_forward,hefs_forward,hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,obs_forward_all_leads_hind)
gc()

hefs_fwd <- shefs_fwd[1,,,]
hefs_idx <- ixx_obs_forward%in%ixx_hefs
hefs_fwd[,hefs_idx,] <- hefs_fwd_sset[,ixx_hefs%in%ixx_obs_forward,]

if(has_86==T){
hefs86_idx <- ixx_obs_forward%in%ixx_hefs86
hefs_fwd[,hefs86_idx,] <- hefs_fwd_86
ixx_keep = ixx_obs_forward%in%c(ixx_hefs86,ixx_hefs)}

if(has_86==F){
  ixx_keep = ixx_obs_forward%in%c(ixx_hefs)}

source('./src/forecast_verification_functions.R')

ixx_obs_forward_wy <- wy_fun(ixx_obs_forward)

#calculate climatology array
climo_farray <- climo_forecast(ixx_obs_forward,hefs_fwd,obs_fwd)

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#calculate date indices
hefs_eval <- hefs_fwd[,ixx_keep,]
shefs_eval<- shefs_fwd[,,ixx_keep,]
climo_eval <- climo_farray[,ixx_keep,]
obs_in <- obs[ixx_obs%in%ixx_obs_forward,cur_site]
obs_samp <- obs[ixx_obs%in%ixx_obs_forward,idx_site]
obs_eval <- obs_in[ixx_keep]
obs_key <- obs_samp[ixx_keep]
rm(hefs_fwd,shefs_fwd,obs,obs_in,obs_samp);gc()

ixx_eval <- ixx_obs_forward[ixx_keep]

saveRDS(obs_eval,paste('./data/',loc,'-',disp_site,'_obs-eval.rds',sep=''))
saveRDS(obs_key,paste('./data/',loc,'-',opt_site,'_obs-key.rds',sep=''))
saveRDS(ixx_eval,paste('./data/',loc,'-',disp_site,'_ixx-eval.rds',sep=''))

#calculate specified subset of the data for calculation
lds<-1:leads  
n_evts<- round((1-disp_pcnt) * length(obs_eval))

obs_date_loc <- order(obs_key,decreasing=TRUE)[1:n_evts]  #index for maximum observation
rmv_idx <- which(obs_date_loc<=leads)
if(length(rmv_idx)>0){
  pool <- order(obs_key,decreasing=TRUE)[(n_evts+1):(n_evts+100)]
  pool <- pool[pool>leads]
  obs_date_loc <- obs_date_loc[-c(rmv_idx)]
  obs_date_loc <- c(obs_date_loc,pool[1:length(rmv_idx)])}

obs_events <- obs_eval[obs_date_loc]
samps <- dim(shefs_eval)[1]

#calculate eCRPS and rank histogram data
hefs_ecrps_vec<-array(NA,c(n_evts,length(lds)))
shefs_ecrps_vec<-array(NA,c(samps,n_evts,length(lds)))
climo_ecrps_vec<-array(NA,c(n_evts,length(lds)))

hefs_rank_vec <- array(NA,c(n_evts,length(lds)))
shefs_rank_vec <- array(NA,c(samps,n_evts,length(lds)))

#HEFS
for(ld in 1:length(lds)){
  for(i in 1:n_evts){
    hefs_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    HEFS <- hefs_eval[,hefs_idx,lds[ld]]
    CLIM <- climo_eval[,hefs_idx,lds[ld]]
    hefs_ecrps_vec[i,ld] <- eCRPS(HEFS,obs_events[i])
    climo_ecrps_vec[i,ld] <- eCRPS(CLIM,obs_events[i])
    hefs_rank_vec[i,ld] <- ens_rank(HEFS,obs_events[i])
  }
}

saveRDS(hefs_ecrps_vec,paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-ecrps-vec.rds',sep=''))
saveRDS(climo_ecrps_vec,paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_climo-ecrps-vec.rds',sep=''))
saveRDS(hefs_rank_vec,paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-rank-vec.rds',sep=''))

#sHEFS
for(s in 1:samps){
  for(ld in 1:length(lds)){
    for(i in 1:n_evts){
      syn_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
      SYN_HEFS <- shefs_eval[s,,syn_idx,lds[ld]]
      shefs_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,obs_events[i])
      shefs_rank_vec[s,i,ld] <- ens_rank(SYN_HEFS,obs_events[i])}

  }
  print(paste(s))
}
  
saveRDS(shefs_ecrps_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-vec.rds',sep=''))
saveRDS(shefs_rank_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-rank-vec.rds',sep=''))

#calculate eCRPS skill score
hefs_ecrps_ss <- 1 - hefs_ecrps_vec/climo_ecrps_vec
shefs_ecrps_ss <- array(NA,dim(shefs_ecrps_vec))
for(i in 1:dim(shefs_ecrps_vec)[1]){
  shefs_ecrps_ss[i,,] <- 1 - shefs_ecrps_vec[i,,]/climo_ecrps_vec
}

saveRDS(hefs_ecrps_ss,paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-ecrps-ss.rds',sep=''))
saveRDS(shefs_ecrps_ss,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-ss.rds',sep=''))

#calc mse and sdev stats
hefs_mse_vec<-array(NA,c(n_evts,length(lds)))
shefs_mse_vec<-array(NA,c(samps,n_evts,length(lds)))

hefs_sdev_vec <- array(NA,c(n_evts,length(lds)))
shefs_sdev_vec <- array(NA,c(samps,n_evts,length(lds)))

hefs_pbias_vec <- array(NA,c(n_evts,length(lds)))
shefs_pbias_vec <- array(NA,c(samps,n_evts,length(lds)))

#HEFS
for(ld in 1:length(lds)){
  for(i in 1:n_evts){
    hefs_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    HEFS <- hefs_eval[,hefs_idx,lds[ld]]
    hefs_mse_vec[i,ld] <- (mean(HEFS)-obs_events[i])^2
    hefs_sdev_vec[i,ld] <- sd(HEFS)
    hefs_pbias_vec[i,ld] <- mean(HEFS) / obs_events[i]
  }
}

saveRDS(hefs_mse_vec,paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-mse-vec.rds',sep=''))
saveRDS(hefs_sdev_vec,paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-sdev-vec.rds',sep=''))
saveRDS(hefs_pbias_vec,paste('./data/',loc,'-',disp_site,'_pct=',disp_pcnt,'_hefs-pbias-vec.rds',sep=''))

#sHEFS
for(s in 1:samps){
  for(ld in 1:length(lds)){
    for(i in 1:n_evts){
      syn_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
      SYN_HEFS <- shefs_eval[s,,syn_idx,lds[ld]]
      shefs_mse_vec[s,i,ld] <- (mean(SYN_HEFS)-obs_events[i])^2
      shefs_sdev_vec[s,i,ld] <- sd(SYN_HEFS)
      shefs_pbias_vec[s,i,ld] <- mean(SYN_HEFS) / obs_events[i]
      }
  }
  print(paste(s))
}

saveRDS(shefs_mse_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-mse-vec.rds',sep=''))
saveRDS(shefs_sdev_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-sdev-vec.rds',sep=''))
saveRDS(shefs_pbias_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-pbias-vec.rds',sep=''))

#calculate eCRPS for top 10 events
neval_evts = 10
sep = 15

eval_evts <- declust_evts_extract(obs_key,neval_evts,sep,leads)
obs_evts <- obs_eval[eval_evts]

hefs_ecrps_pk <- array(NA,c(neval_evts,leads))
climo_ecrps_pk <- array(NA,c(neval_evts,leads))
shefs_ecrps_pk <- array(NA,c(dim(shefs_ecrps_vec)[1],neval_evts,leads))

#HEFS
for(ld in 1:length(lds)){
  for(i in 1:neval_evts){
    hefs_idx <- eval_evts[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    HEFS <- hefs_eval[,hefs_idx,lds[ld]]
    CLIM <- climo_eval[,hefs_idx,lds[ld]]
    hefs_ecrps_pk[i,ld] <- eCRPS(HEFS,obs_events[i])
    climo_ecrps_pk[i,ld] <- eCRPS(CLIM,obs_events[i])
  }
}

#sHEFS
for(s in 1:samps){
  for(ld in 1:length(lds)){
    for(i in 1:neval_evts){
      syn_idx <- eval_evts[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
      SYN_HEFS <- shefs_eval[s,,syn_idx,lds[ld]]
      shefs_ecrps_pk[s,i,ld] <- eCRPS(SYN_HEFS,obs_events[i])}
  }
}

saveRDS(hefs_ecrps_pk,paste('./data/',loc,'-',disp_site,'_hefs-ecrps-peak10.rds',sep=''))
saveRDS(shefs_ecrps_pk,paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-peak10.rds',sep=''))

#calc CRPS-SS
hefs_ecrps_ss_pk <- 1 - hefs_ecrps_pk/climo_ecrps_pk
shefs_ecrps_ss_pk <- array(NA,dim(shefs_ecrps_pk))
for(i in 1:dim(shefs_ecrps_vec)[1]){
  shefs_ecrps_ss_pk[i,,] <- 1 - shefs_ecrps_pk[i,,]/climo_ecrps_pk
}

saveRDS(hefs_ecrps_ss_pk,paste('./data/',loc,'-',disp_site,'_hefs-ecrps-ss-peak10.rds',sep=''))
saveRDS(shefs_ecrps_ss_pk,paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-ss-peak10.rds',sep=''))

print(paste('calc end',Sys.time()))

rm(list=ls());gc()


#####################################################END###################################################