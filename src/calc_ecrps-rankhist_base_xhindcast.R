#Script to calculate verification statistics for the x-HINDCAST period (outside of HINDCAST period [+ or -])

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
ixx_hefs86 <- ixx_hefs}

load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))
cur_site <- which(site_names==disp_site)
idx_site <- which(site_names==opt_site)

syn_hefs_forward <- readRDS(paste(path,'out/',loc,'/syn_hefs_forward_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',opt_site,'_',cal_val_setup,'.rds',sep=''))
shefs_fwd <- syn_hefs_forward[,cur_site,,,]
obs_fwd <- obs_forward_all_leads[cur_site,,]

rm(syn_hefs_forward,hefs_forward,hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,obs_forward_all_leads_hind)
gc()

if(has_86==T){
  ixx_keep = !ixx_obs_forward%in%c(ixx_hefs,ixx_hefs86)
  nevt_ref = c(ixx_hefs,ixx_hefs86)}

if(has_86==F){
  ixx_keep = !ixx_obs_forward%in%c(ixx_hefs)
  nevt_ref = ixx_hefs}

source('./src/forecast_verification_functions.R')

#calculate climatology array
climo_farray <- climo_forecast(ixx_obs_forward,shefs_fwd[1,,,],obs_fwd)

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#calculate date indices
shefs_eval<- shefs_fwd[,,ixx_keep,]
climo_eval <- climo_farray[,ixx_keep,]
obs_in <- obs[ixx_obs%in%ixx_obs_forward,cur_site]
obs_samp <- obs[ixx_obs%in%ixx_obs_forward,idx_site]
obs_eval <- obs_in[ixx_keep]
obs_key <- obs_samp[ixx_keep]
rm(shefs_fwd,obs,obs_in,obs_samp);gc()

ixx_eval <- ixx_obs_forward[ixx_keep]

saveRDS(obs_eval,paste('./data/',loc,'-',disp_site,'_obs-eval_xhc.rds',sep=''))
saveRDS(obs_key,paste('./data/',loc,'-',opt_site,'_obs-key_xhc.rds',sep=''))
saveRDS(ixx_eval,paste('./data/',loc,'-',disp_site,'_ixx-eval_xhc.rds',sep=''))

hefs_len <- length(unique(ixx_hefs$year))
syn_len <- length(unique(ixx_eval$year))

#calculate indices for events based on specified percentile
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


#resample x-hindcast period at lengths equal to HEFS dataset
#NOTE: If x-hindcast period is longer than HEFS, this is done without replacement from the x-hindcast record
#If x-hindcast period is shorter than HEFS, this is done with replacement from x-hindcast record
obs_events_mat <- matrix(rep(obs_events,samps),ncol=samps,byrow=F)
obs_dloc_mat <- matrix(rep(obs_date_loc,samps),ncol=samps,byrow=F)

n_evts<- round((1-disp_pcnt) * length(ixx_hefs))
obs_events_mat <- array(NA,c(n_evts,samps))
obs_dloc_mat <- array(NA,c(n_evts,samps))
for(i in 1:samps){
  if(syn_len > hefs_len){
    yr_samp <- sample(unique(ixx_eval$year),hefs_len,replace=F)}
  if(syn_len <= hefs_len){
    yr_samp <- sample(unique(ixx_eval$year),hefs_len,replace=T)}
  obs_sset <- obs_key[ixx_eval$year%in%yr_samp]
  yr_idx <- c(1:length(obs_key))[ixx_eval$year%in%yr_samp]
  obs_idxs_sset <- order(obs_sset,decreasing=TRUE)[1:n_evts]  #index for maximum observation
  obs_dates <- ixx_eval[yr_idx][obs_idxs_sset]
  obs_events <- obs_eval[ixx_eval%in%obs_dates]
  obs_date_loc <- c(1:length(obs_key))[ixx_eval%in%obs_dates] 
  rmv_idx <- which(obs_date_loc<=leads)
  obs_dates <- as.character(obs_dates)
  #remove any dates that are within the lead window of the end of 1986 period to ensure no cutoffs
  idx_close_86 <- seq(as.Date(tail(ixx_hefs86,1)+(60*60*24)),as.Date(tail(ixx_hefs86,1)+(60*60*24*leads)),by='day')
  if(any(obs_dates%in%idx_close_86==T)){
    rmv_idx <- c(rmv_idx,(1:length(obs_dates))[obs_dates%in%idx_close_86])}
  #remove any dates that are within the lead window of the end of the HEFS period to ensure no cutoffs
  idx_close_hefsend <- seq(as.Date(tail(ixx_hefs,1)+(60*60*24)),as.Date(tail(ixx_hefs,1)+(60*60*24*leads)),by='day')
  if(any(obs_dates%in%idx_close_hefsend==T)){
    rmv_idx <- c(rmv_idx,(1:length(obs_dates))[obs_dates%in%idx_close_hefsend])}
  if(length(rmv_idx)>0){
    pool <- order(obs_key,decreasing=TRUE)[(n_evts+1):(n_evts+100)]
    pool <- pool[pool>leads]
    obs_date_loc <- obs_date_loc[-c(rmv_idx)]
    obs_date_loc <- c(obs_date_loc,pool[1:length(rmv_idx)])}
    
  obs_events_mat[,i] <- obs_events
  obs_dloc_mat[,i] <- obs_date_loc
}


#calculate the eCRPS statistic and rank histogram data for the specified subset
shefs_ecrps_vec<-array(NA,c(samps,n_evts,length(lds)))
climo_ecrps_vec<-array(NA,c(samps,n_evts,length(lds)))
shefs_rank_vec <- array(NA,c(samps,n_evts,length(lds)))

for(s in 1:samps){
  for(ld in 1:length(lds)){
    for(i in 1:n_evts){
      syn_idx <- obs_dloc_mat[i,s]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
      SYN_HEFS <- shefs_eval[s,,syn_idx,lds[ld]]
      CLIM <- climo_eval[,syn_idx,lds[ld]]
      shefs_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,obs_events_mat[i,s])
      climo_ecrps_vec[s,i,ld] <- eCRPS(CLIM,obs_events_mat[i,s])
      shefs_rank_vec[s,i,ld] <- ens_rank(SYN_HEFS,obs_events_mat[i,s])}

  }
  print(paste(s))
}
  
saveRDS(shefs_ecrps_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-vec_xhc.rds',sep=''))
saveRDS(shefs_ecrps_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_climo-ecrps-vec_xhc.rds',sep=''))
saveRDS(shefs_rank_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-rank-vec_xhc.rds',sep=''))

#calc eCRPS skill score
shefs_ecrps_ss <- 1 - shefs_ecrps_vec/climo_ecrps_vec

saveRDS(shefs_ecrps_ss,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-ss_xhc.rds',sep=''))

#calc mse and sdev stats
shefs_mse_vec<-array(NA,c(samps,n_evts,length(lds)))
shefs_sdev_vec <- array(NA,c(samps,n_evts,length(lds)))
shefs_pbias_vec <- array(NA,c(samps,n_evts,length(lds)))

for(s in 1:samps){
  for(ld in 1:length(lds)){
    for(i in 1:n_evts){
      syn_idx <- obs_dloc_mat[i,s]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
      SYN_HEFS <- shefs_eval[s,,syn_idx,lds[ld]]
      shefs_mse_vec[s,i,ld] <- (mean(SYN_HEFS)-obs_events_mat[i,s])^2
      shefs_sdev_vec[s,i,ld] <- sd(SYN_HEFS)
      shefs_pbias_vec[s,i,ld] <- mean(SYN_HEFS) / obs_events_mat[i,s]
      }
  }
  print(paste(s))
}

saveRDS(shefs_mse_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-mse-vec_xhc.rds',sep=''))
saveRDS(shefs_sdev_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-sdev-vec_xhc.rds',sep=''))
saveRDS(shefs_pbias_vec,paste('./data/',loc,'-',disp_site,'_pcntile=',disp_pcnt,'_setup=',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-pbias-vec_xhc.rds',sep=''))

#calculate eCRPS for top 10 events
neval_evts = 10
sep = 15

eval_evts <- declust_evts_extract_xhc(obs_key,neval_evts,sep,leads,ixx_hefs86,ixx_hefs,ixx_eval)
obs_evts <- obs_eval[eval_evts]

shefs_ecrps_pk <- array(NA,c(dim(shefs_ecrps_vec)[1],neval_evts,leads))
climo_ecrps_pk <- array(NA,c(dim(shefs_ecrps_vec)[1],neval_evts,leads))

for(s in 1:samps){
  for(ld in 1:length(lds)){
    for(i in 1:neval_evts){
      syn_idx <- eval_evts[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
      SYN_HEFS <- shefs_eval[s,,syn_idx,lds[ld]]
      CLIM <- climo_eval[,syn_idx,lds[ld]]
      shefs_ecrps_pk[s,i,ld] <- eCRPS(SYN_HEFS,obs_evts[i])
      climo_ecrps_pk[s,i,ld] <- eCRPS(CLIM,obs_evts[i])}
  }
}

saveRDS(shefs_ecrps_pk,paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-peak10_xhc.rds',sep=''))
saveRDS(climo_ecrps_pk,paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_climo-ecrps-peak10_xhc.rds',sep=''))

#calculate eCRPS skill score
shefs_ecrps_ss_pk <- 1 - shefs_ecrps_pk/climo_ecrps_pk

saveRDS(shefs_ecrps_ss_pk,paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-ss-peak10_xhc.rds',sep=''))


print(paste('calc end',Sys.time()))

rm(list=ls());gc()

#####################################################END###################################################