
#args = commandArgs(trailingOnly=TRUE)
#print(paste('task #',args[1]))
#idx = as.numeric(args[1])


print(paste('calc start',Sys.time()))

#library(abind)
#library(doParallel)
#parallel::detectCores()
#n.cores <- parallel::detectCores()
#my.cluster<-parallel::makeCluster(n.cores,type = 'FORK',methods=F,useXDR=F)
#my.cluster<-parallel::makeCluster(n.cores,type = 'PSOCK')
#print(my.cluster)
#doParallel::registerDoParallel(cl = my.cluster)
#foreach::getDoParRegistered()

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Data setup
syn_vers = 2
loc = 'SOD'
opt_site = 'SRWC1'
disp_site = 'SRWC1'
opt_pcnt = 0.9901
cal_val_setup = 'cal' # 'cal' '5fold' '5fold-test'
obj_pwr = 0
opt_strat = 'ecrps-dts'
has_86 = F

#plot setup
disp_pcnt <- 0.995

path = paste('../Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

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
  ixx_keep = !ixx_obs_forward%in%c(ixx_hefs,ixx_hefs_86)
  nevt_ref = c(ixx_hefs,ixx_hefs_86)}

if(has_86==F){
  ixx_keep = !ixx_obs_forward%in%c(ixx_hefs)
  nevt_ref = ixx_hefs}

source('./src/forecast_verification_functions.R')

#calculate climatology array
climo_farray <- climo_forecast(ixx_obs_forward,shefs_fwd[1,,,],obs_fwd)

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Ensemble plots
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

###############eCRPS + Rank Histogram####################################
lds<-1:leads  #specify leads (no more than 5 for plotting constraints)
n_evts<- round((1-disp_pcnt) * length(obs_eval))

###eCRPS###
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
  if(length(rmv_idx)>0){
    pool <- order(obs_key,decreasing=TRUE)[(n_evts+1):(n_evts+100)]
    pool <- pool[pool>leads]
    obs_date_loc <- obs_date_loc[-c(rmv_idx)]
    obs_date_loc <- c(obs_date_loc,pool[1:length(rmv_idx)])}
    
  obs_events_mat[,i] <- obs_events
  obs_dloc_mat[,i] <- obs_date_loc
}

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

#calc CRPS-SS
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

eval_evts <- declust_evts_extract(obs_key,neval_evts,sep,leads)
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

#calc CRPS-SS
shefs_ecrps_ss_pk <- 1 - shefs_ecrps_pk/climo_ecrps_pk

saveRDS(shefs_ecrps_ss_pk,paste('./data/',loc,'-',disp_site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_shefs-ecrps-ss-peak10_xhc.rds',sep=''))


print(paste('calc end',Sys.time()))

rm(list=ls());gc()


#####################################################END###################################################