

rm(list=ls())
library(fields)
library(scales)
library(zoo)
#----------------------------------------
syn_vers = 2
loc = 'YRS'
site = 'ORDC1'
opt_pcnt = 0.99
cal_val_setup = '5fold' # 'cal' 'val' 'wet' '5fold'
v1_param = 'a'

#plot setup
disp_pcnt <- 0.99
n_samp <- 10
save_plots <- F  # T to save plots to png files, F to display in R

path = paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))

if(syn_vers==1){
  syn_hefs_forward <- readRDS(paste(path,'out/',loc,'/syn_hefs_forward-',v1_param,'_',site,'_',cal_val_setup,'_plot-ens.rds',sep=''))} #use only a slice of larger sample array
if(syn_vers==2){
  syn_hefs_forward <- readRDS(paste(path,'out/',loc,'/syn_hefs_forward_pcnt=',opt_pcnt,'_',site,'_',cal_val_setup,'_plot-ens.rds',sep=''))} #use only a slice of larger sample array

ixx_gen <- readRDS(paste(path,'out/',loc,'/ixx_gen.rds',sep='')) 

source('./src/forecast_verification_functions.R')


if(cal_val_setup=='cal' | cal_val_setup=='5fold' | cal_val_setup=='5fold-test'){
  leave_out_years = 48:119 + 1900
}

if(cal_val_setup=='val'){
  lv_out_samps = 6          #how many validation years to leave out
  opt_leave_out = readRDS(paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/data/',loc,'/opt_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
  leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
}

if(cal_val_setup=='wet'){
  lv_out_samps = 6          #how many validation years to leave out
  opt_leave_out = readRDS(paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/data/',loc,'/wet_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
  leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
}

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

ixx_hefs_wy <- wy_fun(ixx_hefs)
ixx_obs_wy <- wy_fun(ixx_obs)
ixx_obs_fwd_wy <- wy_fun(ixx_obs_forward)
ixx_gen_wy <- wy_fun(ixx_gen)

hefs_idx <- ixx_hefs_wy$year%in%(leave_out_years-1900)
obs_idx <- ixx_obs_wy$year%in%(leave_out_years-1900)
obss_idx <- ixx_obs[obs_idx]%in%ixx_hefs[hefs_idx]
obs_fwd_idx <- ixx_obs_fwd_wy$year%in%(leave_out_years-1900)
obss_fwd_idx <- ixx_obs_forward[obs_fwd_idx]%in%ixx_hefs[hefs_idx]
gen_idx <- ixx_gen_wy$year%in%(leave_out_years-1900)
genss_idx <- ixx_gen[gen_idx]%in%ixx_hefs[hefs_idx]
#ixx_val <- ixx_hefs[hefs_idx]
ixx_val <- ixx_obs[obs_idx]
##############Ensemble Density Plots###############

cur_site <- which(site_names==site)

obs <- obs[obs_idx,cur_site]
#obs <- obs[obss_idx]
obs_fwd <- obs_forward_all_leads[cur_site,obs_fwd_idx,]
#obs_fwd <- obs_fwd[obss_fwd_idx,]
hefs_forward <- hefs_forward[,,hefs_idx,,drop=F]
syn_hefs_forward <- syn_hefs_forward[,,,gen_idx,,drop=F]
#syn_hefs_forward <- syn_hefs_forward[,,,genss_idx,,drop=F]

obs_rank <- 1   #pick which obs event to plot (1 largest, 2 second largest, etc)
seed<-5
set.seed(seed)
#get date for plot
#obs_date_loc <- match('1997-01-02',as.character(ixx_obs))
obs_date_loc <- which(obs==sort(obs,decreasing=TRUE)[obs_rank])  #index for maximum observation
obs_date <- ixx_val[obs_date_loc]

if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_event-plot_rnk=',obs_rank,'_gray-clean2_seed=',seed,'.png',sep=''),height = 2500,width = 4096,res=300)}
num_plot_rows <- min(n_samp+1,4)
par(mfcol=c(num_plot_rows-1,4),mar=c(3.5,4,0,0),mgp=c(1,.5,0))
#ymax <- .5*max(obs[,cur_site],hefs_forward[cur_site,,,],syn_hefs_forward[,cur_site,,,])
ymax <- 1.75*obs[obs_date_loc]

for (ld in c(1,3,5,10)) {
  plt_idx <- obs_date_loc-ld #ref index for ensemble plots
  
  #1) Forward looking forecast plots
  #hefs
  #plot(0:leads,obs[plt_idx:(plt_idx+leads)],type='l',lwd=3,main="",ylim=c(0,ymax),
       #xlab='',ylab='',axes=F)
  #axis(2,at=seq(0,15,5),labels=F)
  #if(ld == 1){
    #axis(2,at=seq(0,15,5),labels=seq(0,15,5),cex.axis=1.5)
  #}
  #axis(1,at=seq(0,leads,2),labels=seq(0,leads,2),cex.axis=1.5)
  #text(ld+4.5,ymax-1,obs_date,col='darkorange3')
  #if (length(obs_date_loc)>0) {
    #cur_sum <- round(sum(apply(hefs_forward[cur_site,,plt_idx,],2,FUN=mean)))
    #tt <- paste('HEFS ld',ld,'sum',cur_sum)
    #title(tt)
    #for(i in 1:n_ens){
      #lines(0:leads,c(obs[plt_idx],hefs_forward[cur_site,i,plt_idx,]),col='gray80')
    #}
    #lines(0:leads,obs[plt_idx:(plt_idx+leads)],lwd=3)
    #if(ld==1){
    #mtext('Flow (kcfs)',side=2,line=2,cex=1.5)}
    #mtext('Lead (days)',side=1,line=2.5,cex=1.5)
    #abline(v=ld,col='darkorange3',lty=2)
  #}
  
  #synthetic
  samples_to_plot <- sample(1:n_samp,(num_plot_rows-1),replace=FALSE)
  for (samp in samples_to_plot) {
    plot(0:leads,obs[plt_idx:(plt_idx+leads)],type='l',lwd=3,main="",ylim=c(0,ymax),
         xlab='',ylab='',axes=F)
    axis(2,at=seq(0,30,5),labels=F)
    if(ld == 1){
      axis(2,at=seq(0,30,5),labels=seq(0,30,5),cex.axis=1.5)
    }
    axis(1,at=seq(0,leads,2),labels=seq(0,leads,2),cex.axis=1.5)
    #text(ld+4.5,ymax-1,obs_date,col='darkorange3')
    if (length(obs_date_loc)>0) {
      cur_sum <- round(sum(apply(syn_hefs_forward[samp,cur_site,,plt_idx,],2,FUN=mean)))
      tt <- paste('syn',samp,'ld',ld,'sum',cur_sum)
      #title(tt)
      for(i in 1:n_ens){
        lines(0:leads,c(obs[plt_idx],syn_hefs_forward[samp,cur_site,i,plt_idx,]),col='gray80')
      }
      lines(0:leads,obs[plt_idx:(plt_idx+leads)],lwd=3)
      if(ld==1){
        mtext('Flow (kcfs)',side=2,line=2,cex=1.5)}
      mtext('Lead (days)',side=1,line=2.5,cex=1.5)
      abline(v=ld,col='darkorange3',lty=2)
    }
  }
}
if(save_plots==T){
  dev.off()}
#################################################################

if(syn_vers==2 & cal_val_setup!='5fold-test'){
  out_DE <- readRDS(paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/out/',loc,'/DEopt_',cal_val_setup,'_pcnt=',opt_pcnt,'_',site,'.rds',sep=''))
  opt_out <- readRDS(paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',opt_pcnt,'_',site,'.rds',sep=''))
  
  ecrps_sse = opt_out[1]
  pwr = opt_out[2]
  hi = opt_out[3]
  lo = opt_out[4]
  sig_a = opt_out[5]
  sig_b = opt_out[6]
  
  decay <- ((1:leads)^pwr/ sum((1:leads)^pwr)) 
  
  dcy <- (decay-min(decay)) * (hi-lo)/(max(decay)-min(decay)) + lo
  
  sigmoid_fun <- function(x,a,b){out <- 1/(1+exp(-x*a+b));return(out)}
  
  if(save_plots==T){
    png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_optim-plot.png',sep=''),height = 1024,width = 3072,res=300)}
  par(mfrow=c(1,3))
  plot(1:leads,dcy,ylim=c(1,1.15*max(dcy)),type='l',lwd=2,col='black',xlab='leads',ylab='ratio threshold',main='Ratio threshold decay')
  text(3,max(dcy) + 0.1*max(dcy),paste('hi:',round(hi,digits=2)),cex=1.5)
  text((leads-2),dcy[leads-2] + 0.1*max(dcy),paste('lo:',round(lo,digits=2)),cex=1.5)
  text(8,dcy[6] + 0.1*max(dcy),paste('pwr:',round(pwr,digits=2)),cex=1.5)
  abline(h=1,lty=2,lwd=2,col='gray')
  
  plot(seq(-5,5,0.01),sigmoid_fun(seq(-5,5,0.01),sig_a,sig_b),xlim=c(-5,5),ylim=c(0,1),type='l',lwd=2,col='black',xlab='scaled obs',ylab='activation value',main='Sigmoid obs dependence')
  text(-3,0.9,paste('slope:',round(sig_a,digits=2)),cex=1.5)
  text(-3,0.8,paste('loc:',round(sig_b,digits=2)),cex=1.5)
  
  conv_x <- 1:out_DE$optim$iter
  conv_y <- out_DE$member$bestvalit
  
  plot(conv_x,conv_y,type='l',ylim=c(0,1.05*max(conv_y)),lwd=2,col='black',xlab='iterations',ylab='eCRPS SSE',main='Optimization convergence')
  points(max(conv_x),conv_y[length(conv_y)],pch=18,cex=3)
  text(0.6*max(conv_x),0.9*max(conv_y),paste('final SSE:',round(ecrps_sse,digits=3)),cex=1.5)
  if(save_plots==T){
    dev.off()}
}

if(syn_vers==2 & cal_val_setup=='5fold-test'){
  for(i in 1:5){
  out_DE <- readRDS(paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/out/',loc,'/DEopt_',cal_val_setup,'_pcnt=',opt_pcnt,'_',site,'-',i,'.rds',sep=''))
  opt_out <- readRDS(paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',opt_pcnt,'_',site,'-',i,'.rds',sep=''))
  
  ecrps_sse = opt_out[1]
  pwr = opt_out[2]
  hi = opt_out[3]
  lo = opt_out[4]
  sig_a = opt_out[5]
  sig_b = opt_out[6]
  
  decay <- ((1:leads)^pwr/ sum((1:leads)^pwr)) 
  
  dcy <- (decay-min(decay)) * (hi-lo)/(max(decay)-min(decay)) + lo
  
  sigmoid_fun <- function(x,a,b){out <- 1/(1+exp(-x*a+b));return(out)}
  
  if(save_plots==T){
    png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_optim-plot-',i,'.png',sep=''),height = 1024,width = 3072,res=300)}
  par(mfrow=c(1,3))
  plot(1:leads,dcy,ylim=c(1,1.15*max(dcy)),type='l',lwd=2,col='black',xlab='leads',ylab='ratio threshold',main='Ratio threshold decay')
  text(3,max(dcy) + 0.1*max(dcy),paste('hi:',round(hi,digits=2)),cex=1.5)
  text((leads-2),dcy[leads-2] + 0.1*max(dcy),paste('lo:',round(lo,digits=2)),cex=1.5)
  text(8,dcy[6] + 0.1*max(dcy),paste('pwr:',round(pwr,digits=2)),cex=1.5)
  abline(h=1,lty=2,lwd=2,col='gray')
  
  plot(seq(-5,5,0.01),sigmoid_fun(seq(-5,5,0.01),sig_a,sig_b),xlim=c(-5,5),ylim=c(0,1),type='l',lwd=2,col='black',xlab='scaled obs',ylab='activation value',main='Sigmoid obs dependence')
  text(-3,0.9,paste('slope:',round(sig_a,digits=2)),cex=1.5)
  text(-3,0.8,paste('loc:',round(sig_b,digits=2)),cex=1.5)
  
  conv_x <- 1:out_DE$optim$iter
  conv_y <- out_DE$member$bestvalit
  
  plot(conv_x,conv_y,type='l',ylim=c(0,1.05*max(conv_y)),lwd=2,col='black',xlab='iterations',ylab='eCRPS SSE',main='Optimization convergence')
  points(max(conv_x),conv_y[length(conv_y)],pch=18,cex=3)
  text(0.6*max(conv_x),0.9*max(conv_y),paste('final SSE:',round(ecrps_sse,digits=3)),cex=1.5)
  if(save_plots==T){
    dev.off()}
}
}

##############Ensemble Density Plots###############

cur_site <- which(site_names==site)

obs <- obs[obs_idx,cur_site]
obs <- obs[obss_idx]
obs_fwd <- obs_forward_all_leads[cur_site,obs_fwd_idx,]
obs_fwd <- obs_fwd[obss_fwd_idx,]
hefs_forward <- hefs_forward[,,hefs_idx,,drop=F]
syn_hefs_forward <- syn_hefs_forward[,,,gen_idx,,drop=F]
syn_hefs_forward <- syn_hefs_forward[,,,genss_idx,,drop=F]

obs_rank <- 1   #pick which obs event to plot (1 largest, 2 second largest, etc)
#get date for plot
#obs_date_loc <- match('1997-01-02',as.character(ixx_obs))
obs_date_loc <- which(obs==sort(obs,decreasing=TRUE)[obs_rank])  #index for maximum observation
obs_date <- ixx_val[obs_date_loc]

if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_event-plot_rnk=',obs_rank,'.png',sep=''),height = 3072,width = 4096,res=300)}
num_plot_rows <- min(n_samp+1,4)
par(mfcol=c(num_plot_rows,4),mar=c(3,3,1,1),mgp=c(3,.4,0))
#ymax <- .5*max(obs[,cur_site],hefs_forward[cur_site,,,],syn_hefs_forward[,cur_site,,,])
ymax <- 2*obs[obs_date_loc]

for (ld in c(1,4,7,10)) {
  plt_idx <- obs_date_loc-ld #ref index for ensemble plots
  
  #1) Forward looking forecast plots
  #hefs
  plot(0:leads,obs[plt_idx:(plt_idx+leads)],type='l',lwd=3,main="",ylim=c(0,ymax),
       xlab='',ylab='')
  abline(v=ld,col='darkorange3',lty=2)
  text(ld+4.5,ymax-1,obs_date,col='darkorange3')
  if (length(obs_date_loc)>0) {
    cur_sum <- round(sum(apply(hefs_forward[cur_site,,plt_idx,],2,FUN=mean)))
    tt <- paste('HEFS ld',ld,'sum',cur_sum)
    title(tt)
    for(i in 1:n_ens){
      lines(0:leads,c(obs[plt_idx],hefs_forward[cur_site,i,plt_idx,]),col=i)
    }
    lines(0:leads,obs[plt_idx:(plt_idx+leads)],lwd=3)
    mtext('Flow (kcfs)',side=2,line=1.8)
    mtext('Days from fore. init.',side=1,line=1.8)
  }
    
  #synthetic
  samples_to_plot <- sample(1:n_samp,(num_plot_rows-1),replace=FALSE)
  for (samp in samples_to_plot) {
    plot(0:leads,obs[plt_idx:(plt_idx+leads)],type='l',lwd=3,main="",ylim=c(0,ymax),
         xlab='',ylab='')
    abline(v=ld,col='darkorange3',lty=2)
    text(ld+4.5,ymax-1,obs_date,col='darkorange3')
    if (length(obs_date_loc)>0) {
      cur_sum <- round(sum(apply(syn_hefs_forward[samp,cur_site,,plt_idx,],2,FUN=mean)))
      tt <- paste('syn',samp,'ld',ld,'sum',cur_sum)
      title(tt)
      for(i in 1:n_ens){
        lines(0:leads,c(obs[plt_idx],syn_hefs_forward[samp,cur_site,i,plt_idx,]),col=i)
      }
      lines(0:leads,obs[plt_idx:(plt_idx+leads)],lwd=3)
      mtext('Flow (kcfs)',side=2,line=1.8)
      mtext('Days from fore. init.',side=1,line=1.8)
    }
  }
}
if(save_plots==T){
  dev.off()}

########################################################



#################Aggregate validation#######################################
lds <- c(1,3,5,7)

if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_auto-correlation.png',sep=''),height = 3072,width = 3072,res=300)}
#check autocorrelation for specific lead times (for one sample, ensemble member and site)
par(mfrow=c(2,2),mar=c(2,2,1,1),oma=c(3,3,1,1))
for (ld in lds) {
  HEFS_acf <- rowMeans(apply(hefs_forward[cur_site,,,ld],1,function(x) {acf(x,plot=FALSE)[[1]]}))[2:11]
  sHEFS_acf <- array(NA,c(length(HEFS_acf),n_samp))
  for (s in 1:n_samp) {
    sHEFS_acf[,s] <- rowMeans(apply(syn_hefs_forward[s,cur_site,,,ld],1,function(x) {acf(x,plot=FALSE)[[1]]}))[2:11]
  }
  sHEFS_acf <- rowMeans(sHEFS_acf)
  plot(HEFS_acf,type="b",ylim=c(0,1))
  lines(sHEFS_acf,type="b",col="red")
  title(paste("lead",ld))
}
if(save_plots==T){
  dev.off()}


#check cross correlation across lead times
if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_lead-cross-correlation.png',sep=''),height = 2048,width = 3072,res=300)}
par(mfcol=c(1,2),mar=c(2,2,1,1),oma=c(3,3,1,1))
HEFS_lead_cor <- lapply(1:n_ens,function(x) {cor(hefs_forward[cur_site,x,,])})
arr <- array( unlist(HEFS_lead_cor) , c(leads,leads,n_ens))
HEFS_lead_cor_ensmean <- apply(arr,1:2,mean) #  Get mean of third dimension
#do the same for the synthetic
sHEFS_lead_cor_ensmean <- array(NA,c(dim(HEFS_lead_cor_ensmean),n_samp))
for (s in 1:n_samp) {
  sHEFS_lead_cor <- lapply(1:n_ens,function(x) {cor(syn_hefs_forward[s,cur_site,x,,])})
  arr <- array( unlist(sHEFS_lead_cor) , c(leads,leads,n_ens))
  sHEFS_lead_cor_ensmean[,,s] <- apply(arr,1:2,mean) #  Get mean of third dimension  
}
sHEFS_lead_cor_ensmean_final <- apply(sHEFS_lead_cor_ensmean,c(1,2),mean)
#plot the correlations across lead times
image.plot(1:leads,1:leads,HEFS_lead_cor_ensmean,zlim=c(0,1))
image.plot(1:leads,1:leads,sHEFS_lead_cor_ensmean_final,zlim=c(0,1))
if(save_plots==T){
  dev.off()}



#check correlation across ensemble members for larger events
num_events <- round((1-disp_pcnt) * length(obs))
obs_date_loc <- order(obs,decreasing=TRUE)[1:num_events]  #index for maximum observation

if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_ens-correlation_quant=',disp_pcnt,'.png',sep=''),height = 3072,width = 3072,res=300)}
par(mfrow=c(4,2),mar=c(2,2,1,1),oma=c(3,3,1,1),mgp=c(3,.4,0))
for (ld in c(1,4,7,10)) {
  hefs_idx <- sapply(obs_date_loc-ld,function(x){max(x,1)})  #the max(,1) just ensures we dont get negative indices
  syn_idx <- sapply(obs_date_loc-ld,function(x) {max(x,1)}) #the max(,1) just ensures we dont get negative indices
  HEFS_ens_cor <- cor(t(hefs_forward[cur_site,,hefs_idx,ld]))
  sHEFS_ens_cor <- array(NA,c(dim(HEFS_ens_cor),n_samp))
  for (s in 1:n_samp) {
    sHEFS_ens_cor[,,s] <- cor(t(syn_hefs_forward[s,cur_site,,syn_idx,ld]))
  }
  sHEFS_ens_cor_final <- apply(sHEFS_ens_cor,c(1,2),mean)
  image.plot(1:n_ens,1:n_ens,HEFS_ens_cor,zlim=c(0,1))
  mtext(paste('HEFS',round(mean(HEFS_ens_cor),2)),side=3,line=0)
  image.plot(1:n_ens,1:n_ens,sHEFS_ens_cor_final,zlim=c(0,1))
  mtext(paste('syn-HEFS',round(mean(sHEFS_ens_cor_final),2)),side=3,line=0)
}
if(save_plots==T){
  dev.off()}
########################################################



###############eCRPS + Rank Histogram####################################
lds<-c(1,3,7,10)  #specify leads (no more than 5 for plotting constraints)
num_events_crps <- round((1-disp_pcnt) * length(obs))
num_events_rankhist <- round((1-disp_pcnt) * length(obs))


###CRPS###
obs_date_loc <- order(obs,decreasing=TRUE)[1:num_events_crps]  #index for maximum observation
obs_major_events <- obs[obs_date_loc]

hefs_ecrps_vec<-array(NA,c(num_events_crps,length(lds)))
shefs_ecrps_vec<-array(NA,c(dim(syn_hefs_forward)[1],num_events_crps,length(lds)))
for(ld in 1:length(lds)){
  for(i in 1:num_events_crps){
    hefs_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    syn_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    HEFS <- hefs_forward[cur_site,,hefs_idx,lds[ld]]
    hefs_ecrps_vec[i,ld] <- eCRPS(HEFS,obs_major_events[i])
    for(s in 1:dim(syn_hefs_forward)[1]){
      SYN_HEFS <- syn_hefs_forward[s,cur_site,,syn_idx,lds[ld]]
      shefs_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,obs_major_events[i])
    }
  }
}

###Cumulative Rank Histogram##
obs_date_loc <- order(obs,decreasing=TRUE)[1:num_events_rankhist]  #index for maximum observation
obs_major_events <- obs[obs_date_loc]

hefs_rank_vec <- array(NA,c(num_events_rankhist,length(lds)))
shefs_rank_vec <- array(NA,c(dim(syn_hefs_forward)[1],num_events_rankhist,length(lds)))
for(ld in 1:length(lds)){
  for(i in 1:num_events_rankhist){
    hefs_idx <- max(obs_date_loc[i]-lds[ld],1) #need to back up by lds[ld] because forecasts are in 'forward' format. the max(,1) ensures we dont get a negative index
    syn_idx <- max(obs_date_loc[i]-lds[ld],1) #need to back up by lds[ld] because forecasts are in 'forward' format. the max(,1) ensures we dont get a negative index
    HEFS <- hefs_forward[cur_site,,hefs_idx,lds[ld]]
    hefs_rank_vec[i,ld] <- ens_rank(HEFS,obs_major_events[i])
    for(s in 1:dim(syn_hefs_forward)[1]){
      SYN_HEFS <- syn_hefs_forward[s,cur_site,,syn_idx,lds[ld]]
      shefs_rank_vec[s,i,ld] <- ens_rank(SYN_HEFS,obs_major_events[i])
    }
  }
}

if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_ecrps-rankhist_quant=',disp_pcnt,'.png',sep=''),height = 2048,width = 3072,res=300)}
par(mfrow=c(2,length(lds)))
for(ld in 1:length(lds)){
  CRPSS <- apply(shefs_ecrps_vec[,,ld],1,function(x){1 - mean(x)/mean(hefs_ecrps_vec[,ld])})
  #boxplot(CRPSS)
  boxplot(hefs_ecrps_vec[,ld],as.vector(shefs_ecrps_vec[,,ld]),
          boxwex=0.5,
          main=paste('ld =',lds[ld],'; top',num_events_crps, 'events'),
          names=c('HEFS','sHEFS'),
          col=c('dodgerblue3','tomato3'),
          ylab='eCRPS',
          range=0,ylim=c(0,median(hefs_ecrps_vec[,ld])+3*IQR(hefs_ecrps_vec[,ld])))
}

my.org=alpha('tomato3',alpha=0.2)
my.gray=alpha('gray',alpha=0.4)
for(i in 1:length(lds)){
  hefs_count <- hist(hefs_rank_vec[,i],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
  hefs_cumul_frac <- roll_sum(hefs_count$counts)/length(hefs_rank_vec[,i])
  plot(0:(n_ens+1),c(0,hefs_cumul_frac),type='l',lwd=3,col='dodgerblue3',
       main=paste('ld =',lds[i],'; top',num_events_rankhist, 'events'),
       xlab='ensemble rank',ylab='cumulative fraction',
       ylim=c(0,1))
  abline(0,1/42,col=my.gray,lwd=2,lty=2)
  legend('topleft',c('HEFS','sHEFS'),lwd=c(3,2),col=c('dodgerblue3','tomato3'))
  for(s in 1:dim(syn_hefs_forward)[1]){
    shefs_count<-hist(shefs_rank_vec[s,,i],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
    shefs_cumul_frac<-roll_sum(shefs_count$counts)/length(hefs_rank_vec[,i])
    lines(0:(n_ens+1),c(0,shefs_cumul_frac),lwd=2,col=my.org)
  }
}
if(save_plots==T){
  dev.off()}
########################################################





##################Event based statistics################################
if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_event-plot_top10.png',sep=''),height = 2048,width = 3072,res=300)}

if (n_samp<2) {
  {stop("cant make these plots without 2 samples or more")}
} else {
  #pick which obs event to plot (1 largest, 2 second largest, etc)
  all_leads <- c(1,4,7,10)
  num_events <- 10
  syn_HEFS_cumul_mean <- array(NA,c(n_samp,num_events,length(all_leads)))
  HEFS_cumul_mean <- array(NA,c(num_events,length(all_leads)))
  obs_cumul_mean <- array(NA,c(num_events,length(all_leads)))
  syn_HEFS_cumul_sd <- array(NA,c(n_samp,num_events,length(all_leads)))
  HEFS_cumul_sd <- array(NA,c(num_events,length(all_leads)))
  
  for (i in 1:num_events) {   
    #get date for plot
    obs_date_loc <- which(obs==sort(obs,decreasing=TRUE)[i])  #index for maximum observation
    obs_date <- ixx_val[obs_date_loc]

  
    for (ld_indx in 1:length(all_leads)) {
      ld <- all_leads[ld_indx]
      obs_idx <- obs_date_loc-ld #ref index for ensemble plots
      obs_forward_idx <- obs_date_loc-ld
      hefs_idx <- obs_date_loc-ld
      syn_idx <- obs_date_loc-ld
    
      #ensemble mean
      HEFS_cumul_mean[i,ld_indx] <- mean(hefs_forward[cur_site,,hefs_idx,])
      syn_HEFS_cumul_mean[,i,ld_indx] <- apply(syn_hefs_forward[,cur_site,,syn_idx,],1,FUN=mean)
      obs_cumul_mean[i,ld_indx] <- mean(obs_events <- obs_fwd[obs_forward_idx,])
      #ensemble standard deviation
      HEFS_cumul_sd[i,ld_indx] <- sd(apply(hefs_forward[cur_site,,hefs_idx,],FUN=mean,1))
      syn_HEFS_cumul_sd[,i,ld_indx] <- apply(apply(syn_hefs_forward[,cur_site,,syn_idx,],c(1,2),FUN=mean),1,FUN=sd)
  
    }
  }
  
  par(mfcol=c(2,4),mar=c(2,2,1,1))
  for (ld_indx in 1:length(all_leads)) {
    ymin <- min(syn_HEFS_cumul_mean[,,ld_indx],HEFS_cumul_mean[,ld_indx],obs_cumul_mean[,ld_indx])
    ymax <- max(syn_HEFS_cumul_mean[,,ld_indx],HEFS_cumul_mean[,ld_indx],obs_cumul_mean[,ld_indx])
    boxplot(syn_HEFS_cumul_mean[,,ld_indx],main=paste("Ens. Mean, Lead",all_leads[ld_indx]),ylim=c(ymin,ymax))
    points(HEFS_cumul_mean[,ld_indx],col="red",pch=16)
    points(obs_cumul_mean[,ld_indx],col="blue",pch=17)
    ymin <- min(syn_HEFS_cumul_sd[,,ld_indx],HEFS_cumul_sd[,ld_indx])
    #ymin <- min(HEFS_cumul_sd[,ld_indx])
    ymax <- max(syn_HEFS_cumul_sd[,,ld_indx],HEFS_cumul_sd[,ld_indx])
    #ymax <- max(HEFS_cumul_sd[,ld_indx])
    boxplot(syn_HEFS_cumul_sd[,,ld_indx],main=paste("Ens. SD, Lead",all_leads[ld_indx]),ylim=c(ymin,ymax))
    points(HEFS_cumul_sd[,ld_indx],col="red",pch=16)
  }
}

if(save_plots==T){
  dev.off()}
########################################################################




##################Skill metrics across events################################

all_leads <- 1:leads
num_events <- round((1-disp_pcnt) * length(obs))

MSE_HEFS <- array(length(all_leads))
MSE_syn_HEFS <- array(length(all_leads))

if (n_samp<2) {
  {stop("cant make these plots without 2 samples or more")}
} else {
  #pick which obs to get skill on (1 largest, 2 second largest, etc)
  syn_HEFS_mean <- array(NA,c(n_samp,num_events,length(all_leads)))
  HEFS_mean <- array(NA,c(num_events,length(all_leads)))
  obs_mean <- array(NA,c(num_events,length(all_leads)))

  #get date for plot
  obs_date_loc <- order(obs,decreasing=TRUE)[1:num_events]  #index for maximum observation
  obs_date <- ixx_val[obs_date_loc]
  
  for (ld_indx in 1:length(all_leads)) {
    ld <- all_leads[ld_indx]
    obs_forward_idx <- obs_date_loc-ld
    hefs_idx <- obs_date_loc-ld
    syn_idx <- obs_date_loc-ld
    
    #ensemble mean
    HEFS_ens_mean <- apply(hefs_forward[cur_site,,hefs_idx,ld],FUN=mean,2)
    syn_HEFS_ens_mean <- apply(syn_hefs_forward[,cur_site,,syn_idx,ld],c(1,3),FUN=mean)
    obs_events <- obs_fwd[obs_forward_idx,ld]

    MSE_HEFS[ld_indx] <- mean(sqrt((HEFS_ens_mean - obs_events)^2))
    MSE_syn_HEFS[ld_indx] <- mean(apply(syn_HEFS_ens_mean,1,function(x) {mean(sqrt((x - obs_events)^2))}))    
  }
}  

if(save_plots==T){
  png(paste('./figs/v',syn_vers,'/',loc,'/',site,'/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_mse_quant=',disp_pcnt,'.png',sep=''),height = 2048,width = 3072,res=300)}
par(mfrow=c(1,1),mar=c(2,2,1,1))
ymax <- 1.2*max(rbind(MSE_HEFS,MSE_syn_HEFS))
barplot(rbind(MSE_HEFS,MSE_syn_HEFS),beside=T,col=c("aquamarine3","coral"),ylim=c(0,ymax))
title("MSE")
legend("topleft", c("HEFS","synHEFS"), pch=15, 
       col=c("aquamarine3","coral"), 
       bty="n")

if(save_plots==T){
  dev.off()}


##################Skill metrics across events################################

all_leads <- 1:leads
num_events <- 30

obs_mean_ratio_HEFS <- array(NA,c(num_events,leads))
obs_sd_ratio_HEFS <- array(NA,c(num_events,leads))
obs_arr <- array(NA,c(num_events,leads))

#get date for plot
obs <- obs[ixx_obs%in%ixx_hefs]
obs_fwd <- obs_fwd[ixx_obs_forward%in%ixx_hefs,]

obs_date_loc <- order(obs,decreasing=TRUE)[1:num_events]  #index for maximum observation
obs_date <- ixx_hefs[obs_date_loc]
  
for (ld_indx in 1:length(all_leads)) {
  ld <- all_leads[ld_indx]
  obs_forward_idx <- obs_date_loc-ld
  hefs_idx <- obs_date_loc-ld
    
  #ensemble mean
  HEFS_ens_mean <- apply(hefs_forward[cur_site,,hefs_idx,ld],FUN=mean,2)
  HEFS_ens_sd <- apply(hefs_forward[cur_site,,hefs_idx,ld],FUN=sd,2)
  obs_events <- obs_fwd[obs_forward_idx,ld]
    
  obs_mean_ratio_HEFS[,ld_indx] <- obs_events/HEFS_ens_mean
  obs_sd_ratio_HEFS[,ld_indx] <- obs_events/HEFS_ens_sd
  obs_arr[,ld_indx] <- obs_events   
}
  

if(save_plots==T){
  png(paste('./figs/comparison/',loc,'/',site,'/',loc,'-',site,'_obs-ensmean_ratios.png',sep=''),height = 2048,width = 3072,res=300)}
par(mfrow=c(3,5),mar=c(2.5,2.5,1,1),mgp=c(1.5,0.5,0))
for(i in 1:leads){
  plot(obs_arr[,i],obs_mean_ratio_HEFS[,i],xlab='obs (kcfs)',ylab='obs/ens-mean ratio',main=paste('ld',i))
  points(obs_arr[1,i],obs_mean_ratio_HEFS[1,i],pch=18,cex=2,col='red')
  #legend('topleft',paste(obs_date[1]),col='red',pch=18,cex=1.5,bty='n')
}
if(save_plots==T){
  dev.off()}

if(save_plots==T){
  png(paste('./figs/comparison/',loc,'/',site,'/',loc,'-',site,'_obs-ensmean_ratios_hist.png',sep=''),height = 2048,width = 3072,res=300)}
par(mfrow=c(3,5),mar=c(2.5,2.5,1,1),mgp=c(1.5,0.5,0))
for(i in 1:leads){
  hist(obs_mean_ratio_HEFS[,i],xlab='obs/ens-mean ratio',ylab='frequency',main=paste('ld',i))
  points(obs_mean_ratio_HEFS[1,i],0,pch=18,cex=2,col='red')
  #legend('topright',paste(obs_date[1]),col='red',pch=18,cex=1.5,bty='n')
}
if(save_plots==T){
  dev.off()}



if(save_plots==T){
  png(paste('./figs/comparison/',loc,'/',site,'/',loc,'-',site,'_obs-enssdev_ratios.png',sep=''),height = 2048,width = 3072,res=300)}
par(mfrow=c(3,5),mar=c(2.5,2.5,1,1),mgp=c(1.5,0.5,0))
ymax <- max(obs_mean_ratio_HEFS)
for(i in 1:leads){
  plot(obs_arr[,i],obs_sd_ratio_HEFS[,i],xlab='obs (kcfs)',ylab='obs/ens-sdev ratio',main=paste('ld',i))
  points(obs_arr[1,i],obs_sd_ratio_HEFS[1,i],pch=18,cex=2,col='red')
  #legend('topleft',paste(obs_date[1]),col='red',pch=18,cex=1.5,bty='n')
}
if(save_plots==T){
  dev.off()}

if(save_plots==T){
  png(paste('./figs/comparison/',loc,'/',site,'/',loc,'-',site,'_obs-enssdev_ratios_hist.png',sep=''),height = 2048,width = 3072,res=300)}
par(mfrow=c(3,5),mar=c(2.5,2.5,1,1),mgp=c(1.5,0.5,0))
for(i in 1:leads){
  hist(obs_sd_ratio_HEFS[,i],xlab='obs/ens-sdev ratio',ylab='frequency',main=paste('ld',i))
  points(obs_sd_ratio_HEFS[1,i],0,pch=18,cex=2,col='red')
  #legend('topright',paste(obs_date[1]),col='red',pch=18,cex=1.5,bty='n')
}
if(save_plots==T){
  dev.off()}

########################################################################

