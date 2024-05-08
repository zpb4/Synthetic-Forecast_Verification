

rm(list=ls())
library(fields)
library(scales)
library(zoo)
#----------------------------------------
syn_vers = 2
loc = 'YRS'
site = 'NBBC1'
parm = 'q'

path = paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))
syn_hefs_forward <- readRDS(paste(path,'out/',loc,'/syn_hefs_forward-',parm,'_plot-ens.rds',sep='')) #use only a slice of larger sample array
ixx_sim <- readRDS(paste(path,'out/',loc,'/ixx_gen.rds',sep='')) 
#n_samp <- readRDS(paste(path,'out/',loc,'/n_samp.rds',sep=''))

n_samp <- 10
source('./src/forecast_verification_functions.R')
#obs<-obs/1000

#################################################################


##############Ensemble Density Plots###############

cur_site <- which(site_names==site)


obs_rank <- 3   #pick which obs event to plot (1 largest, 2 second largest, etc)
#get date for plot
#obs_date_loc <- match('1997-01-02',as.character(ixx_obs))
obs_date_loc <- which(obs[,cur_site]==sort(obs[,cur_site],decreasing=TRUE)[obs_rank])  #index for maximum observation
obs_date <- ixx_obs[obs_date_loc]
hefs_date_loc <- match(obs_date,ixx_hefs)
syn_hefs_date_loc <- match(obs_date,ixx_sim)

num_plot_rows <- min(n_samp+1,4)
par(mfcol=c(num_plot_rows,4),mar=c(3,3,1,1),mgp=c(3,.4,0))
#ymax <- .5*max(obs[,cur_site],hefs_forward[cur_site,,,],syn_hefs_forward[,cur_site,,,])
ymax <- 1.5*max(obs[,cur_site])

for (ld in c(1,4,7,10)) {
  obs_idx <- obs_date_loc-ld #ref index for ensemble plots
  hefs_idx <- hefs_date_loc-ld
  syn_idx <- syn_hefs_date_loc-ld
  
  #1) Forward looking forecast plots
  #hefs
  plot(0:leads,obs[obs_idx:(obs_idx+leads),cur_site],type='l',lwd=3,main="",ylim=c(0,ymax),
       xlab='',ylab='')
  abline(v=ld,col='darkorange3',lty=2)
  text(ld+4.5,ymax-1,obs_date,col='darkorange3')
  if (length(hefs_date_loc)>0) {
    cur_sum <- round(sum(apply(hefs_forward[cur_site,,hefs_idx,],2,FUN=mean)))
    tt <- paste('HEFS ld',ld,'sum',cur_sum)
    title(tt)
    for(i in 1:n_ens){
      lines(0:leads,c(obs[obs_idx,cur_site],hefs_forward[cur_site,i,hefs_idx,]),col=i)
    }
    lines(0:leads,obs[obs_idx:(obs_idx+leads),cur_site],lwd=3)
    mtext('Flow (kcfs)',side=2,line=1.8)
    mtext('Days from fore. init.',side=1,line=1.8)
  }
    
  #synthetic
  samples_to_plot <- sample(1:n_samp,(num_plot_rows-1),replace=FALSE)
  for (samp in samples_to_plot) {
    plot(0:leads,obs[obs_idx:(obs_idx+leads),cur_site],type='l',lwd=3,main="",ylim=c(0,ymax),
         xlab='',ylab='')
    abline(v=ld,col='darkorange3',lty=2)
    text(ld+4.5,ymax-1,obs_date,col='darkorange3')
    if (length(syn_hefs_date_loc)>0) {
      cur_sum <- round(sum(apply(syn_hefs_forward[samp,cur_site,,syn_idx,],2,FUN=mean)))
      tt <- paste('syn',samp,'ld',ld,'sum',cur_sum)
      title(tt)
      for(i in 1:n_ens){
        lines(0:leads,c(obs[obs_idx,cur_site],syn_hefs_forward[samp,cur_site,i,syn_idx,]),col=i)
      }
      lines(0:leads,obs[obs_idx:(obs_idx+leads),cur_site],lwd=3)
      mtext('Flow (kcfs)',side=2,line=1.8)
      mtext('Days from fore. init.',side=1,line=1.8)
    }
  }
}


########################################################



#################Aggregate validation#######################################
lds <- c(1,3,5,7)

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


#check cross correlation across lead times
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




#check correlation across ensemble members for larger events
num_events <- 50
keep <- ixx_obs%in%ixx_hefs
obs_date_loc <- order(obs[keep,cur_site],decreasing=TRUE)[1:num_events]  #index for maximum observation
obs_date <- ixx_obs[keep][obs_date_loc]
hefs_date_loc <- match(obs_date,ixx_hefs)
syn_hefs_date_loc <- match(obs_date,ixx_sim)

par(mfrow=c(4,2),mar=c(2,2,1,1),oma=c(3,3,1,1),mgp=c(3,.4,0))
for (ld in c(1,4,7,10)) {
  hefs_idx <- sapply(hefs_date_loc-ld,function(x){max(x,1)})  #the max(,1) just ensures we dont get negative indices
  syn_idx <- sapply(syn_hefs_date_loc-ld,function(x) {max(x,1)}) #the max(,1) just ensures we dont get negative indices
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
########################################################



###############eCRPS + Rank Histogram####################################
lds<-c(1,3,7,10)  #specify leads (no more than 5 for plotting constraints)
num_events_crps <- 30
num_events_rankhist <- 100


###CRPS###
keep <- ixx_obs%in%ixx_hefs
obs_date_loc <- order(obs[keep,cur_site],decreasing=TRUE)[1:num_events_crps]  #index for maximum observation
obs_date <- ixx_obs[keep][obs_date_loc]
hefs_date_loc <- match(obs_date,ixx_hefs)
syn_hefs_date_loc <- match(obs_date,ixx_sim)
obs_major_events <- obs[keep,cur_site][obs_date_loc]

hefs_ecrps_vec<-array(NA,c(num_events_crps,length(lds)))
shefs_ecrps_vec<-array(NA,c(dim(syn_hefs_forward)[1],num_events_crps,length(lds)))
for(ld in 1:length(lds)){
  for(i in 1:num_events_crps){
    hefs_idx <- hefs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    syn_idx <- syn_hefs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    HEFS <- hefs_forward[cur_site,,hefs_idx,lds[ld]]
    hefs_ecrps_vec[i,ld] <- eCRPS(HEFS,obs_major_events[i])
    for(s in 1:dim(syn_hefs_forward)[1]){
      SYN_HEFS <- syn_hefs_forward[s,cur_site,,syn_idx,lds[ld]]
      shefs_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,obs_major_events[i])
    }
  }
}

###Cumulative Rank Histogram##
keep <- ixx_obs%in%ixx_hefs
obs_date_loc <- order(obs[keep,cur_site],decreasing=TRUE)[1:num_events_rankhist]  #index for maximum observation
obs_date <- ixx_obs[keep][obs_date_loc]
hefs_date_loc <- match(obs_date,ixx_hefs)
syn_hefs_date_loc <- match(obs_date,ixx_sim)
obs_major_events <- obs[keep,cur_site][obs_date_loc]

hefs_rank_vec <- array(NA,c(num_events_rankhist,length(lds)))
shefs_rank_vec <- array(NA,c(dim(syn_hefs_forward)[1],num_events_rankhist,length(lds)))
for(ld in 1:length(lds)){
  for(i in 1:num_events_rankhist){
    hefs_idx <- max(hefs_date_loc[i]-lds[ld],1) #need to back up by lds[ld] because forecasts are in 'forward' format. the max(,1) ensures we dont get a negative index
    syn_idx <- max(syn_hefs_date_loc[i]-lds[ld],1) #need to back up by lds[ld] because forecasts are in 'forward' format. the max(,1) ensures we dont get a negative index
    HEFS <- hefs_forward[cur_site,,hefs_idx,lds[ld]]
    hefs_rank_vec[i,ld] <- ens_rank(HEFS,obs_major_events[i])
    for(s in 1:dim(syn_hefs_forward)[1]){
      SYN_HEFS <- syn_hefs_forward[s,cur_site,,syn_idx,lds[ld]]
      shefs_rank_vec[s,i,ld] <- ens_rank(SYN_HEFS,obs_major_events[i])
    }
  }
}



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
          range=0,ylim=c(0,1.5*max(hefs_ecrps_vec[,ld])))
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
########################################################





##################Event based statistics################################

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
    keep <- ixx_obs%in%ixx_hefs
    obs_date_loc <- which(obs[keep,cur_site]==sort(obs[keep,cur_site],decreasing=TRUE)[i])  #index for maximum observation
    obs_date <- ixx_obs[keep][obs_date_loc]
    obs_forward_date_loc <- which(ixx_obs_forward==obs_date)
    hefs_date_loc <- which(ixx_hefs==obs_date)
    syn_hefs_date_loc <- which(ixx_sim==obs_date)
  
    for (ld_indx in 1:length(all_leads)) {
      ld <- all_leads[ld_indx]
      obs_idx <- obs_date_loc-ld #ref index for ensemble plots
      obs_forward_idx <- obs_forward_date_loc-ld
      hefs_idx <- hefs_date_loc-ld
      syn_idx <- syn_hefs_date_loc-ld
    
      #ensemble mean
      HEFS_cumul_mean[i,ld_indx] <- mean(hefs_forward[cur_site,,hefs_idx,])
      syn_HEFS_cumul_mean[,i,ld_indx] <- apply(syn_hefs_forward[,cur_site,,syn_idx,],1,FUN=mean)
      obs_cumul_mean[i,ld_indx] <- mean(obs_forward_all_leads[cur_site,obs_forward_idx,])
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
########################################################################




##################Skill metrics across events################################

all_leads <- 1:leads
num_events <- 100

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
  keep <- ixx_obs%in%ixx_hefs
  obs_date_loc <- order(obs[keep,cur_site],decreasing=TRUE)[1:num_events]  #index for maximum observation
  obs_date <- ixx_obs[keep][obs_date_loc]
  obs_forward_date_loc <- which(ixx_obs_forward%in%obs_date)
  hefs_date_loc <- which(ixx_hefs%in%obs_date)
  syn_hefs_date_loc <- which(ixx_sim%in%obs_date)
  
  for (ld_indx in 1:length(all_leads)) {
    ld <- all_leads[ld_indx]
    obs_forward_idx <- obs_forward_date_loc-ld
    hefs_idx <- hefs_date_loc-ld
    syn_idx <- syn_hefs_date_loc-ld
    
    #ensemble mean
    HEFS_ens_mean <- apply(hefs_forward[cur_site,,hefs_idx,ld],FUN=mean,2)
    syn_HEFS_ens_mean <- apply(syn_hefs_forward[,cur_site,,syn_idx,ld],c(1,3),FUN=mean)
    obs_events <- obs_forward_all_leads[cur_site,obs_forward_idx,ld]

    MSE_HEFS[ld_indx] <- mean(sqrt((HEFS_ens_mean - obs_events)^2))
    MSE_syn_HEFS[ld_indx] <- mean(apply(syn_HEFS_ens_mean,1,function(x) {mean(sqrt((x - obs_events)^2))}))    
  }
}  
  

par(mfrow=c(1,1),mar=c(2,2,1,1))
ymax <- 1.2*max(rbind(MSE_HEFS,MSE_syn_HEFS))
barplot(rbind(MSE_HEFS,MSE_syn_HEFS),beside=T,col=c("aquamarine3","coral"),ylim=c(0,ymax))
title("MSE")
legend("topleft", c("HEFS","synHEFS"), pch=15, 
       col=c("aquamarine3","coral"), 
       bty="n")

########################################################################

