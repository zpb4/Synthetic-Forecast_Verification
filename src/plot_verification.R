#script to plot forecast verification metrics for synthetic forecasts against 
#observed data and actual HEFS forecasts
#----------------------------------------
rm(list=ls());gc()
setwd('z:/Synthetic-Forecast_Verification/')

library(fGarch)
library(stringr)
library(scales)
library(dplyr)
source('./src/forecast_verification_functions.R')

syn_vers<-'v2'
syn_forecast_path <- paste('z:/Synthetic-Forecast-',syn_vers,'-FIRO-DISES/',sep='')
v1_mod<-'bvar1-sged' #specify model type for the 'v1' approach

load(paste(syn_forecast_path,'out/model-fit_rdata.RData',sep=''))

if(syn_vers=='v1'){
  syn_hefs_forward <- readRDS(paste(syn_forecast_path,'out/syn_hefs_forward_',v1_mod,'.rds',sep=''))
  syn_err<-readRDS(paste(syn_forecast_path,'out/syn_err_',v1_mod,'.rds',sep=''))
  sep_est<-readRDS(paste(syn_forecast_path,'out/sep-fit_',v1_mod,'.rds',sep=''))
  var_est<-readRDS(paste(syn_forecast_path,'out/var-fit_',v1_mod,'.rds',sep=''))
  source(paste(syn_forecast_path,'src/functions/sep_function.R',sep=''))
  dist<-str_split(v1_mod,'-')[[1]][2]}

if(syn_vers=='v2'){
  syn_hefs_forward <- readRDS(paste(syn_forecast_path,'out/syn_hefs_forward.rds',sep=''))
}

site<-'NHGC1'
cur_site <- which(site_names==site)
mth_sset<-c(12,1,2,3) #specify monthly subset to look at in a vector
lds<-c(1,4,7,10)  #specify leads (no more than 5 for plotting constraints)
ens_samp<-1
bins<-20 # no of bins for binned spread-mse diagrams
pcnt<-c(0.9,1) # quantile (lwr,upr) for cumul rank histograms and eCRPS
data_sset<-'cal' # 'cal', 'val', or 'test'
plt_out<-T # T to show plots, F to save png files

ixx_gen <- readRDS(paste(syn_forecast_path,'out/ixx_gen.rds',sep='')) 
#ixx_sim_shefs <- as.POSIXlt(seq(as.Date('1979-10-02'),as.Date('2021-09-16'),by='day'),tz = "UTC")
n_samp <- readRDS(paste(syn_forecast_path,'out/n_samp.rds',sep='')) 
ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'),tz = "UTC")
cal_idx <- !(ixx_fit_all$year%in%(leave_out_years-1900))
ixx_cal <- ixx_fit_all[cal_idx] #fit years excluding leave out years
val_idx <- ixx_fit_all$year%in%(leave_out_years-1900)
ixx_val <- ixx_fit_all[val_idx] #fit years excluding leave out years

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
if(length(mth_sset)==1){mth_lab<-mths[mth_sset]}
if(length(mth_sset)>1){mth_lab<-paste(mths[mth_sset[1]],'-',mths[mth_sset[length(mth_sset)]])}

if(data_sset=='cal'){
  sset_hefs<-ixx_hefs%in%ixx_cal
  sset_shefs<-ixx_gen%in%ixx_cal
  sset_obs<-ixx_obs%in%ixx_cal
  ixx_plt_hefs<-ixx_hefs[sset_hefs]
  ixx_plt_shefs<-ixx_gen[sset_shefs]
  ixx_plt_obs<-ixx_obs[sset_obs]
  mth_sset_idx_hefs<-which(ixx_plt_hefs$mo%in%(mth_sset-1))
  mth_sset_idx_shefs<-which(ixx_plt_shefs$mo%in%(mth_sset-1))
  mth_sset_idx_obs<-which(ixx_plt_obs$mo%in%(mth_sset-1))
  obs_sset<-obs[sset_obs,cur_site][mth_sset_idx_obs]
  obs_sset_hefs<-obs_sset}

if(data_sset=='val'){
  sset_hefs<-ixx_hefs%in%ixx_val
  sset_shefs<-ixx_gen%in%ixx_val
  sset_obs<-ixx_obs%in%ixx_val
  ixx_plt_hefs<-ixx_hefs[sset_hefs]
  ixx_plt_shefs<-ixx_gen[sset_shefs]
  ixx_plt_obs<-ixx_obs[sset_obs]
  mth_sset_idx_hefs<-which(ixx_plt_hefs$mo%in%(mth_sset-1))
  mth_sset_idx_shefs<-which(ixx_plt_shefs$mo%in%(mth_sset-1))
  mth_sset_idx_obs<-which(ixx_plt_obs$mo%in%(mth_sset-1))
  obs_sset<-obs[sset_obs,cur_site][mth_sset_idx_obs]
  obs_sset_hefs<-obs_sset}

if(data_sset=='test'){
  sset_hefs<-ixx_hefs%in%ixx_val #use HEFS val data for comparison in test period
  sset_shefs<-!(ixx_gen%in%ixx_fit_all)
  sset_obs<-ixx_obs%in%ixx_gen
  sset_obs_hefs<-ixx_obs%in%ixx_val
  ixx_plt_hefs<-ixx_hefs[sset_hefs]
  ixx_plt_shefs<-ixx_gen[sset_shefs]
  #ixx_plt_obs<-ixx_obs[sset_obs]
  ixx_plt_obs_hefs<-ixx_obs[sset_obs_hefs]
  mth_sset_idx_hefs<-which(ixx_plt_hefs$mo%in%(mth_sset-1))
  mth_sset_idx_shefs<-which(ixx_plt_shefs$mo%in%(mth_sset-1))
  #mth_sset_idx_obs<-which(ixx_plt_obs$mo%in%(mth_sset-1))
  mth_sset_idx_obs_hefs<-which(ixx_plt_obs_hefs$mo%in%(mth_sset-1))
  obs_gen<-obs[sset_obs,cur_site][sset_shefs]
  obs_sset<-obs_gen[mth_sset_idx_shefs]
  obs_sset_hefs<-obs[sset_obs_hefs,cur_site][mth_sset_idx_obs_hefs]}

cur_dir <- paste("./figs/",syn_vers,sep='')
if (!dir.exists(cur_dir)) {
  dir.create(cur_dir)
}
#-------------------------------------------------------------------
#plot a_t fit comparison
#only plot if 'v1' approach

if(syn_vers=='v1' & (data_sset=='cal'|data_sset=='val')){
  par(mfrow=c(length(mth_sset),length(lds)))
  xlm=c(-3,3)

  for(j in 1:length(mth_sset)){
    mth_sset_ix<-which(ixx_plt_hefs$mo%in%(mth_sset[j]-1))
    for(i in 1:length(lds)){
      at_sset<-var_est[[1]][ens_samp,sset_hefs,lds[i]]
      at<-at_sset[mth_sset_ix]
      par<-sep_est[ens_samp,mth_sset[j],lds[i],]
      brks<-c(min(at),seq(-5,5,0.2),max(at))
      hist(at,main=paste('mth =',mths[mth_sset[j]],'; ld =',lds[i]),xlab='a_t',sub=data_sset,breaks = brks,xlim=xlm,freq = FALSE,col='gray90')
      if(dist=='sged'){
        lines(seq(-5,5,0.1),dsged(seq(-5,5,0.1),mean=par[1],sd=par[2],nu=par[3],xi=par[4]),col='tomato3',lwd=2)
      }
      if(dist=='sep'){
        lines(seq(-5,5,0.1),sep_pdf(seq(-5,5,0.1),theta=par[1],sigma=par[2],beta=par[3],alpha=par[4]),col='tomato3',lwd=2)
      }
    }
  }
}

#----------------------------------------------------------------------
#plot raw error autocorrelations
if(syn_vers=='v1'){
  par(mfrow=c(1,length(lds)))
  syn_err_fit_idx<-syn_err[,,sset_shefs,]
  hefs_err_fit_idx<-err_mat[,sset_hefs,]
  lagmax<-10
  my.gray=alpha('gray',alpha=0.4)

  for(i in 1:length(lds)){
    err_emp<-err_mat[ens_samp,mth_sset_idx_hefs,lds[i]]
    err_syn<-syn_err_fit_idx[,ens_samp,mth_sset_idx_shefs,lds[i]]
    err_syn_acf<-apply(err_syn,1,function(x){out=acf(x,lag.max = lagmax,plot=F);return(out$acf)})
    max_acf<-apply(err_syn_acf,1,max)
    min_acf<-apply(err_syn_acf,1,min)
    acf(err_emp,lag.max = lagmax,main=paste('mth =',mth_lab,'; ld =',lds[i]),sub=data_sset)
    polygon(x=c(0:lagmax,lagmax:0),y=c(min_acf,rev(max_acf)),col=my.gray,border=NA)
  }
}

#---------------------------------------------------------------------
#plot cumulative rank histogram
my.org=alpha('tomato3',alpha=0.2)
my.gray=alpha('gray',alpha=0.4)

syn_forc_obs_synch<-array(NA,dim(syn_hefs_forward[,cur_site,,sset_shefs,]))
for(i in 1:dim(syn_forc_obs_synch)[1]){
  syn_forc_obs_synch[i,,,]<-fwd_forecast_rearrange(syn_hefs_forward[i,cur_site,,sset_shefs,])
}
syn_forc_obs_synch<-syn_forc_obs_synch[,,mth_sset_idx_shefs,]

hefs_obs_synch<-fwd_forecast_rearrange(hefs_forward[cur_site,,,])
hefs <- hefs_obs_synch[,sset_hefs,]
hefs <-hefs[,mth_sset_idx_hefs,]

sort_obs_idx<-sort(obs_sset,index.return=TRUE)
upr_idx<-round(pcnt[2]*length(sort_obs_idx$ix))
lwr_idx<-max(1,round(pcnt[1]*length(sort_obs_idx$ix)))
samp_idx<-sort_obs_idx$ix[lwr_idx:upr_idx]

sort_obs_idx_hefs<-sort(obs_sset_hefs,index.return=TRUE)
upr_idx_hefs<-round(pcnt[2]*length(sort_obs_idx_hefs$ix))
lwr_idx_hefs<-max(1,round(pcnt[1]*length(sort_obs_idx_hefs$ix)))
samp_idx_hefs<-sort_obs_idx_hefs$ix[lwr_idx_hefs:upr_idx_hefs]

hefs_rank_vec<-array(NA,c(length(samp_idx_hefs),length(lds)))
shefs_rank_vec<-array(NA,c(dim(syn_forc_obs_synch)[1],length(samp_idx),length(lds)))

for(k in 1:length(lds)){
  for(i in 1:length(samp_idx_hefs)){
    hefs_rank_vec[i,k]<-ens_rank(hefs[,samp_idx_hefs[i],lds[k]],obs_sset_hefs[samp_idx_hefs[i]])
  }
}

for(k in 1:length(lds)){
  for(i in 1:length(samp_idx)){
    for(s in 1:dim(syn_forc_obs_synch)[1]){
      shefs_rank_vec[s,i,k]<-ens_rank(syn_forc_obs_synch[s,,samp_idx[i],lds[k]],obs_sset[samp_idx[i]])
    }
  }
}


#plot
par(mfrow=c(1,length(lds)))

for(i in 1:length(lds)){
  hefs_count<-hist(hefs_rank_vec[,i],breaks=seq(0.5,dim(hefs)[1]+1.5),plot = FALSE)
  hefs_cumul_frac<-roll_sum(hefs_count$counts)/length(hefs_rank_vec[,i])
  plot(0:dim(hefs)[1],hefs_cumul_frac,type='l',lwd=3,col='dodgerblue3',main=c(paste('mth =',mth_lab,'; ld =',lds[i]),paste('quantile =',pcnt[1],'-',pcnt[2])),
       xlab='ensemble rank',ylab='cumulative fraction',sub=data_sset)
  abline(0,1/42,col=my.gray,lwd=2,lty=2)
  legend('topleft',c('HEFS','sHEFS'),lwd=c(3,2),col=c('dodgerblue3','tomato3'))
  for(s in 1:dim(syn_forc_obs_synch)[1]){
    shefs_count<-hist(shefs_rank_vec[s,,i],breaks=seq(0.5,dim(hefs)[1]+1.5),plot = FALSE)
    shefs_cumul_frac<-roll_sum(shefs_count$counts)/length(hefs_rank_vec[,i])
    lines(0:dim(hefs)[1],shefs_cumul_frac,lwd=2,col=my.org)
  }
  
}

#------------------------------------------------------------------------------
#plot eCRPS

hefs_ecrps_vec<-array(NA,c(length(samp_idx_hefs),length(lds)))
shefs_ecrps_vec<-array(NA,c(dim(syn_forc_obs_synch)[1],length(samp_idx),length(lds)))

for(k in 1:length(lds)){
  for(i in 1:length(samp_idx_hefs)){
    hefs_ecrps_vec[i,k]<-eCRPS(hefs[,samp_idx_hefs[i],lds[k]],obs_sset_hefs[samp_idx_hefs[i]])
  }
}

for(k in 1:length(lds)){
  for(i in 1:length(samp_idx)){
    for(s in 1:dim(syn_forc_obs_synch)[1]){
      shefs_ecrps_vec[s,i,k]<-eCRPS(syn_forc_obs_synch[s,,samp_idx[i],lds[k]],obs_sset[samp_idx[i]])
    }
  }
}

par(mfrow=c(1,length(lds)))
for(i in 1:length(lds)){
  boxplot(hefs_ecrps_vec[,i],as.vector(shefs_ecrps_vec[,,i]),
          boxwex=0.5,
          main=c(paste('mth =',mth_lab,'; ld =',lds[i]),paste('quantile =',pcnt[1],'-',pcnt[2])),
          names=c('HEFS','sHEFS'),
          sub=data_sset,
          col=c('dodgerblue3','tomato3'),
          ylab='eCRPS',
          ylim=c(0,median(hefs_ecrps_vec[,i])+2*IQR(hefs_ecrps_vec[,i])),
          range=0)
}


#--------------------------------------------------------------------------
#bin-spread diagrams
ens_mns_hefs<-apply(hefs[,,],c(2,3),mn_ens)
ens_vars_hefs<-apply(hefs[,,],c(2,3),st2_ens)
bin_mse_hefs<-array(NA,c(length(lds),bins))
bin_spd_hefs<-array(NA,c(length(lds),bins))

ens_mns_shefs<-apply(syn_forc_obs_synch[,,,],c(1,3,4),mn_ens)
ens_vars_shefs<-apply(syn_forc_obs_synch[,,,],c(1,3,4),st2_ens)
bin_mse_shefs<-array(NA,c(dim(syn_forc_obs_synch)[1],length(lds),bins))
bin_spd_shefs<-array(NA,c(dim(syn_forc_obs_synch)[1],length(lds),bins))

for(s in 1:dim(syn_forc_obs_synch)[1]){
  for(k in 1:length(lds)){
    ct_hefs<-cut(rank(ens_vars_hefs[,k],ties.method = 'first'),20,labels=1:20)
    ct_shefs<-cut(rank(ens_vars_shefs[s,,k],ties.method = 'first'),20,labels=1:20)
    for(j in 1:bins){
      hefs_var_sset<-ens_vars_hefs[which(ct_hefs==j),k]
      hefs_mn_sset<-ens_mns_hefs[which(ct_hefs==j),k]
      bin_spd_hefs[k,j]<-bin_spread(hefs_var_sset,dim(hefs)[1])
      bin_mse_hefs[k,j]<-bin_mse(hefs_mn_sset,obs_sset_hefs[which(ct_hefs==j)],dim(hefs)[1])
      
      shefs_var_sset<-ens_vars_shefs[s,which(ct_shefs==j),k]
      shefs_mn_sset<-ens_mns_shefs[s,which(ct_shefs==j),k]
      bin_spd_shefs[s,k,j]<-bin_spread(shefs_var_sset,dim(hefs)[1])
      bin_mse_shefs[s,k,j]<-bin_mse(shefs_mn_sset,obs_sset[which(ct_shefs==j)],dim(hefs)[1])
    }
  }
}

par(mfrow=c(1,length(lds)))
for(i in 1:length(lds)){
  xlm<-c(0,1.25*max(bin_spd_hefs[i,],bin_mse_hefs[i,]))
  ylm<-c(0,1.25*max(bin_spd_hefs[i,],bin_mse_hefs[i,]))
  plot(bin_spd_hefs[i,],bin_mse_hefs[i,],type='l',lwd=3,col='dodgerblue3',main=paste('mth =',mth_lab,'; ld =',lds[i]),
    xlab='ensemble spread',ylab='ensemble mse',xlim=xlm,ylim=ylm,sub=data_sset)
  abline(0,1,lwd=2,lty=2,col=my.gray)
  legend('topleft',c('HEFS','sHEFS'),lwd=c(3,2),col=c('dodgerblue3','tomato3'))
  for(s in 1:dim(syn_forc_obs_synch)[1]){
    lines(bin_spd_shefs[s,i,],bin_mse_shefs[s,i,],lwd=2,col=my.org)
  }
}


####################################END####################################