
load("z:/Synthetic-Forecast-v1-FIRO-DISES/out/model-fit_rdata.RData")

fwd_forecast_rearrange<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,(i+1):dim(forecast)[2],i]<-forecast[,1:(dim(forecast)[2]-i),i]
  }
  return(forecast_out)
}


hefs_mat = hefs_fit[2,,,]
obs=obs_fit$NHGC1
hefs_mat[hefs_mat>0]<-1

frac_fun<-function(x){out = sum(x)/length(x);return(out)}

hefs_nzero_frac = apply(hefs_mat,c(2,3),frac_fun)

q = quantile(obs,probs = seq(0, 1, 0.1))
qvec = unique(q)
qvec[length(qvec)]=Inf

#look at returns for different specific values
cutt = cut(0,breaks=qvec,include.lowest = T)


ct = cut(obs,breaks=qvec,include.lowest = T)
ct_idx = cut(obs,breaks=qvec,include.lowest = T,labels=1:(length(qvec)-1))
ct_idx = as.numeric(ct_idx)

ld = 1

plt_arr <- array(NA,c(2,(length(qvec)-1)))
for(i in 1:(length(qvec)-1)){
  ix = which(ct_idx==i)
  plt_arr[1,i]=median(obs[ix])
  plt_arr[2,i]=mean(hefs_nzero_frac[ix,ld])
}

plot(obs,hefs_nzero_frac[,ld])
lines(plt_arr[1,],plt_arr[2,],col='green')
