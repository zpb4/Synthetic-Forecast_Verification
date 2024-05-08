
out = array(NA,c(100,10000))
for(s in 1:100){
  out[s,] = arima.sim(n=10000,list(ar = c(0.5,0.25), ma = c(0)))
}

out_bounds_frac=c()

for(k in 1:100){
  out_bounds<-c()
  for(j in 1:10000){
    if(out[k,j]>max(out[-k,j])|out[k,j]<min(out[-k,j])){out_bounds[j]<-1}
    else{out_bounds[j]<-0}
  }
  out_bounds_frac[k]=sum(out_bounds)/length(out_bounds)
}

#saveRDS(out_bounds_frac,'./src/out_frac.rds')

#out_bounds_frac<-readRDS('z:/Synthetic-Forecast_Verification/src/out_frac.rds')
hist(out_bounds_frac)
mean(out_bounds_frac)
