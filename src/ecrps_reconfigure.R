
ens = 1:42


ecrps_for = function(ens){
  out = matrix(0,ncol=length(ens),nrow=length(ens))
  for(i in 1:(length(ens)-1)){
    for(k in (i+1):length(ens)){
      out[i,k] = abs(ens[i]-ens[k])
    }
  }
  return(out)
}

ec = sum(ecrps_for(ens))


ecrps_mat = function(ens){
  mat1 = matrix(rep(ens,each=length(ens)),ncol=length(ens),byrow=F)
  mat2 = matrix(rep(ens,length(ens)),ncol=length(ens),byrow=F)
  diff = abs(mat1[lower.tri(mat1)]-mat2[lower.tri(mat2)])
  return(diff)
}

ec2 = sum(ecrps_mat(ens))