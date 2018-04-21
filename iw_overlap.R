iw_overlap<-function(ps,E,a,b){
  ps1<-ps[which(E==1)]
  ps0<-ps[which(E==0)]
  ps1<-ps1[order(ps1)]
  ps0<-ps0[order(ps0)]
  
  RO<-rep(0,length(ps))
  for (k in ps){
    cte<-0
    for (e in 0:1){
      temp<-get(paste('ps',e,sep=''))
      fooless<-temp[which(temp<k)]
      foomore<-temp[which(temp>k)]
      if (length(fooless)>=b & length(foomore)>=b){
        
        allcheck<-c(fooless[(length(fooless)-(b-1)):length(fooless)],k,foomore[1:b])
        
      } else if (length(fooless)>=b & length(foomore)<b){
        
        if (k %in% temp){
          allcheck<-c(fooless[(length(fooless)-(b-1)):length(fooless)],k,foomore,rep(Inf,b-length(foomore)))
        } else{
          allcheck<-c(fooless[(length(fooless)-(b-1)):length(fooless)],k,foomore,max(ps),rep(Inf,b-(length(foomore)+1)))
        }
        
      } else if (length(fooless)<b & length(foomore)>=b){
        
        if (k %in% temp){
          allcheck<-c(rep(-Inf,b-length(fooless)),fooless,k,foomore[1:b])
        } else{
          allcheck<-c(rep(-Inf,b-(length(fooless)+1)),min(ps),fooless,k,foomore[1:b])
        }
        
      } else{
        #allcheck<-c(rep(-Inf,b-length(fooless)),fooless,k,foomore,rep(Inf,b-length(foomore)))
        stop(paste('Less than b propensity scores above and below PS',k,'in exposure group',e,', choose a smaller b value.',sep=' '))
      }

      cth<-0
      if (k %in% temp){
        for (h in 1:(b+1)){
          if (abs(allcheck[h]-allcheck[h+b])<a) cth<-cth+1
        }
      } else{
        for (h in 1:b){
          if (abs(allcheck[h]-allcheck[h+b+1])<a) cth<-cth+1
        }
      }
      
      if (cth>0) cte<-cte+1
    }
    if (cte==2) RO[which(ps==k)]<-1
  }
  
  return(RO)
}