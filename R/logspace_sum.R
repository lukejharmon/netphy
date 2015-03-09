logspace_sum <-
function(logx) {
  r<-logx[1]
  if(length(logx)>1)
    for(i in 2:length(logx))
      r<-logspace_add(r, logx[i])
  r  
}
