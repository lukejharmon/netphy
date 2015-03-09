branchLike <-
function(tip.like, bl, q)
    
  {
    
    nb.states<-length(tip.like)
    
    r<-rep(0, nb.states)
    
    p<-MatrixExp.eig(q*bl)
    
    for(i in 1:nb.states)
      r[i]<-logspace_sum(log(p[i,])+tip.like)
    
    return(r)
    
  }
