

## Code lifted from ape function "ace"
MatrixExp.eig<-
  function(Q)
  {
    tmp <- eigen(Q, symmetric = FALSE)
    P1 <- tmp$vectors %*% diag(exp(tmp$values)) %*% solve(tmp$vectors)
    return(P1)
    
  }

logspace_add<-function(logx, logy) {
  if(logx==-Inf) return(logy) else max(logx, logy) + log1p(exp (-abs (logx - logy)));
}

logspace_sum<-function(logx) {
  r<-logx[1]
  if(length(logx)>1)
    for(i in 2:length(logx))
      r<-logspace_add(r, logx[i])
  r  
}

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