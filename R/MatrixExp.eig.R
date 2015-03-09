MatrixExp.eig <-
function(Q)
  {
    tmp <- eigen(Q, symmetric = FALSE)
    P1 <- tmp$vectors %*% diag(exp(tmp$values)) %*% solve(tmp$vectors)
    return(P1)
    
  }
