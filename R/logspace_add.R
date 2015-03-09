logspace_add <-
function(logx, logy) {
  if(logx==-Inf) return(logy) else max(logx, logy) + log1p(exp (-abs (logx - logy)));
}
