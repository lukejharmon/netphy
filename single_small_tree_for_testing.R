set.seed(42)

small.short.tree <- TreeSim::sim.bd.taxa.age(n = 4, numbsim = 1, lambda = 0.4, mu = 0, age = 2, mrca = TRUE)
small.long.tree <- TreeSim::sim.bd.taxa.age(n = 4, numbsim = 1, lambda = 0.4, mu = 0, age = 10, mrca = TRUE)

small.short.net <- simPhyloNetwork(small.short.tree[[1]], qRate = 0.1, sProb = 0.5)
small.long.net <- simPhyloNetwork(small.long.tree[[1]], qRate = 0.1, sProb = 0.5)


foo.short <- function(par) {
    q01.small <- exp(par[1])
    q10.small <- exp(par[2])
    pSpec <- exp(par[3])
    fitPhyloNetwork(small.short.tree[[1]], small.short.net, c(q01.small, q01.small, pSpec))
}

foo.long <- function(par) {
    q01.small <- exp(par[1])
    q10.small <- exp(par[2])
    pSpec <- exp(par[3])
    fitPhyloNetwork(small.long.tree[[1]], small.long.net, c(q01.small, q01.small, pSpec))
}



fit.small.short <- optim(log(c(0.1, 0.1, 0.5)), fn = foo.short, control = list(trace = 6, reltol = .Machine$double.eps^.25))
fit.small.long <- optim(log(c(0.1, 0.1, 0.5)), fn = foo.long, control = list(trace = 6, reltol = .Machine$double.eps^.25))

exp(fit.small.short$par)
exp(fit.small.long$par)




### From examples.R

tt<-birthdeath.tree(b=0.1, d=0, taxa.stop=5)
tt<-drop.tip(tt, "s5")
xx<-simPhyloNetwork(tt, qRate=0.1, sProb=0.5)
fitPhyloNetwork(tt, xx, c(q01=0.1, q10=0.1, pSpec=0.5))

foo<-function(x) {
  q01<-exp(x[1])
  q10<-exp(x[2])
  pSpec<-exp(x[3])
  fitPhyloNetwork(tt, xx, c(q01, q10, pSpec))
}

foo(log(c(0.1, 0.1, 0.5)))
out<-optim(log(c(0.5, 0.5, 0.5)), foo, control=list(trace=6))
exp(out$par)

