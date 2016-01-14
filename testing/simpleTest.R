library(devtools)
install_github("lukejharmon/netphy")


library(geiger)
library(np)
phy<-read.tree(text="(((a:1.0,b:1.0):2.0,(c:0.5,d:0.5):2.5):0.5,e:3.5);")
network<-rbind(c(0, 1, 0, 0, 1),
               c(1, 0, 0, 0, 1),
               c(0, 0, 0, 0, 0),
               c(0, 0, 0, 0, 1),
               c(1, 1, 0, 1, 0))
rownames(network)<-c("a","b","c","d","e")			   
colnames(network)<-c("a","b","c","d","e")	

q01<-0.1
q10<-0.1
pSpec<-0.5

fitPhyloNetwork(phy, network, pars=c(0.1, 0.1, 0.7))

foo<-function(x) {
  q01<-exp(x[1])
  q10<-exp(x[2])
  pSpec<-exp(x[3])
  fitPhyloNetwork(phy, network, pars=c(q01, q10, pSpec))
}

foo(log(c(0.1, 0.1, 0.7)))

out<-optim(log(c(0.5, 0.5, 0.5)), foo, control=list(trace=6))
exp(out$par)

tt<-sim.bdtree(b=1, d=0, stop="taxa", n=60)
tt<-drop.tip(tt, "60")
xx<-simPhyloNetwork(tt, qRate=0.1, sProb=0.5)

fitPhyloNetwork(tt, xx, c(0.1, 0.1, 0.5))
foo<-function(x) {
  q01<-exp(x[1])
  q10<-exp(x[2])
  pSpec<-exp(x[3])
  fitPhyloNetwork(tt, xx, q01, q10, pSpec)
}

foo(log(c(0.1, 0.1, 0.5)))
out<-optim(log(c(0.5, 0.5, 0.5)), foo, control=list(trace=6))
exp(out$par)

lnlSurf<-matrix(nrow=9, ncol=9)
for(fr in 1:9)
  for(br in 1:9)
    lnlSurf[fr, br]<-fitPhyloNetwork(tt, xx, q01=fr/45, q10=br/45, pSpec=0.5)

contour(lnlSurf, levels=320:340, col="red", x=1:9/45, y=1:9/45)	 
contour(lnlSurf, , x=1:9/45, y=1:9/45, add=T)	 


# MCMC

nGen<-1000
parStart<-0.5
propWidth<-0.1
upperLimit<-1
lowerLimit<-0

lnL_save<-numeric(nGen)
lnL_save[]<-NA
frPost<-numeric(nGen)
brPost<-numeric(nGen)
pPost<-numeric(nGen)

frCurr<-parStart
brCurr<-parStart
pCurr<-parStart

lnL_save[1]<- - fitPhyloNetwork(tt, xx, q01= frCurr, q10= brCurr, pSpec= pCurr)
frPost[1]<-frCurr
brPost[1]<-brCurr
pPost[1]<-pCurr

plot(1:1000, lnL_save, pch=19, cex=0.3, ylim=c(-300, -500))
for(i in 2: nGen) {
  frProp<-frCurr+runif(1, min=-propWidth, max=propWidth)
  brProp<-brCurr+runif(1, min=-propWidth, max=propWidth)
  pProp<-pCurr+runif(1, min=-propWidth, max=propWidth)
  
  if(frProp > upperLimit) frProp <-upperLimit
  if(brProp > upperLimit) brProp <-upperLimit
  if(pProp > 1.0) pProp <-1.0
  
  if(frProp < lowerLimit) frProp <-lowerLimit
  if(brProp < lowerLimit) brProp <-lowerLimit
  if(pProp < lowerLimit) pProp <-lowerLimit
  
  
  
  lnLpc<-lnL_save[i-1]
  lnLpp<- - fitPhyloNetwork(tt, xx, q01= frProp, q10= brProp, pSpec= pProp)
  LR<-exp(lnLpp-lnLpc)
  rr<-runif(1)
  if(rr<LR) { # accept
    frCurr<-frProp
    brCurr<-brProp
    pCurr<-pProp
    lnL_save[i]<-lnLpp
  } else {
    lnL_save[i]<-lnLpc
  }
  frPost[i]<-frCurr
  brPost[i]<-brCurr
  pPost[i]<-pCurr
  points(i, lnL_save[i], cex=0.3, pch=19)
  cat(i, lnL_save[i], frCurr, brCurr, "\n")
}

library(lessR)

pdf("/Users/lukeh/Documents/nuismerGrant/heuristicModel/mcmcNetworkFit.pdf", width=8, height=8)
layout(matrix(1:4, nrow=2, ncol=2))
plot(lnL_save, type="l", lwd=3, xlab="Generation", ylab="lnL")
hist(pPost[-(1:100)], xlim=c(0,1), breaks=0:20/20, freq=F, xlab="P at speciation", main="")
lines(c(0, 0, 1,1), c(0, 1, 1, 0), col="grey")
lines(c(0.5, 0.5), c(0,10),  lwd=2, lty=2)
hist(frPost[-(1:100)], xlim=c(0,1), breaks=0:20/20, freq=F, xlab="Forward rate", main="")
lines(c(0, 0, 1,1), c(0, 1, 1, 0), col="grey")
lines(c(0.1, 0.1), c(0,100),  lwd=2, lty=2)

hist(brPost[-(1:100)], xlim=c(0,1), breaks=0:20/20, freq=F, xlab="Backward rate", main="")
lines(c(0, 0, 1,1), c(0, 1, 1, 0), col="grey")
lines(c(0.1, 0.1), c(0,100), lwd=2, lty=2)
dev.off()


library(geiger)
library(network)
phy<-birthdeath.tree(b=1, d=0, taxa.stop=41)
phy<-drop.tip(phy, "41")

nw<-simPhyloNetwork(phy, qRate=0.03, sProb=0.03)

nn<-network(nw, directed=F, vertex.attrnames=rownames(nw))

plot(phy)
plot.network(nn, label=network.vertex.names(nn), mode="circle")
plot.network(nn, label=network.vertex.names(nn), mode="fruchtermanreingold")
plot.network(nn, label=network.vertex.names(nn), mode="kamadakawai")

