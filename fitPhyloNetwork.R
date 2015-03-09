fitPhyloNetwork<-function(phy, network, q01, q10, pSpec) {
	if(sum(phy$edge.length==0)>0) cat("Function probably will not work when tree has zero-length branches\n\n")
	
	mm<-match(phy$tip.label, rownames(network))
	network<-network[mm,mm]
	nNodes<-max(phy$edge)
	nTips<-length(phy$tip.label)
	lik<-list()
	
	# should deal with this
	if(pSpec>1) pSpec<-1
	
	qMatrix<-rbind(c(-q01, q01), c(q10, -q10))
	
	# make tip matrices
	for(i in 1:nTips) {
		theTip<-phy$tip.label[i]
		theRow<-which(rownames(network)==theTip)
		interact<-network[theRow,]
		rr<-matrix(nrow=nTips, ncol=2)
		rownames(rr)<-phy$tip.label
		colnames(rr)<-c(0,1)
		
		for(j in 1:nTips) 
			if(interact[j]==0) rr[j,]<-c(0,-Inf) else rr[j,]<-c(-Inf,0) 
		
		lik[[i]]<-rr
		names(lik)[i]<-theTip	
		
		
	}
	
	currTree<-phy
	cumlnL<-0
	# begin algorithm

	while(1) {	
		# get most recent node and daughter species
		bt<-branching.times(currTree)
		recentNode<-names(bt)[which(bt==min(bt))][1]
		daughter<-node.leaves(currTree,recentNode)
		timeInterval<-bt[recentNode]
		
		# Calculations for merging species
		d1<-which(names(lik)==daughter[1])
		d1L<-lik[[d1]]
		d2<-which(names(lik)==daughter[2])
		d2L<-lik[[d2]]
		
		dCon<-d1L[daughter[2],]
		# should be the same!
		#d2L[daughter[1],]
		
		# likelihood of states at node
		dNode<-branchLike(dCon, timeInterval, qMatrix)
		# times prob dist at speciation
		dd<-dNode + log(c(1-pSpec, pSpec))
		# add to running total of likelihood
		cumlnL<-cumlnL + -logspace_sum(dd)
			
		if(length(currTree$tip.label)==2) break;
		
		# prune out one species
		currTree<-drop.tip(currTree, daughter[1])
			
		oldLik<-lik
		lik<-list()
		
		# shorten all other branches and recalculate states
		nn<-length(currTree$tip.label)
		for(i in 1: nn) {
			theTip<-currTree$tip.label[i]
			rr<-matrix(nrow=nn, ncol=2)
			rownames(rr)<-currTree$tip.label
			colnames(rr)<-c(0,1)
			
			# either this is one of the two that were merged
			if(theTip %in% daughter) {
				# get both interaction matrices
				xx<-which(names(oldLik)==daughter[1])
				ol1<-oldLik[[xx]]
				xx<-which(names(oldLik)==daughter[2])
				ol2<-oldLik[[xx]]
				
				for(j in 1:nrow(ol1)) {
					# daughter 1 has been pruned and is redundant
					if(rownames(ol1)[j] != daughter[1]) {
						l1<-branchLike(ol1[j,], timeInterval, qMatrix)
						l2<-branchLike(ol2[j,], timeInterval, qMatrix)
						rr[rownames(ol1)[j],]<-l1+l2
						
					}
					
				}
							
			} else { # or not
				xx<-which(names(oldLik)==theTip)
				theOldLik<-oldLik[[xx]]
				
				# deal with merged thingies
				lMerg<-theOldLik[daughter,]
				l1<-branchLike(lMerg[1,], timeInterval, qMatrix)
				l2<-branchLike(lMerg[2,], timeInterval, qMatrix)
				rr[daughter[2],]<-l1+l2
				
				#deal with not-merged thingies
				theRest<-rownames(rr)[which(!(rownames(rr) %in% daughter))]
				for(tt in theRest) {
					lold<-theOldLik[tt,]
					lnew<-branchLike(lold, timeInterval, qMatrix)
					rr[tt,]<-lnew
				}
			}
			
			rr[theTip,]<-c(0,-Inf)
			lik[[i]]<-rr
			names(lik)[i]<-theTip
			nEdge<-which(currTree$edge[,2]==i)
			currTree$edge.length[nEdge]<-currTree$edge.length[nEdge]-timeInterval
				
				
			}
	
			if(length(lik)==80) break;
	
	}	
	
	return(cumlnL)	
	
}

library(geiger)
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

fitPhyloNetwork(phy, network, q01=0.1, q10=0.1, pSpec=0.7)

foo<-function(x) {
	q01<-exp(x[1])
	q10<-exp(x[2])
	pSpec<-exp(x[3])
	fitPhyloNetwork(phy, network, q01, q10, pSpec)
}

foo(log(c(0.1, 0.1, 0.7)))

out<-optim(log(c(0.5, 0.5, 0.5)), foo, control=list(trace=6))
exp(out$par)

tt<-birthdeath.tree(b=1, d=0, taxa.stop=60)
tt<-drop.tip(tt, "60")
xx<-simPhyloNetwork(tt, qRate=0.1, sProb=0.5)
fitPhyloNetwork(tt, xx, q01=0.1, q10=0.1, pSpec=0.5)
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