#' Fit a model of network evolution to network data and a phylogenetic tree
#'
#' This function returns the likelihood
#'
#' @param phy An object of class 'phylo'
#' @param network A network, specified by an n x n matrix filled with 0s and 1s 
#' @param pars Model parameters: q01, q10, pSpec
#' @export
fitPhyloNetwork<-function(phy, network, pars) {
	if(sum(phy$edge.length==0)>0) cat("Function probably will not work when tree has zero-length branches\n\n")
	
	mm<-match(phy$tip.label, rownames(network))
	network<-network[mm,mm]
	nNodes<-max(phy$edge)
	nTips<-length(phy$tip.label)
	lik<-list()
	
  q01<-pars[1]
  q10<-pars[2]
  pSpec<-pars[3]
  
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
	
} #
