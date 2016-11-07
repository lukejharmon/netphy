library(geiger)
library(plyr)
library(doMC)
library(foreach)
library(dplyr)
lapply(paste0("./R/",as.list(list.files("./R/"))), source)

## sProb <- c(0,2)  # pSpec = probability of speciation within the network
## qProb <- c(0,1)  # transition probability between interacting and not-interacting
## tProb <- c(10,100)  # total number of taxa in the tree
## spProb <- c(0,2)  # lambda for tree simulation

##### Symmetrical transitions

## tree.num <- 1:1000
## pSpec <- runif(1000, sProb[1], sProb[2])
## q01 <- runif(1000, qProb[1], qProb[2])
## lambda <- runif(1000, spProb[1], spProb[2])
## ntaxa <- round(runif(1000, tProb[1], tProb[2]))
## trees <- as.list(rep(NA, 1000))
## nets <- as.list(rep(NA, 1000))

## inputData <- data.frame(tree.num = tree.num, q01.sim = q01, q10.sim = q01, pSpec.sim = pSpec, lambda.sim = lambda, ntaxa = ntaxa)

## write.table(inputData, "./inputData.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
inputDataSym <- read.csv("./inputData.csv")

registerDoMC(56)


cat(c("tree.number", "q01.sim", "q10.sim", "pSpec.sim", "lambda.sim", "ntaxa.sim", "q01.fit", "q10.fit", "pSpec.fit", "llik.fit"), "\n", append = FALSE, file = "results_symtrans.txt", sep = "\t")

result.symtrans <- foreach(i = 1:dim(inputDataSym)[1]) %dopar% {
    n <- inputDataSym$tree.num[i]
    q01 <- inputDataSym$q01.sim[i]
    q10 <- inputDataSym$q10.sim[i]
    pSpec <- inputDataSym$pSpec.sim[i]
    lambda <- inputDataSym$lambda.sim[i]
    ntaxa <- inputDataSym$ntaxa[i]

    tt <- sim.bdtree(b = lambda, d = 0, stop = "taxa", n = ntaxa)
    tt <- drop.tip(tt, paste0("s", ntaxa))
    write.tree(tt, file = paste0("./trees/symtrans/tree_", n, "_symtrans.txt"))
    nets <- simPhyloNetwork(tt, qRate = q01, sProb = pSpec)
    write.table(nets, file = paste0("./nets/symtrans/net_", n, "_symtrans.csv"), sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)

    foo <- function(par) {
        q01 <- exp(par[1])
        q10 <- exp(par[2])
        pSpec <- exp(par[3])
        fitPhyloNetwork(tt, nets, c(q01, q10, pSpec))
    }

    boxconstrain=function (f, lower, upper, fail.value)
	{
		function(x) {
			if (any(x < lower | x > upper)) fail.value else f(x)
		}
	}
    f = boxconstrain(foo, -10, 10, fail.value=1e200)


    out <- tryCatch(optim(log(c(q01, q10, pSpec)), f, control = list(trace = 6, , reltol = .Machine$double.eps^.25)), error = function(x) return(data,frame(q01.fit = NA, q10.fit = NA, pSpec.fit = NA, llik.fit = NA)))
    data.frame(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)
    cat(sprintf("%s\t", c(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)), "\n", file = "./results_symtrans_boxconstrain.txt", append = TRUE)
}


##### Asymmetrical transitions

## tree.num <- 1:1000
## pSpec <- runif(1000, sProb[1], sProb[2])
## q01 <- runif(1000, qProb[1], qProb[2])
## q10 <- runif(1000, qProb[1], qProb[2])
## lambda <- runif(1000, spProb[1], spProb[2])
## ntaxa <- round(runif(1000, tProb[1], tProb[2]))
## trees <- as.list(rep(NA, 1000))
## nets <- as.list(rep(NA, 1000))

## inputDataAsym <- data.frame(tree.num = tree.num, q01.sim = q01, q10.sim = q10, pSpec.sim = pSpec, lambda.sim = lambda, ntaxa = ntaxa)

## write.table(inputDataAsym, "./inputDataAsym.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")

## inputDataAsym <- read.csv("./inputDataAsym.csv")

## registerDoMC(56)

## cat(c("tree.number", "q01.sim", "q10.sim", "pSpec.sim", "lambda.sim", "ntaxa.sim", "q01.fit", "q10.fit", "pSpec.fit", "llik.fit"), "\n", append = FALSE, file = "results_asymtrans.txt", sep = "\t")

## result.asymtrans <- foreach(i = 1:dim(inputDataAsym)[1]) %dopar% {
##     n <- inputDataAsym$tree.num[i]
##     q01 <- inputDataAsym$q01.sim[i]
##     q10 <- inputDataAsym$q10.sim[i]
##     pSpec <- inputDataAsym$pSpec.sim[i]
##     lambda <- inputDataAsym$lambda.sim[i]
##     ntaxa <- inputDataAsym$ntaxa[i]

##     tt <- sim.bdtree(b = lambda, d = 0, stop = "taxa", n = ntaxa)
##     tt <- drop.tip(tt, paste0("s", ntaxa))
##     write.tree(tt, file = paste0("./trees/asymtrans/tree_", n, "_asymtrans.txt"))
##     nets <- simPhyloNetwork(tt, qRate = q01, sProb = pSpec)
##     write.table(nets, file = paste0("./nets/asymtrans/net_", n, "_asymtrans.csv"), sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)

##     foo <- function(par) {
##         q01 <- exp(par[1])
##         q10 <- exp(par[2])
##         pSpec <- exp(par[3])
##         fitPhyloNetwork(tt, nets, c(q01, q10, pSpec))
##     }

##     out <- tryCatch(optim(log(c(0.5, 0.5, 0.5)), foo, control = list(trace = 6)), error = function(x) return(data,frame(q01.fit = NA, q10.fit = NA, pSpec.fit = NA, llik.fit = NA)))
##     data.frame(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)
##     cat(sprintf("%s\t", c(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)), "\n", file = "./results_asymtrans.txt", append = TRUE)
## }
