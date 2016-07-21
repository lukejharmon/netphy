library(geiger)
library(plyr)
library(doMC)
library(foreach)
library(dplyr)
lapply(paste0("./R/",as.list(list.files("./R/"))), source)

sProb <- c(0,2)  # pSpec = probability of speciation within the network
qProb <- c(0,1)  # transition probability between interacting and not-interacting
tProb <- c(10,100)  # total number of taxa in the tree
spProb <- c(0,2)  # lambda for tree simulation

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
inputData <- read.csv("./inputData.csv")

registerDoMC(56)
:dim(inputData)[1]

cat(c("tree.number", "q01.sim", "q10.sim", "pSpec.sim", "lambda.sim", "ntaxa.sim", "q01.fit", "q10.fit", "pSpec.fit", "llik.fit"), "\n", append = FALSE, file = "results_symtrans.txt", sep = "\t")

result.symtrans <- foreach(i = c(7,11)) %dopar% {
    n <- inputData$tree.num[i]
    q01 <- inputData$q01.sim[i]
    q10 <- inputData$q10.sim[i]
    pSpec <- inputData$pSpec.sim[i]
    lambda <- inputData$lambda.sim[i]
    ntaxa <- inputData$ntaxa[i]

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

    #foo(log(c(q01[x], q10[x], pSpec[x])))
    out <- tryCatch(optim(log(c(0.5, 0.5, 0.5)), foo, control = list(trace = 6)), error = function(x) return(data,frame(q01.fit = NA, q10.fit = NA, pSpec.fit = NA, llik.fit = NA)))
    data.frame(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)
    cat(sprintf("%s\t", c(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)), "\n", file = "./results_symtrans.txt", append = TRUE)
    #names(res) <- c("q01.fit", "q10.fit", "pSpec.fit", "llik.fit")
    #return(res)
}

#result.symtrans <- adply(.data = inputData, .margins = 1, .fun = sym.trans, .parallel = TRUE)
#results.symtrans <- bind_rows(result.symtrans)
#names(results.symtrans) <- c("tree.number", "q01.sim", "q10.sim", "pSpec.sim", "lambda.sim", "ntaxa.sim", "q01.fit", "q10.fit", "pSpec.fit", "llik.fit")
#write.table(result.symtrans, "./results_symtrans.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

##### Asymmetrical transitions

## pSpec <- runif(1, sProb[1], sProb[2])
## q01 <- runif(1, qProb[1], qProb[2])
## q10 <- runif(1, qProb[1], qProb[2])
## lambda <- runif(1, spProb[1], spProb[2])
## ntaxa <- round(runif(1, tProb[1], tProb[2]))


## tt <- sim.bdtree(b = lambda, d = 0, stop = "taxa", n = ntaxa)
## tt <- drop.tip(tt, paste0("s", ntaxa))
## net <- simPhyloNetwork(tt, qRate = q01, sProb = pSpec)


## foo <- function(par) {
##     q01 <- exp(par[1])
##     q10 <- exp(par[2])
##     pSpec <- exp(par[3])
##     fitPhyloNetwork(tt, net, c(median(qProb), median(qProb), median(sProb)))
## }

## foo(log(c(q01, q10, pSpec)))
## out <- optim(log(c(0.5, 0.5, 0.5)), foo, control = list(trace = 6))
## exp(out$par)
