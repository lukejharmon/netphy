library(geiger)
library(plyr)
library(doMC)
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

sym.trans <- function(pars) {
    n <- pars[, 1]
    q01 <- pars[, 2]
    q10 <- pars[, 3]
    pSpec <- pars[, 4]
    lambda <- pars[, 5]
    ntaxa <- pars[, 6]

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
    out <- optim(log(c(0.5, 0.5, 0.5)), foo, control = list(trace = 6))
    res <- c(exp(out$par), out$value)
    names(res) <- c("q01.fit", "q10.fit", "pSpec.fit", "llik.fit")
    return(res)
}

registerDoMC(50)
result.symtrans <- adply(.data = inputData, .margins = 1, .fun = sym.trans, .parallel = TRUE)
write.table(result.symtrans, "./results_symtrans.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

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
