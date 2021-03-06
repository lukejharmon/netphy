#### Tweaking optimization parameters to get better estimates

### Target tree is #100

source("figs.R")
lapply(paste0("./R/",as.list(list.files("./R/"))), source)

inputDataSym <- read.csv("./inputData.csv")
i <- 100

n <- inputDataSym$tree.num[i]
q01 <- inputDataSym$q01.sim[i]
q10 <- inputDataSym$q10.sim[i]
pSpec <- inputDataSym$pSpec.sim[i]
lambda <- inputDataSym$lambda.sim[i]
ntaxa <- inputDataSym$ntaxa[i]

tt <- read.tree("./trees/symtrans/tree_100_symtrans.txt")
#tt <- drop.tip(tt, paste0("s", ntaxa))
#write.tree(tt, file = paste0("./trees/symtrans/tree_", n, "_symtrans.txt"))
nets <- read.csv("./nets/symtrans/net_100_symtrans.csv")
#write.table(nets, file = paste0("./nets/symtrans/net_", n, "_symtrans.csv"), sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)

foo <- function(par) {
    q01 <- exp(par[1])
    q10 <- exp(par[2])
    pSpec <- exp(par[3])
    fitPhyloNetwork(tt, nets, c(q01, q10, pSpec))
}

minTrans<--500
maxTrans<- log(100*Nedge(tt))-log(sum(tt$edge.length))
mn=rep(0, 3)
mx=rep(1, 3)
out <- tryCatch(optim(log(c(0.5, 0.5, 0.5)), foo, method = "L-BFGS-B", control = list(trace = 6), lower = mn, upper = mx))#, error = function(x) return(data,frame(q01.fit = NA, q10.fit = NA, pSpec.fit = NA, llik.fit = NA)))
data.frame(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)



### Box constraining
foo.constr <- function(par) {
    ll <- NA
    q01 <- exp(par[1])
    q10 <- exp(par[2])
    pSpec <- exp(par[3])
    if(par[1] > -5 | par[1] < 5 | par[2] > -5 | par[2] < 5 | par[3] > -5 | par[3] < 5)
    {
        ll <- fitPhyloNetwork(tt, nets, c(q01, q10, pSpec))
    }
    ll
}

minTrans<--500
	maxTrans<- log(100*Nedge(tt))-log(sum(tt$edge.length))
	mn=rep(-10, 3)
	mx=rep(10, 3)
out <- optim(log(c(0.5, 0.5, 0.5)), foo.constr, control = list(trace = 6))#, error = function(x) return(data,frame(q01.fit = NA, q10.fit = NA, pSpec.fit = NA, llik.fit = NA)))
data.frame(n, q01, q10, pSpec, lambda, ntaxa, exp(out$par[1]), exp(out$par[2]), exp(out$par[3]), out$value)


### Starting with simulated values
## NEEDS TO BE DONE




### Testing boxconstrain from diversitree - based on geiger's traits.R line 691

boxconstrain=function (f, lower, upper, fail.value)
	{
		function(x) {
			if (any(x < lower | x > upper)) fail.value else f(x)
		}
	}
f = boxconstrain(foo, -100, 100, fail.value=1e200)

out.box2 <- optim(log(c(q01, q10, pSpec)), f, control = list(trace = 6, reltol = .Machine$double.eps^.25, parscale = rep(0.1, 3)))
