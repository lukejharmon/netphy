\name{simPhyloNetwork}
\alias{simPhyloNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simPhyloNetwork(phy, qRate, sProb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
}
  \item{qRate}{
%%     ~~Describe \code{qRate} here~~
}
  \item{sProb}{
%%     ~~Describe \code{sProb} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (phy, qRate, sProb) 
{
    bt <- sort(branching.times(phy), decr = T)
    nb <- dim(phy$edge)[1]
    interactionMatrix <- matrix(nrow = nb, ncol = nb)
    interactionMatrix[] <- 0
    currentEdges <- numeric()
    rootNumber <- as.numeric(names(bt)[1])
    qMatrix <- rbind(c(-1, 1), c(1, -1)) * qRate
    for (i in 1:length(bt)) {
        thisEdge <- as.numeric(names(bt)[i])
        ancestorRow <- which(phy$edge[, 2] == thisEdge)
        descendantRow <- which(phy$edge[, 1] == thisEdge)
        if (length(ancestorRow) == 0) {
            r <- runif(1)
            if (r < sProb) {
                interactionMatrix[descendantRow, descendantRow] <- 1
                diag(interactionMatrix) <- 0
            }
            currentEdges <- c(currentEdges, descendantRow)
        }
        else {
            timeSpan <- bt[i - 1] - bt[i]
            tProb <- geiger:::MatrixExp.eig(qMatrix * timeSpan)
            cn <- length(currentEdges)
            for (j in 1:(cn - 1)) for (k in (j + 1):cn) {
                e1 <- currentEdges[j]
                e2 <- currentEdges[k]
                currState <- interactionMatrix[e1, e2]
                p0 <- tProb[currState + 1, 1]
                r <- runif(1)
                if (r < p0) {
                  interactionMatrix[e1, e2] <- 0
                  interactionMatrix[e2, e1] <- 0
                }
                else {
                  interactionMatrix[e1, e2] <- 1
                  interactionMatrix[e2, e1] <- 1
                }
            }
            toCut <- which(currentEdges == ancestorRow)
            currentEdges <- currentEdges[-toCut]
            currentEdges <- c(currentEdges, descendantRow)
            for (j in 1:length(descendantRow)) {
                interactionMatrix[, descendantRow[j]] <- interactionMatrix[, 
                  ancestorRow]
                interactionMatrix[descendantRow[j], ] <- interactionMatrix[ancestorRow, 
                  ]
                diag(interactionMatrix) <- 0
            }
            r <- runif(1)
            if (r < sProb) {
                interactionMatrix[descendantRow, descendantRow] <- 1
            }
        }
    }
    nTaxa <- length(phy$tip.label)
    tips <- which(phy$edge[, 2] <= nTaxa)
    oo <- phy$edge[tips, 2]
    res <- interactionMatrix[tips, tips]
    rownames(res) <- phy$tip.label[oo]
    colnames(res) <- phy$tip.label[oo]
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
