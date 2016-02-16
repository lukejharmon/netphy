setwd("/Users/lukeharmon/Documents/networksOnTreesPaper/netphy/manuscriptAnalyses/networkPropertySims/data/")

dfiles<-dir()
theFoodWebs<-grep("fw", dfiles, value=T)
fw_data<-list()
for(i in 1:length(theFoodWebs)) {
  fw_data[[i]]<-read.table(theFoodWebs[i], sep="\t", row.names=1, header=T)
}

library(igraph)

networkStatistics<-matrix(nrow=length(theFoodWebs), ncol=7)
colnames(networkStatistics)<-c("nnode", "density", "avg_degree", "avg_path_length", "clustering_coefficient", "cohesion", "centr_degree")
rownames(networkStatistics)<-theFoodWebs

pdf("../allFwPlots.pdf")
for(i in 1:length(theFoodWebs)) {
  gg<-graph.adjacency(as.matrix(fw_data[[i]]))
  plot(gg)
  networkStatistics[i,1]<-vcount(gg)
  networkStatistics[i,2]<-graph.density(gg)
  networkStatistics[i,3]<-mean(degree(gg))
  networkStatistics[i,4]<-average.path.length(gg)
  networkStatistics[i,5]<-transitivity(gg)
  networkStatistics[i,6]<-graph.cohesion(gg)
  networkStatistics[i,7]<-centr_degree(gg)$centralization
}
dev.off()

hist(networkStatistics[,1])

library(geiger)


#library(devtools)
#install_github("lukejharmon/netphy")
library(np)

treeSize<-50
qq<-1:10/10
sp<-1:10/10
simStatistics<-matrix(nrow=100*10*10, ncol=9)
colnames(simStatistics)<-c("qRate", "sProb", "nnode", "density", "avg_degree", "avg_path_length", "clustering_coefficient", "cohesion", "centr_degree")
simStatistics[,1]<-rep(qq, each=100, 10)
simStatistics[,2]<-rep(sp, each=1000)

for(i in 1:10000) {
  tt<-birthdeath.tree(b=1, d=0, taxa.stop=50)
  tt<-drop.tip(tt, "50")
  xx<-simPhyloNetwork(tt, qRate=simStatistics[i,1], sProb=simStatistics[i,2])
  gg<-graph.adjacency(xx)
  simStatistics[i,3]<-vcount(gg)
  simStatistics[i,4]<-graph.density(gg)
  simStatistics[i,5]<-mean(degree(gg))
  simStatistics[i,6]<-average.path.length(gg)
  simStatistics[i,7]<-transitivity(gg)
  simStatistics[i,8]<-graph.cohesion(gg)
  simStatistics[i,9]<-centr_degree(gg)$centralization
  cat(i, "\n")
}

colnames(simStatistics)
plot(simStatistics[,c(2,4)])

pdf("../density_sims.pdf", width=24, height=12)
layout(matrix(1:10, nrow=2, ncol=5, byrow=T))
for(i in 1:10) {
  rr<-which(simStatistics[,2]==sp[i])
  xx<-as.factor(c(simStatistics[rr,1], rep("data", length(theFoodWebs))))
  yy<-c(simStatistics[rr,4], networkStatistics[,2])
  plot(xx,yy,main=paste("sp = ", sp[i]), xlab="q", ylab="density")
}
dev.off()

pdf("../avg_degree_sims.pdf", width=12, height=6)
layout(matrix(1:10, nrow=2, ncol=5, byrow=T))
for(i in 1:10) {
  rr<-which(simStatistics[,2]==sp[i])
  xx<-as.factor(c(simStatistics[rr,1], rep("data", length(theFoodWebs))))
  yy<-c(simStatistics[rr,5], networkStatistics[,3])
  plot(xx,yy,main=paste("sp = ", sp[i]), xlab="q", ylab="degree")
}
dev.off()

pdf("../avg_pathlength_sims.pdf", width=12, height=6)
layout(matrix(1:10, nrow=2, ncol=5, byrow=T))
for(i in 1:10) {
  rr<-which(simStatistics[,2]==sp[i])
  xx<-as.factor(c(simStatistics[rr,1], rep("data", length(theFoodWebs))))
  yy<-c(simStatistics[rr,6], networkStatistics[,4])
  plot(xx,yy,main=paste("sp = ", sp[i]), xlab="q", ylab="path length")
}
dev.off()


pdf("../avg_clustering_coefficient_sims.pdf", width=12, height=6)
layout(matrix(1:10, nrow=2, ncol=5, byrow=T))
for(i in 1:10) {
  rr<-which(simStatistics[,2]==sp[i])
  xx<-as.factor(c(simStatistics[rr,1], rep("data", length(theFoodWebs))))
  yy<-c(simStatistics[rr,7], networkStatistics[,5])
  plot(xx,yy,main=paste("sp = ", sp[i]), xlab="q", ylab="clustering_coefficient")
}
dev.off()

pdf("../graphcohesion_sims.pdf", width=12, height=6)
layout(matrix(1:10, nrow=2, ncol=5, byrow=T))
for(i in 1:10) {
  rr<-which(simStatistics[,2]==sp[i])
  xx<-as.factor(c(simStatistics[rr,1], rep("data", length(theFoodWebs))))
  yy<-c(simStatistics[rr,8], networkStatistics[,6])
  plot(xx,yy,main=paste("sp = ", sp[i]), xlab="q", ylab="graph cohesion")
}
dev.off()

pdf("../centralization_sims.pdf", width=12, height=6)
layout(matrix(1:10, nrow=2, ncol=5, byrow=T))
for(i in 1:10) {
  rr<-which(simStatistics[,2]==sp[i])
  xx<-as.factor(c(simStatistics[rr,1], rep("data", length(theFoodWebs))))
  yy<-c(simStatistics[rr,9], networkStatistics[,6])
  plot(xx,yy,main=paste("sp = ", sp[i]), xlab="q", ylab="centralization")
}
dev.off()

