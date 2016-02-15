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

for(i in 4392:length(simStatistics)) {
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

library(ggplot2)
type<-rep(c("data", "sim"), times=c(length(theFoodWebs), 100))
rr<-rbind(networkStatistics, simStatistics)
allData<-data.frame(type, rr)


pdf("../comparisons.pdf")
plot(allData[,c("type", "density")])
plot(allData[,c("type", "avg_degree")])
plot(allData[,c("type", "avg_path_length")])
plot(allData[,c("type", "clustering_coefficient")])
plot(allData[,c("type", "cohesion")])
plot(allData[,c("type", "centr_degree")])
dev.off()
