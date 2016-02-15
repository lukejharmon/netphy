setwd("/Users/lukeharmon/Documents/networksOnTreesPaper/netphy/manuscript/networkPropertySims/data/")

dfiles<-dir()
theFoodWebs<-grep("fw", dfiles, value=T)
fw_data<-list()
for(i in 1:length(theFoodWebs)) {
  fw_data[[i]]<-read.table(theFoodWebs[i], sep="\t", row.names=1, header=T)
}

library(igraph)

centr_fw<-numeric(length(theFoodWebs))

pdf("../allFwPlots.pdf")
for(i in 1:length(theFoodWebs)) {
  gg<-graph.adjacency(as.matrix(fw_data[[i]]))
  plot(gg)
  centr_fw[i]<-centr_degree(gg)$centralization
}
dev.off()

hist(centr_fw)


