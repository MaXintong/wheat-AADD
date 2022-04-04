tpm=read.table("../AADD_TPM.txt",header = T,row.names = 1)
tpm=tpm[,c(28:36,52:54,37:51,1:9,25:27,10:24)]
tpm$gene=rownames(tpm)

homoeos=read.table("../homoeologs_expr_9642.txt")
homoe=data.frame(gene=matrix(ncol = 1,nrow = 9642*2))

homoe[1:9642,]=homoeos[1:9642,1]
homoe[9643:19284,]=homoeos[1:9642,2]

expr=data.frame(merge(tpm,homoe,by="gene",sort = F),row.names = 1)
expr=expr[apply(expr,1,var)>0.5,]

datExprT=t(expr)
data_mix=datExprT[grep("mix",rownames(datExprT)),]
data_at2=datExprT[grep("AT2",rownames(datExprT)),]

nSets = 2
setLabels = c("mix" ,"at2")
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = data_mix);
multiExpr[[2]] = list(data = data_at2);

pdf("SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off()
save(multiExpr,nSets, setLabels, file = "dataInput.RData")

library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
load("dataInput.RData")

powers = c(c(1:10), seq(from = 12, to=40, by=2))
pickSoftThreshold(multiExpr[[1]]$data, powerVector=powers, verbose = 2, networkType = "signed")

power=18
for(set in 1:nSets){
  subData=multiExpr[[set]]$data
  subData=apply(subData,2,as.numeric)
  power=power
  net=blockwiseModules(subData, 
                       power = power, blocks = NULL, maxBlockSize = 17278,
                       networkType = "signed",TOMType = "signed", corType = "pearson",
                       deepSplit = 2, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,saveTOMFileBase=paste0(setLabels[set],"_TOM"),
                       verbose = 3)
  mergedColors = labels2colors(net$colors)
  pdf(paste0(setLabels[set],"_modules.pdf"),width = 10,height = 5)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$goodGenes],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  assign(paste0(setLabels[set],"_net"),net)
}
save(list = grep(".+net",ls(),value = T),file = "individual-Network.RData")