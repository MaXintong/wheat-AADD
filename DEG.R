library(DESeq2)

count=read.table("../allgenes.count",header = T,row.names = 1)
coldata<-data.frame( sample = factor(sub("_\\d+","", names(count))))

dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)

batch=rbind(
  c("AT2_lf","mix_lf"),
  c("AT2pri","mix_pri"),
  c("AT2_rt","mix_rt"),
  c("AT2_st","mix_st"),
  
  c("AT2s05","mix_s05"),
  c("AT2s10","mix_s10"),
  c("AT2s15","mix_s15"),
  c("AT2s20","mix_s20"),
  
  c("AT2s25","mix_s25")
)
pairwiseDE<-function(dds, contrast,savePath){

  print(contrast)
  ddsPW <-dds[,dds$sample %in% contrast]
  ddsPW$sample<-droplevels(ddsPW$sample)
  res <- results(DESeq(ddsPW))
  print( summary(res,alpha=.05) ) 
  write.table(res,quote = F, file=paste0(paste0(contrast, collapse="vs"),".txt"))
}  
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = "G:\\wheat_AADD_wgcna\\DEG"))