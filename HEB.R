library(reshape2)
library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(ggradar)
library(ggpubr)

data=read.table("../../allgenes_noHE.count",header = 1,row.names = 1)
data$geneA=rownames(data)
data$geneD=rownames(data)

homoeo=read.table("../homoeologs_expr_9642.txt")
colnames(homoeo)=c("geneA","geneD")

countsA=merge(homoeo,data,by="geneA")
colnames(countsA)[2:56]=c("geneD",paste0("A_",colnames(countsA)[3:56]))
counts=merge(countsA,data,by="geneD")
counts=counts[c(2,1,3:56,58:111)]
colnames(counts)[c(1,57:110)]=c("geneA",paste0("D_",colnames(counts)[57:110]))
##
len=read.table("../AADD_longest_exon.length")

colnames(len)=c("geneA","len")
lenA=data.frame(merge(counts,len,by="geneA",sort=F)[c(1,111)],row.names = 1)

colnames(len)=c("geneD","len")
lenD=data.frame(merge(counts,len,by="geneD",sort=F)[c(1,111)],row.names = 1)

deg=list()
for(i in seq(3,56,3)){
  dd=counts[,c(i:(i+2),(i+54):(i+56))]
  rownames(dd)=counts$geneA
  condition=factor(sub("_\\d+","",colnames(dd)))
  n=sub("A_|D_","",condition[1])
  
  coldata=data.frame(row.names = colnames(dd),condition)
  dds=DESeqDataSetFromMatrix(countData = dd,colData = coldata,design = ~condition)
  dds=estimateSizeFactors(dds)
  sf=as.numeric(sizeFactors(dds))
  
  normfactors_A1=lenA*sf[1]/1000
  normfactors_A2=lenA*sf[2]/1000
  normfactors_A3=lenA*sf[3]/1000
  normfactors_D1=lenD*sf[4]/1000
  normfactors_D2=lenD*sf[5]/1000
  normfactors_D3=lenD*sf[6]/1000
  normFactors = as.matrix(cbind(normfactors_A1,normfactors_A2,normfactors_A3,normfactors_D1,normfactors_D2,normfactors_D3))
  normFactors = normFactors / exp(rowMeans(log(normFactors)))
  normalizationFactors(dds) = normFactors
  
  print(n)
  DE=DESeq(dds)
  res=as.data.frame(results(DE,lfcThreshold=1,alpha = 0.05))
  
  for(j in 1:nrow(res)){
    if(!is.na(res[j,2])&!is.na(res[j,6])){
      if(res[j,6]<0.05&res[j,2]>=1&res[j,2]<2){
        res[j,7]="Dbias_2fold"
      }else if(res[j,6]<0.05&res[j,2]>=2){
        res[j,7]="Dbias_4fold"
      }else if(res[j,6]<0.05&res[j,2]<(-1)&res[j,2]>(-2)){
        res[j,7]="Abias_2fold"
      }else if(res[j,6]<0.05&res[j,2]<=(-2)){
        res[j,7]="Abias_4fold"
      }else{
        res[j,7]="equal"
      }
    }else{
      res[j,7]="Sil"
    }
  }
  
  colnames(res)[7]=n
  res$gene=rownames(res)
  
  deg[[i/3]]=res[c(8,7)]
}

bias=Reduce(merge,deg)
write.table(bias,"bias_deseq_strict.txt",sep = "\t",quote = F,row.names = F)