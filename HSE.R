data=read.table("../AADD_homoeolog.tpm",header = T,row.names = 2)[-1]

df=data.frame(matrix(nrow = nrow(data),ncol = 27))
for(i in 1:nrow(df)){
  for(j in seq(1,27,3)){
    tA=data[i,j:(j+2)]
    mA=data[i,(j+27):(j+29)]
    tD=data[i,(j+54):(j+56)]
    mD=data[i,(j+81):(j+83)]
    lgmA_mD=log2((mA+0.01)/(mD+0.01))
    lgtA_tD=log2((tA+0.01)/(tD+0.01))
    
    P=t.test(mA,mD)$p.value
    H=t.test(tA,tD)$p.value
    T=t.test(lgmA_mD,lgtA_tD)$p.value
    lgmix=t.test(lgmA_mD,lgtA_tD)$estimate[1]
    lgat=t.test(lgmA_mD,lgtA_tD)$estimate[2]
    
    if(is.na(P)==FALSE&is.na(H)==FALSE&is.na(T)==FALSE){
      if(P<0.05&H<0.05&T>0.05){
        df[i,j]="cis only"
      }else if(P<0.05&H>0.05&T<0.05){
        df[i,j]="trans only"
      }else if(P<0.05&H<0.05&T<0.05){
        if(lgmix>0&lgat>0){
          df[i,j]="cis + trans"
        }else if(lgmix<0&lgat<0){
          df[i,j]="cis + trans"
        }else if(lgmix>0&lgat<0){
          df[i,j]="cis X trans"
        }else if(lgmix<0&lgat>0){
          df[i,j]="cis X trans"
        }
      }else if(P>0.05&H<0.05&T<0.05){
        df[i,j]="compensatory"
      }else if(T>0.05&(H>0.05|P>0.05)){
        df[i,j]="conserved"
      }else{
        df[i,j]="ambiguous"
      }
    }else{
      df[i,j]="ambiguous"
    }
  }
}

cisTrans=df[seq(1,27,3)]
rownames(cisTrans)=rownames(data)
colnames(cisTrans)=c("leaf","pri.root","root","spike0.5","spike1.0","spike1.5","spike2.0","spike2.5","shoot")
cisTrans=cisTrans[c(1:3,9,4:8)]

write.table(cisTrans,"cisTrans.txt",sep = "\t",quote = F)