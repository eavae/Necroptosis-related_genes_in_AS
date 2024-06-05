

library(limma)     
expFile="normalize.txt"    
geneFile="gene.txt"        
setwd("C:\\Users\\ghs\\Desktop\\06.FRGexp")    


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)


gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]


out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="NCGexp.txt", sep="\t", quote=F, col.names=F)
