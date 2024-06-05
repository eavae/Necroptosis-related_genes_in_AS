

library(limma)
library(pheatmap)

expFile="FRGexp.txt"     
setwd("C:\\Users\\ghs\\Desktop\\07.diff")    


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
conNum=length(Type[Type=="Control"])
treatNum=length(Type[Type=="Treat"])


sigVec=c()
sigGeneVec=c()
diffTab=data.frame()
for(i in row.names(data)){
	test=wilcox.test(data[i,] ~ Type)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	conMean=mean(data[i,1:conNum])
	treatMean=mean(data[i,(conNum+1):ncol(data)])
	logFC=treatMean-conMean
	if(pvalue<0.05){
		sigVec=c(sigVec, paste0(i, Sig))
		sigGeneVec=c(sigGeneVec, i)
		if(logFC>0){diffTab=rbind(diffTab, cbind(Gene=i,conMean,treatMean,pvalue,Type="Up"))}
		if(logFC<0){diffTab=rbind(diffTab, cbind(Gene=i,conMean,treatMean,pvalue,Type="Down"))}
	}
}

write.table(diffTab, file="diff.txt", sep="\t", quote=F, row.names=F)

data=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec


names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

