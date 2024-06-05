
library(ConsensusClusterPlus)     
expFile="diffGeneExp.txt"          
workDir=""    
setwd(workDir)       


data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)


group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="treat"]


maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="png")



clusterNum=2        
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("m6Acluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$m6Acluster))
cluster$m6Acluster=letter[match(cluster$m6Acluster, uniqClu)]
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="m6Acluster.txt", sep="\t", quote=F, col.names=F)



library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

gene="TARDBP"                          
clusterFile="m6Acluster.txt"         
ssgseaFile="CIBERSORT-Results.txt"       



cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)


ssgsea=read.table(ssgseaFile, header=T, sep="\t", check.names=F, row.names=1)
ssgsea=t(ssgsea)


sameSample=intersect(row.names(cluster), row.names(ssgsea))
cluster=cluster[sameSample,gene,drop=F]
ssgsea=ssgsea[sameSample,,drop=F]
rt=cbind(ssgsea, cluster)
rt[,gene]=ifelse(rt[,gene]>median(rt[,gene]), "High", "Low")
rt[,gene]=factor(rt[,gene], levels=c("Low", "High"))

data=melt(rt, id.vars=c(gene))
colnames(data)=c("Gene", "Immune", "Fraction")


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"Gene"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="Gene",
     xlab="",
     ylab="Immune infiltration",
     legend.title=gene,
     palette=bioCol)
p=p+rotate_x_text(50)

pdf(file="boxplot.pdf", width=8, height=6.5)
p+stat_compare_means(aes(group=Gene),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()

