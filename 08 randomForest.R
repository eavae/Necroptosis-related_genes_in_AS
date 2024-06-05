
library(randomForest)
library(limma)
library(ggpubr)
set.seed(123)

inputFile="normalize.txt"      
setwd("")

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
colnames(data)=gsub("-", "afaf", colnames(data))

afcon=as.matrix(table(sample[,1]))[1,1]
afcon=as.vector(afcon)
group=c(rep("con",afcon),rep("treat",nrow(data)-afcon))


rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()


optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)


importance=importance(x=rf2)
importance=as.data.frame(importance)
importance$size=rownames(importance)
importance=importance[,c(2,1)]
names(importance)=c("Gene","importance")


af=importance[1:20,]
p=ggdotchart(af, x = "Gene", y = "importance",
             color = "importance", 
             sorting = "descending",                     
             add = "segments",                             
             add.params = list(color = "lightgray", size = 2), 
             dot.size = 6,                        
             font.label = list(color = "white", size = 9,
                               vjust = 0.5),              
             ggtheme = theme_bw()         ,               
             rotate=TRUE                                       )
p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +
  grids()   

pdf(file="importance.pdf", width=6, height=6)
print(p1)
dev.off()

rfGenes=importance[order(importance[,"importance"], decreasing = TRUE),]
write.table(rfGenes, file="rfGenes.xls", sep="\t", quote=F, col.names=T, row.names=F)
