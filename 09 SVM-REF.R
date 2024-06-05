


setwd("")
library(tidyverse)
library(glmnet)
source('msvmRFE.R')  
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
library(limma)



inputFile="normalize.txt"      

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
afcon=as.matrix(table(sample[,1]))[1,1]
afcon=as.vector(afcon)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(as.numeric(group))
rownames(group)=rownames(data)
colnames(group)="Type"
input <- as.data.frame(cbind(group,data))
input$Type=as.factor(input$Type)

svmRFE(input, k = 10, halve.above = 100) 
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) 
top.features = WriteFeatures(results, input, save=F) 
head(top.features)

write.csv(top.features,"feature_svm.csv")


featsweep = lapply(1:18, FeatSweep.wrap, results, input) 


no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
pdf("svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()

pdf("svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) 
dev.off()


which.min(errors) 
