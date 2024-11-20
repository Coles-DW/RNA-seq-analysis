setwd("C:/Users/30059564/OneDrive - Western Sydney University/Myrtle rust/Apsidii control paper/Plant RNA-seq/filtering counts/deseq2/plsda")
#pls-da or pca analysis in R
#required the mixOmics package in R
BiocManager::install("mixOmics")
library(mixOmics)
#data entry, need a data input sheet with samples in rows and data in columns, code following
#is if you are using an rld (rlogtransformed) transcriptomic data set from deseq2
rlog1 <- as.data.frame(assay(rld))
rlogdata <- as.data.frame(t(rlog1))
myData <- rlogdata[-c(1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 19, 20, 21), ]
#need a condition vector with the classification of each sample
conditionvector <- as.data.frame(condition)
write.csv(conditionvector, "definegroupsrvss.csv")
conditionvector <- read.csv("definegroupsrvss.csv")
#conditionvector for all samples
conditionvector2 <- as.data.frame(condition)
#ADD grouping column if doing plsda
#generate pca
mypca1<-pca(rlogdata)
mypca2<-pca(myData)
#or for plsda
myplsda<-plsda(myData,conditionvector$condition)
#plot pca and add in confidence ellipses and legend
plotIndiv(mypca1, ellipse=TRUE, ind.names=FALSE, legend=TRUE, group=conditionvector2$condition, pch=as.factor(conditionvector2$condition), legend.title="Treatment", title="PCA", cex=0.5)
plotIndiv(mypca2, col = color.mixo(4:5), ellipse=TRUE, ind.names=FALSE, legend=TRUE, group=conditionvector$condition, pch=as.factor(conditionvector$condition), legend.title="Treatment", title="PCA", cex=0.5)
#plot plsda and add in confidence ellipses and legend
plotIndiv(myplsda, col = color.mixo(4:5), ellipse=TRUE, ind.names=FALSE, legend=TRUE, group=conditionvector$condition, pch=as.factor(conditionvector$condition), legend.title="Treatment", legend.title.pch="control grouping", title="PLS-DA", cex=0.5)
#extracting the loadings
loadings_rvss<-as.data.frame(mypca[["loadings"]][["X"]])
write.table(loadings_rvss, file="pcaloadings_rvss.txt", sep="\t")
#extracting the coordinates
coordinates_rvss<-as.data.frame(mypca[["variates"]][["X"]])
write.table(coordinates_rvss, file="pcacoordinates_rvss.txt", sep="\t")

myData2 <- t(myData)
write.csv(myData2, "rldcountsrvss.csv")
