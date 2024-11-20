####DESEQ2 OZMYC data############
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("reshape")
BiocManager::install("ggplot2")
##################################

library("DESeq2")
library("ggplot2")
library('reshape')
workingDir = "C:/Users/30059564/OneDrive - Western Sydney University/Myrtle rust/Apsidii control paper/Plant RNA-seq/filtering counts/deseq2"
setwd(workingDir)
data=read.table("apsidii counts_filtered.txt", header = TRUE, row.names = 1)
colnames(data)
sample=colnames(data)
#condition = c("ContS","ContS","ContS","InocS","InocS","InocS","HypS","HypS","HypS",
              #"InocHypS","InocHypS","InocHypS","ContR","ContR","ContR","InocR","InocR","InocR")

condition = c("ContS","ContS","ContS","HypS","HypS","HypS","ContR","ContR","ContR")

condition = c("InocS","InocS","InocS","InocHypS","InocHypS","InocHypS","InocR","InocR","InocR")
condition = c("Myrtle_Rust_Susceptible","Myrtle_Rust_Susceptible","Myrtle_Rust_Susceptible","Myrtle_Rust_Hyperparasite_Susceptible","Myrtle_Rust_Hyperparasite_Susceptible","Myrtle_Rust_Hyperparasite_Susceptible","Myrtle_Rust_Resistant","Myrtle_Rust_Resistant","Myrtle_Rust_Resistant")
colData=cbind(colnames(data), sample, condition)
######verification des col -treatment###
colData=data.frame(colData)
write.table(colData, file="MRcol_Samples_Conditions_myrtle.txt", sep="\t")
######boxplot of data transnformation
epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(data + epsilon)), breaks=100, col="blue", border="white",
     main="Expressed_Log2-transformed counts per gene_myrtle", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
dev.copy(png,'Log2-transformed_Histogram_myrtle.png')
dev.off()
condition_1 <- as.data.frame(t(condition))
boxplot(log2(data + epsilon), condition_1$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(counts +1)")
dev.copy(png,'boxplot_expressed_transformeddata_myrtle.png')
dev.off()
######dds=data deseq
dds= DESeqDataSetFromMatrix(data, colData, ~condition)
dds= DESeq(dds)
results(dds)
str(dds)
######sortie donnees normalises############
norm_counts=counts(dds, normalized=TRUE)
write.table(norm_counts, file="myrtle_normalized_expressed.txt", sep="\t")
##rlog tranformation
rld<- rlogTransformation(dds, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)

install.packages("gplots")
install.packages("ggplot2")
install.packages("RColorBrewer")

library("gplots")
library('RColorBrewer')
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(10, 10))
dev.copy(png,'Expressed_myrtle_sampletosample_heatmap.png')
dev.off()
pca_1<-print(plotPCA(rld))
pcaData <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

install.packages("ggforce")
library(ggforce)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=1) + 
  ggforce::geom_mark_ellipse(aes(fill = condition,
                                 color = condition)) +
  ylim(-10, 10) +
  xlim(-20, 20) +
  xlab(paste0("\n PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance \n")) + 
  theme_classic()
dev.copy(png, 'priming_PCA.png')
dev.off()
#####savoir la dimension du fichier, nombre col et lignes)###
dim(norm_counts)
str(dds)
#####Filtering 10 reads per sample repetition->excel########
##############Hclust##########
## transfo en log2(x+1) ##
norm_counts_log=log2(norm_counts+1)
## calcul correlation de Pearson##
r = cor(norm_counts_log, method = "pearson")
write.table(r, file="myrtle_correlation_values.txt", sep="\t")
## calcul de la dissimilarité ##
d=1-r
## Clustering ##
h = hclust(as.dist(d), method = "complete", members = NULL)
plot(h)
dev.copy(pdf, 'Expressed_myrtle_hclust.pdf')
dev.off()

####comparaison####
####Pairwise comparison#####
res1 <- results(dds, contrast=c("condition","ContR","ContS"))
res2 <- results(dds, contrast=c("condition","HypS","ContS"))

res1 <- results(dds, contrast=c("condition","InocR","InocS"))
res2 <- results(dds, contrast=c("condition","InocHypS","InocS"))

res1 <- results(dds, contrast=c("condition","Myrtle_Rust_Resistant","Myrtle_Rust_Susceptible"))
res2 <- results(dds, contrast=c("condition","Myrtle_Rust_Hyperparasite_Susceptible","Myrtle_Rust_Susceptible"))
######Mise en forme#########
res1=res1[,c(2,6)]
res2=res2[,c(2,6)]


#Fichier Final
colnames(res1)=c("Log2Foldchange_InocRvsInocS","padj_InocRvsInocS")
colnames(res2)=c("Log2Foldchange_InocHypSvsInocS","padj_InocHypSvsInocS")

colnames(res1)=c("Log2Foldchange_Myrtle_Rust_ResistantvsMyrtle_Rust_Susceptible","padj_Myrtle_Rust_ResistantvsMyrtle_Rust_Susceptible")
colnames(res2)=c("Log2Foldchange_Myrtle_Rust_Hyperparasite_SusceptiblevsMyrtle_Rust_Susceptible","padj_Myrtle_Rust_Hyperparasite_SusceptiblevsMyrtle_Rust_Susceptible")

myrtlegeneexpressionresults=cbind(res1, res2)
write.table(myrtlegeneexpressionresults, file="Differentialexpression_myrtle pathogen.txt", sep="\t")

