library(abSNF)
library(rdist)
library(Seurat)
library(ggplot2)

setwd("C:/Users/alexw/School/R-work/BREM")
data <- Read10X("./filtered_feature_bc_matrix") 
pbmc.rna <- CreateAssayObject(counts = data$`Gene Expression`)
pbmc.sal <- CreateAssayObject(counts = data$`Antibody Capture`)

pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.sal <- NormalizeData(pbmc.sal)

pbmc.rna <- FindVariableFeatures(pbmc.rna, selection.method = "vst", nfeatures = 5000)
pbmc.rna <- ScaleData(pbmc.rna)

var_genes <- VariableFeatures(pbmc.rna)
rmatrix <- as.matrix(GetAssayData(pbmc.rna)[var_genes,])
smatrix <- as.matrix(pbmc.sal@counts)

rmatrix <- as.matrix(pbmc.rna@scale.data)
smatrix <- as.matrix(pbmc.sal@data)

rmatrix <- t(rmatrix)
smatrix <- t(smatrix)

rna.dist <- rdist::pdist(rmatrix, metric = "euclidean", p = 2L)
sal.dist <- rdist::pdist(smatrix, metric = "euclidean", p = 2L)

rna.similarity2 <- affinityMatrix(rna.dist, K = 20, sigma=.5)
sal.similarity <- affinityMatrix(sal.dist, K = 20, sigma=.5)

save(sal.similarity, file = "salSimilarity.Rdata")

hist(rna.similarity[rna.similarity >= 0], col = 'lightblue', freq = TRUE)
hist(rna.similarity[rna.similarity > .02], col = 'lightblue', freq = TRUE)

mean(rna.similaity[, 2])
ggplot(as.data.frame(rna.similarity), row.names = 'Genes', col.names='Cells') +
aes(x=rna.similarity[, 1]) + 
  geom_histogram(binwidth=1, color="black", fill="red") + 
  labs(y = "Frequencies", title = "Top 1000 variable gene frequencies on a random barcode") +
  geom_vline(aes(xintercept=mean(rna.similarity[, 1])),
             color="blue", linetype="dashed", size=1)
library(Hmisc)
hist.data.frame(sal.similarity)
qplot(6000, data = as.data.frame(rmatrix), geom = 'histogram', binwidth=1,) + 
  labs(y = "Frequencies", title = "Top 1000 variable gene frequencies on a random barcode")