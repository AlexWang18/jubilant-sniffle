nlog <- function(num) {
  return(-log10(num));  
}

sumColumns <- function(row) {
  Reduce('+', lapply(row, nlog))
}

getWeights <- function(matrix) {
  # nlog(p_k)/ sum of that feature
  apply(matrix, 1, sumColumns)
}

euclidean <- function(a, b) {
  sqrt(sum((a - b)^2)) # * weights
}

library("devtools")
library(dplyr)
library(rdist)
library(Seurat)
library(SeuratDisk) 
library(abSNF)
setwd("C:/Users/alexw/School/R-work/BREM")

data <- Read10X("./filtered_feature_bc_matrix") 
pbmc.rna <- CreateAssayObject(counts = data$`Gene Expression`)
pbmc.sal <- CreateAssayObject(counts = data$`Antibody Capture`)

pbmcCombined <- CreateSeuratObject(counts = pbmc.rna)
pbmcCombined[['SAL']] <- pbmc.sal # add an assay

plot1 <- FeatureScatter(pbmcCombined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmcCombined <- subset(pbmcCombined, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 18000)

pbmc.rna <- NormalizeData(pbmcCombined@assays$RNA)
pbmc.sal <- NormalizeData(pbmcCombined@assays$SAL)

# remove 3 protein markers from the ADT
SALcounts <- GetAssayData(pbmc.sal)
SALcounts <- SALcounts[-(which(rownames(SALcounts) %in% c('CD8a-TotalSeqB','CD16-TotalSeqB','CD127-TotalSeqB'))),]
pbmc.sal <- subset(pbmc.sal, features = rownames(SALcounts))

pbmc.rna <- FindVariableFeatures(pbmcCombined@assays$RNA, selection.method = "vst", nfeatures = 500)
pbmc.rna <- ScaleData(pbmc.rna)

rnaWeight <- getWeights()

rna.matrix <- t(pbmc.rna@scale.data)
sal.matrix <- t(as.matrix(pbmc.sal@data))

rnaDist <- pdist(rna.matrix, metric = "euclidean", p = 2)
salDist <- pdist(sal.matrix, metric = "euclidean", p = 2) 

# Next, construct similarity graphs
graphRNA = affinityMatrix(rnaDist)
graphSAL= affinityMatrix(salDist)

list <- list(graphRNA, graphSAL) # each element of the list should be a square

fused = SNF(list, 20, 20) # 20 neighbor

