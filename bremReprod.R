plotBCells <- function(mat) {
  # top left is B cells
  x <- log(mat[1 , ])
  y <- log(mat[8, ])
  bcell.plot <- plot(x, y, pch=20, ylab = 'log CD19', xlab = ' log CD3', main = 'B cells -')
  abline(v=c(6.25,4.5), col=c("blue","blue"), lty = 2) # cd3
  abline(h=c(5.5,3), col=c("red", "red"), lty= 2) # cd19
  
  bcells <- colnames(mat)[which(x < 4.5 & y > 5.5)] # get the col names that pass the cond
  # which gives the indexes
  tcells <- colnames(mat)[which(x > 6.25 & y < 3)]
  tcells.mat <- mat[, which(x > 6.25 & y < 3)] # cd3+ and cd19- aka t cells
  bcell.mat <- mat[, which(log(mat[1, ]) < 4.5 & log(mat[8, ]) > 5.5)]
}


plotTCells <- function() {
  x <- log(tcells.mat[2, ])
  y <- log(tcells.mat[3, ])
  
  plot(x, y, pch=20, ylab = 'log CD8a', xlab = ' log CD4', main = 'T cells')
  abline(v=c(7,4.5), col=c("blue","blue"), lty = 2) # cd4
  abline(h=c(6.5,5), col=c("red", "red"), lty= 2) # cd8a
  
  # CD 4+ 
  cd4tcells <- colnames(tcells.mat)[which(x > 7 & y < 5)]
  # CD 8+   I guess the remaining 400 ish cells dont get classified. consistent
  cd8tcells <- colnames(tcells.mat)[which(x < 4.5 & y > 6.5)]
}

plotMonoCyteCells <- function() {
  notNeeded <- c(bcells,cd4tcells,cd8tcells)
  remainingCells <- mat[, -which(colnames(mat) %in% notNeeded)] # take subset that are not in notNeeded
  # use %in% instead bc different length vectors
  x <- log(remainingCells[4, ]) # cd14 
  y <- log(remainingCells[6, ]) # cd16
  monocyte.plot <- plot(x, y, pch=20, xlab = 'log CD14', ylab = 'log CD16', main = 'Monocytes -')
  abline(v=c(6,4), col=c("blue","blue"), lty = 2) 
  abline(h=c(6,5), col=c("red", "red"), lty= 2) 
  
  cd14cells <- colnames(remainingCells)[which(x > 6 & y < 5)]
  cd16.mat <- mat[, which(x < 4 & y > 6)]
}

plotNKCells <- function() {
  x <- log(cd16.mat[7, ]) # cd 56
  y <- log(cd16.mat[14, ]) # cd 127
  View(cd16.mat)
  
  nk.plot <- plot(x, y, pch=20, xlab = 'log cd56', ylab = 'log cd127', main = 'NK and CD16 Monocytes')
  abline(v=c(5,3.5), col=c("blue","blue"), lty = 2) 
  abline(h=c(3,4), col=c("red","red"), lty = 2) # idk what the cutoff should be
  
  nkcells <- colnames(cd16.mat)[which(x > 5 & y < 3.5)]
  cd16cells <- colnames(cd16.mat)[which(x < 3.5 & y > 4)]
}

plotDendtritic <- function() {
  DC1 <- pbmc.rna@counts[c('CD1C','FCER1A'), ]
  DC1.cells <- colnames(DC1)[which(DC1['CD1C', ] > 0 & DC1['FCER1A', ] > 0)] # take subset of col names
  
  cd14.neg <- colnames(mat[, which(log(mat[4, ]) < 6)]) # CD14 - cutoff is 6 
  DC1.cells <- intersect(cd14.neg, DC1.cells)
  
  DC2 <- as.matrix(pbmc.rna@counts[c('CD1C','HLA-DRB1'), ])
  DC2.cells <- colnames(DC2)[which(DC2['CD1C', ] > 0 & DC2["HLA-DRB1", ] > 0)]
  cd14.pos <- colnames(mat[, which(log(mat[4, ]) > 6)])  
  DC2.cells <- intersect(cd14.pos, DC2.cells)
  
  intersect(DC1.cells, DC2.cells) # should be zero
  
  pDC <- pbmc.rna@counts['IL3RA', ] # only grab the gene we care abt
  pDC <- pDC[which(pDC > 0)] # this for now dk the plus cutoff
  pDC.cells <- names(pDC)
  
  rm(DC1, DC2, pDC) # clean up a little
  rm(cd14.neg, cd14.pos)
  
}

if(!require(devtools)) install.packages("devtools")

library("devtools")
install_github("mojaveazure/seurat-disk")
install_github("tarot0410/BREMSC")
library(dplyr)
library(Seurat)
library(SeuratDisk) 
library(BREMSC)
library(factoextra)
library(NbClust)
setwd("C:/Users/alexw/School/R-work/BREM")

data <- Read10X("./filtered_feature_bc_matrix") 
pbmc.rna <- CreateAssayObject(counts = data$`Gene Expression`)
pbmc.sal <- CreateAssayObject(counts = data$`Antibody Capture`)

pbmcCombined <- CreateSeuratObject(counts = pbmc.rna)
pbmcCombined[['SAL']] <- pbmc.sal # add an assay

mat <- as.matrix(pbmc.sal@counts)

cluster <- function(SALmatrix) {
  result <- vector("list", length = 7865)
  for(j in 1:ncol(SALmatrix)) {       # for-loop over cells
    for(i in 1: nrow(SALmatrix)) {
      if(data[i,j] > log(400)) {
        result[[j]] <- 1  # B cell
        break;
      }
      data1[i, j] <- data1[i , j] + 10 
    }
  }
}

## Need to cluster the cells based on cell surface markers
# log(CD3+1) is < 4.XX AND log(CD19+1) is > 5.8X, it is B cell.

### NEED TO REPROD THE ARI AND AMI results. USE SAME PREFILTERING. 

plot1 <- FeatureScatter(pbmcCombined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#pbmcCombined <- subset(pbmcCombined, subset = nFeature_RNA > 300 & nFeature_RNA < 5000)

pbmc.rna <- NormalizeData(pbmcCombined@assays$RNA)
pbmc.sal <- NormalizeData(pbmcCombined@assays$SAL)

# remove 3 protein markers from the ADT
SALcounts <- GetAssayData(pbmc.sal)
SALcounts <- SALcounts[-(which(rownames(SALcounts) %in% c('CD8a-TotalSeqB','CD16-TotalSeqB','CD127-TotalSeqB'))),]
pbmc.sal <- subset(pbmc.sal, features = rownames(SALcounts))

### We identified seven cell types based on the biological knowledge of both protein and gene markers as the approximate
#truth, which is illustrated in Supplementary Figure S1. Examples of such cell type identification procedure are shown
#in Supplementary Figure S2. Taken together, >80% of single cells can be assigned to a specific cell type. Cells with
#uncertain cell types (not identified in the ground truth) were
#removed from computing ARIs and AMIs.
###

pbmc.rna <- FindVariableFeatures(pbmcCombined@assays$RNA, selection.method = "vst", nfeatures = 1000)
pbmc.rna <- ScaleData(pbmc.rna)
rnaReduc <- RunPCA(pbmc.rna, features = VariableFeatures(object = pbmc.rna))

var_genes <- VariableFeatures(pbmc.rna)
rmatrix <- t(as.matrix(
                GetAssayData(pbmc.rna)[var_genes,]
                ))

smatrix <- t(as.matrix(pbmc.sal@counts))

rmatrix <- t(rmatrix)
smatrix <- t(smatrix)

# max(rmatrix['JCHAIN', ])

# seven cell types
result3 <- BREMSC(smatrix, rmatrix, K=7, nChains=3, nMCMC=500)
occurences <- table(result$clusterID)

View(result$posteriorProb)

save(result3, file = "result3.Rdata")

# remove the 20% of cells that did not map to one of the 7 clusters.

plot(result2$vecLogLik, type = "l", xlab = "MCMC Iterations", ylab = "Log likelihood")

adjustedRandIndex(result2$clusterID, result3$clusterID)

.test = function() {
  data("dataADT")
  data("dataRNA")
  testRun_BREMSC = BREMSC(dataADT, dataRNA, K=4, nChains=3, nMCMC=100)
  plot(testRun_BREMSC$vecLogLik, type = "l", xlab = "MCMC Iterations", ylab = "Log likelihood")
}

