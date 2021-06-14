length(bcells) + length(cd14cells) + length(cd16cells) + length(cd4tcells) + length(cd8tcells) + length(nkcells) + dendtritic 

stuff <- c(bcells, cd14cells, cd16cells)
stuff2 <- c(cd4tcells, cd8tcells, nkcells)

intersect(cd14cells, cd16cells)
intersect(stuff, stuff2)

clusteredCells <- c(bcells, cd14cells, cd16cells, cd4tcells, cd8tcells, nkcells, DC1.cells, DC2.cells, pDC.cells)
x <- list(bcells, cd14cells, cd16cells, cd4tcells, cd8tcells, nkcells, DC1.cells, DC2.cells, pDC.cells)
length(unique(clusteredCells))

myls <- vector("list", length = 7865) # our approx truth clustering.
# need to reduce the size. should be 6554

intersect(clusteredCells, colnames(pbmc.rna@counts))

# For the clusters: bcells are 1, # cd4t are 2, cd8t are 3 
# cd14 are 4, cd16 are 5,nk are 6, dendritic are 7
un <- unlist(x)
res <- Map(`[`, x, relist(!duplicated(un), skeleton = x))
identical(res, res.list)
#[1] TRUE
removeDups <- function (x) { # does not work
  for(k in 1:(length(x) - 1)) {
    for (l in (k + 1):length(x)) {
      to.remove <- which(x[[k]] %in% x[[l]])
      x[[k]] <- x[[k]][-to.remove]
    }
  }
  return(x)
}

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
  
  for (c in bcells) { # take only clustered
    myls[[count]] <- 1
    count <- count+1 
  }
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
  for (c in cd4tcells) {
    myls[[count]] <- 2
    count <- count+1 
  }
  for (c in cd8tcells) {
    myls[[count]] <- 3
    count <- count+1 
  }
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
  for (c in cd14cells) {
    myls[[count]] <- 4
    count <- count+1 
  }
  
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
  
  for (c in cd16cells) {
    myls[[count]] <- 5
    count <- count+1 
  }
  for (c in nkcells) {
    myls[[count]] <- 6
    count <- count+1 
  }
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
  
  dendtritic <- length(DC1.cells) + length(DC2.cells) + length(pDC.cells)
  for (c in 1:dendtritic) {
    myls[[count]] <- 7
    count <- count+1 
  }
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

# plot1 <- FeatureScatter(pbmcCombined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sub <- matrix(unique(clusteredCells), ncol=1) # take only classified cells
pbmcCombined <- pbmcCombined[, sub[, 1]]

pbmc.rna <- NormalizeData(pbmcCombined@assays$RNA)
pbmc.sal <- NormalizeData(pbmcCombined@assays$SAL)

# remove 3 protein markers from the ADT
SALcounts <- GetAssayData(pbmc.sal)
# select rows besides those listed
SALcounts <- SALcounts[-(which(rownames(SALcounts) %in% c('CD8a-TotalSeqB','CD16-TotalSeqB','CD127-TotalSeqB'))),] 
pbmc.sal <- subset(pbmc.sal, features = rownames(SALcounts))
rm(SALcounts)

pbmc.rna <- FindVariableFeatures(pbmcCombined@assays$RNA, selection.method = "vst", nfeatures = 1000)
pbmc.rna <- ScaleData(pbmc.rna)
rnaReduc <- RunPCA(pbmc.rna, features = VariableFeatures(object = pbmc.rna))

var_genes <- VariableFeatures(pbmc.rna)
rmatrix <- t(as.matrix(
                GetAssayData(pbmc.rna)[var_genes,]
                ))
rm(var_genes)

smatrix <- t(as.matrix(pbmc.sal@counts))
rmatrix <- t(rmatrix)
smatrix <- t(smatrix)

# max(rmatrix['JCHAIN', ])

# seven cell types
result7 <- BREMSC(smatrix, rmatrix, K=7, nChains=3, nMCMC=500)
occurences <- table(result$clusterID)


save(result7, file = "result7.Rdata")

# remove the 20% of cells that did not map to one of the 7 clusters.

plot(result2$vecLogLik, type = "l", xlab = "MCMC Iterations", ylab = "Log likelihood")

adjustedRandIndex(result2$clusterID, result3$clusterID)

.test = function() {
  data("dataADT")
  data("dataRNA")
  testRun_BREMSC = BREMSC(dataADT, dataRNA, K=4, nChains=3, nMCMC=100)
  plot(testRun_BREMSC$vecLogLik, type = "l", xlab = "MCMC Iterations", ylab = "Log likelihood")
}

