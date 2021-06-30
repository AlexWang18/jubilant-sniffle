library(abSNF)
library(rdist)
library(Seurat)
library(ggplot2)

vn_entropy <- function(mat) {
  EV <- eigen(mat)
  EV <- EV$values
  EV <- EV[EV>0]
  log2.EV <- log2(EV)
  result <- -dot(EV, log2.EV)
}


setwd("./research/data")
data <- Read10X("./filtered_feature_bc_matrix") 
pbmc.rna <- CreateAssayObject(counts = data$`Gene Expression`)
pbmc.sal <- CreateAssayObject(counts = data$`Antibody Capture`)

pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.sal <- NormalizeData(pbmc.sal)

# impute missing genes, etc
pbmcCombined <- CreateSeuratObject(counts = pbmc.rna)
pbmcCombined[['SAL']] <- pbmc.sal # add an assay

pbmcCombined = pbmcCombined[rowSums(pbmcCombined) != 0, ] # remove genes never expressed by any cell
pbmc.rna <- FindVariableFeatures(pbmc.rna, selection.method = "vst", nfeatures = 5000)
pbmc.rna <- ScaleData(pbmc.rna)

var_genes <- VariableFeatures(pbmc.rna)
rmatrix <- as.matrix(GetAssayData(pbmc.rna)[var_genes,])
smatrix <- as.matrix(pbmc.sal@counts)

rmatrix <- as.matrix(pbmc.rna@scale.data)
smatrix <- as.matrix(pbmc.sal@data)

rmatrix <- t(rmatrix)
smatrix <- t(smatrix)

#load("~/research/data/rnaDist.Rdata")
# dist2 is with data
rna.dist2 <- rdist::pdist(rmatrix, metric = "euclidean", p = 2L)
save(rna.dist2, file="rnaDist3.Rdata")
sal.dist2 <- rdist::pdist(smatrix, metric = "euclidean", p = 2L)

rna.similarity2 <- affinityMatrix(rna.dist2, K = 20, sigma=.5)
sal.similarity2 <- affinityMatrix(sal.dist, K = 20, sigma=.5) # 2 is with normalized.

qplot(rna.similarity[rna.similarity >= 0], data = as.data.frame(rna.similarity), geom = 'histogram', binwidth=1,) + 
  labs(y = "Frequencies", title = "Top 5000 variable gene frequencies")

r.h <- hist(rna.similarity[rna.similarity <= .005], col = 'lightblue', freq = TRUE, breaks=10, main="RNA Similarity Matrix Histogram")
text(r.h$mids,r.h$counts,labels=r.h$counts, adj=c(0.6, -0.6), cex = .5)

r.h2 <- hist(-log10(rna.similarity[]), col = 'lightblue', freq = TRUE, breaks=10, main="RNA Similarity Matrix Histogram")
text(r.h2$mids,r.h2$counts, labels=r.h2$counts, adj=c(0.6, -0.6), cex = .5)

#  Freedman–Diaconis rule of bin width
r.hFD <- hist(-log10(rna.similarity[]), col = 'lightblue', freq = TRUE, breaks="FD", main="RNA Similarity Matrix Histogram")
text(r.hFD$mids,r.hFD$counts,labels=r.hFD$counts, adj=c(0.6, -0.6), cex = .5)

s.h <- hist(-log10(sal.similarity[]), col = 'lightblue', freq = TRUE, breaks=10, main="SAL Similarity Matrix Histogram")
text(s.h$mids,s.h$counts,labels=s.h$counts, adj=c(0.6, -0.6), cex = .5)
entropy(r.h2$counts, method='ML')
# [1] 1.29663
entropy(s.h$counts, method = 'ML') # estimates maximum liklihood / Shannon entropy
# [1] 1.614044
x1 = runif(10000) # normal distrbution
y1 = discretize(x1, numBins=10, r=c(0,1))
entropy(y1) # highest uncertainity - 2.302

#hist(sal.similarity2[sal.similarity2 < .001])
#s.h <- hist(sal.similarity2[sal.similarity2 <.00001], col = 'lightblue', freq = TRUE, breaks=10, main="SAL Similarity Matrix Histogram")
#text(s.h$mids,s.h$counts,labels=s.h$counts, adj=c(0.6, -0.6), cex = .5)
#length(which(sal.similarity[]<0.000005))
