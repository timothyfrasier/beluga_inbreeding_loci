###############################
# COMMANDS FOR WISH-R ANALYSIS #
# FOR RW DATA                 #
###############################
library(data.table)

#---------------------------#
# Install required packages #
#---------------------------#

install.packages("devtools", dependencies = TRUE)
install.packages("curl", dependencies = TRUE)
install.packages("httr", dependencies = TRUE)
install.packages("doParallel", dependencies = TRUE)
install.packages("foreach", dependencies = TRUE)
install.packages("fastcluster", dependencies = TRUE)
install.packages("Rcpp", dependencies = TRUE)
install.packages("RcppEigen", dependencies = TRUE)
install.packages("data.table", dependencies = TRUE)
install.packages("corrplot", dependencies = TRUE)
install.packages("heatmap3", dependencies = TRUE)
install.packages("flashClust", dependencies = TRUE)
install.packages("bigmemory", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)


#-------------------------------#
# Install BioConductor Packages #         
#-------------------------------#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute", force=TRUE) # This one takes a long time, changed, output asked me to add force=TRUE
BiocManager::install("preprocessCore", force=TRUE)
BiocManager::install("GO.db", force=TRUE)
BiocManager::install("AnnotationDbi", force=TRUE) # changed, output asked me to add force=TRUE

install.packages("WGCNA", dependencies = TRUE)

library("devtools")
install_github("QSG-Group/WISH", force=TRUE)
######################

### Installing devtools
install.packages(c("devtools","curl", "httr"))

### Install the rest of the dependencies
install.packages(c("doParallel", "foreach","fastcluster", "Rcpp", "RcppEigen", "data.table", "corrplot", "heatmap3", "flashClust", "bigmemory", "parallel", "ggplot2"))

### Install WISH
source("https://install-github.me/QSG-Group/WISH")
# or
library("devtools")
install_github("QSG-Group/WISH")

# Install required packages
install.packages("BiocManager")
BiocManager::install(c(
  "AnnotationDbi",
  "GO.db",
  "preprocessCore",
  "impute"
), force=TRUE)
install.packages("WGCNA")
install.packages("devtools")
library(devtools)
library(WISH)
install.packages("remotes")
remotes::install_github("cran/Epistasis", force=TRUE)
install.packages("WGCNA")
library(WGCNA)

install.packages("epistasis")
library(epistasis) ## 

devtools::install_github("QSG-Group/WISH", force=TRUE)
library(data.table)


#------------------#
# Read in Data     #
#------------------#

#--- Infection Data ---#
infectiona_ped  <- fread("filtered12_0.1_infecta.ped", data.table = F)
infectiona_tped <- fread("tfiltered12_0.1_infecta.tped", data.table = F)
infectiong_ped  <- fread("filtered12_0.1_infectg.ped", data.table = F)
infectiong_tped <- fread("tfiltered12_0.1_infectg.tped", data.table = F)
phenotypes <- read.table("pheno_sub.txt", header = TRUE, sep = "\t")
infection_phenos <- phenotypes[, 11]

#--- Neonate Data ---#
nta_ped  <- fread("filtered12_0.1_nta.ped", data.table = F)
nta_tped <- fread("tfiltered12_0.1_nta.tped", data.table = F)
ntg_ped  <- fread("filtered12_0.1_ntg.ped", data.table = F)
ntg_tped <- fread("tfiltered12_0.1_ntg.tped", data.table = F)
nta_phenos <- phenotypes[, 13]

#--- Neoplasia Data ---#
npa_ped  <- fread("filtered12_0.1_npa.ped", data.table = F)
npa_tped <- fread("tfiltered12_0.1_npa.tped", data.table = F)
npg_ped  <- fread("filtered12_0.1_npg.ped", data.table = F)
npg_tped <- fread("tfiltered12_0.1_npg.tped", data.table = F)
npa_phenos <- phenotypes[, 12]

#--- Dystocia Data ---#
da_ped  <- fread("filtered12_0.1_da.ped", data.table = F)
da_tped <- fread("tfiltered12_0.1_da.tped", data.table = F)
dg_ped  <- fread("filtered12_0.1_dg.ped", data.table = F)
dg_tped <- fread("tfiltered12_0.1_dg.tped", data.table = F)
da_phenos <- phenotypes[, 10]

#---------------------#
# Run the analyses    #
#---------------------#
library(WISH)
infecta_genotypes <- generate.genotype(infectiona_ped, infectiona_tped)
infectg_genotypes <- generate.genotype(infectiong_ped, infectiong_tped)
ntg_genotypes   <- generate.genotype(ntg_ped, ntg_tped)
nta_genotypes   <- generate.genotype(nta_ped, nta_tped)
npg_genotypes   <- generate.genotype(npg_ped, npg_tped)
npa_genotypes   <- generate.genotype(npa_ped, npa_tped)
dg_genotypes   <- generate.genotype(dg_ped, dg_tped)
da_genotypes   <- generate.genotype(da_ped, da_tped)

infecta_correlations <- epistatic.correlation(infection_phenos, infecta_genotypes, threads = 4, test=F, simple=F)
infectg_correlations <- epistatic.correlation(infection_phenos, infectg_genotypes, threads = 4, test=F, simple=F)
ntg_correlations    <- epistatic.correlation(nta_phenos, ntg_genotypes, threads = 4, test=F, simple=F)
nta_correlations    <- epistatic.correlation(nta_phenos, nta_genotypes, threads = 4, test=F, simple=F)
npg_correlations    <- epistatic.correlation(npa_phenos, npg_genotypes, threads = 4, test=F, simple=F)
npa_correlations    <- epistatic.correlation(npa_phenos, npa_genotypes, threads = 4, test=F, simple=F)
dg_correlations    <- epistatic.correlation(da_phenos, dg_genotypes, threads = 4, test=F, simple=F)
da_correlations    <- epistatic.correlation(da_phenos, da_genotypes, threads = 4, test=F, simple=F)

library(WGCNA)
library(flashClust) 

#--- Infection allelic ---#
# Based on Correlations
n.snps = dim(infecta_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- infecta_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- infecta_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = sft$powerEstimate, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
infecta_modules_cor <- output   # Select appropriate modules


#--- Infection genotypic ---#
# Based on Correlations
n.snps = dim(infectg_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- infectg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- infectg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = 6, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules", ylim=c(0.7,0.9))
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
infectg_modules_cor <- output   # Select appropriate modules




#--- nta ---# 
# Based on Correlations
n.snps = dim(nta_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- nta_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- nta_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = 6, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30, cutHeight=0.99)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
nta_modules_cor <- output   # Select appropriate modules


#--- ntg ---# 
# Based on Correlations
n.snps = dim(ntg_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- ntg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- ntg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = 6, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
ntg_modules_cor <- output   # Select appropriate modules

#--- npa ---#
# Based on Correlations
n.snps = dim(npa_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- npa_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- npa_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = sft$powerEstimate, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 50)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
npa_modules_cor <- output   # Select appropriate modules


#--- npg ---#
# Based on Correlations
n.snps = dim(npg_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- npg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- npg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = sft$powerEstimate, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
npg_modules_cor <- output   # Select appropriate modules

MElist <- moduleEigengenes(npg_genotypes, colors = moduleColors) 
MEs <- MElist$eigengenes 
MEs

#--- da ---# 
# Based on Correlations
n.snps = dim(da_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- da_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- da_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = 6, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 50)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
da_modules_cor <- output   # Select appropriate modules


#--- dg ---# 
# Based on Correlations
n.snps = dim(dg_correlations$Coefficients)[1]  # Specify appropriate correlations
temp_corr <- dg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr < 0] <- 0
temp_corr <- temp_corr/(max(temp_corr))
corr <- temp_corr
temp_corr <- dg_correlations$Coefficients  # Specify appropriate correlations
temp_corr[temp_corr > 0] <- 0
temp_corr <- temp_corr/(abs(min(temp_corr)))
corr[temp_corr < 0] <- temp_corr[temp_corr < 0]

sft <-  pickSoftThreshold(corr, powerVector = c(seq(1, 10, 0.1), c(12:22)), verbose = 5)
connectivity <- adjacency.fromSimilarity(corr, power = 6, type = "signed")
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(connectivity, xlab = "connectivity")
connectivity[is.nan(connectivity)] <- 0 # Replacing NaN with 0
scaleFreePlot(connectivity)
par(mfrow = c(1, 1))
select.snps <- corr[rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps, rank(-colSums(connectivity), ties.method = "first") <= 
                      n.snps]
select.snps[, c(1:ncol(select.snps))] <- sapply(select.snps[, c(1:ncol(select.snps))], as.numeric)
adjMat <- adjacency.fromSimilarity(select.snps, power = sft$powerEstimate, type = "signed")
dissTOM <- 1 - (TOMsimilarity(adjMat))
genetree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <-cutreeDynamic(dendro = genetree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
moduleColors = labels2colors(dynamicMods)
plotDendroAndColors(genetree, moduleColors, dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
output <- list(select.snps, connectivity, adjMat, dissTOM, 
               genetree, dynamicMods, moduleColors, sft$powerEstimate)
names(output) <- c("SNPs", "connectivity", "adjMat", "dissTom", 
                   "genetree", "modules", "modulecolors", "power.estimate")
dg_modules_cor <- output   # Select appropriate modules

