# WISH-R: epistasis analysis

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
remotes::install_github("cran/Epistasis")
install.packages("WGCNA")
library(WGCNA)

install.packages("epistasis")
library(Epistasis) ## problem with epistasis package

devtools::install_github("QSG-Group/WISH")
library(data.table)

# Neonate mortality - allelic (original effort/method)
## Install datasets
tped_nta <- fread("tfiltered12_0.1_nta.tped", data.table=F)
ped_nta <- fread("filtered12_0.1_nta.ped", data.table=F)

tped_nta_test <- fread("tfiltered_0.1_nta.tped", data.table=F)

pheno_sub <- read.table("pheno_sub.txt")
phenosub <- read.table("phenosub_wish_nta.txt", skip = 1, header = FALSE, colClasses = c("NULL", "numeric"))
pheno <- pheno_sub[-1, -c(1:6)] # Column V13 for neonate mortality

# Generate genotype file (using already reduced based on p-value results from gwas)
genotype <- generate.genotype(ped_nta, tped_nta, gwas.p=NULL, snp.id=NULL)
genotype_test <- generate.genotype(ped_nta, tped_nta_test, gwas.p=NULL, snp.id=NULL)


# tried changing phenotype file
fam <- read.table("filtered_data.fam")
csv <- read.csv("pheno_wish_nta.csv", header=TRUE)

## Keep only rows in .fam
pheno_wish_nta <- csv[csv$IID %in% fam$V2, ]

write.table(pheno_wish_nta, "phenosub_wish_nta.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

## Convert allele pairs to genotype dosage (only need to do if don't use the generate.genotype command)
#geno_raw_nta <- tped_nta[, -c(1:4)] # remove the first four columns
#geno_numeric_nta <- apply(geno_raw_nta, 1, function(row) {
#alleles <- matrix(row, ncol = 2, byrow = TRUE) # convert allele pairs to 0/1/2 where 0=A1A1, 1=A1A2, 2=A2A2
 # rowSums(alleles == unique(alleles)[2])
#})

## Transpose so rows are individuals and columns are SNPs (only need to do if don't use the generate.genotype command)
#geno_numeric_nta <- t(geno_numeric_nta)

## Pruning for linkage disequilibrium (only need to do this if pruning isn't done in PLINK)
#LD_genotype_nta <- LD_blocks(tped_nta)
#tped_nta <- LD_genotype_mta$genotype

## Helpful to make sure dimensions are the same for both files.
dim(genotype_test)
dim(phenosub)

## Get estimates for the model run time
epistatic.correlation(phenosub, genotype_test, threads = 4, test = TRUE)


## Run the model
correlations_nta <- epistatic.correlation(phenosub, genotype_test, threads=4, test=F, glm=T)
#correlations_nta <- epistatic.correlation(pheno$V13, genotype, threads=4, test=F, simple=T)

## Correcting
correlations_nta$Coefficients[(is.na(correlations_nta$Coefficients))]<-0
correlations_nta$Pvalues[(is.na(correlations_nta$Pvalues))]<-1
#cor_mat <- as.matrix(data.frame(correlations)) # tried to view output as a matrix

#correlations_dub <- correlations_nta
#head(correlations_dub)

## Generate SNP gene modules
modules <- generate.modules(correlations_nta) # having issues with this line
ME <- moduleEigengenes(genotype, colors=modules$modulecolors, softPower=modules$power.estimate)

## Calculate association to cause of death
cor(pheno$V13, ME$eigengenes)

## Visualize ** difficulty
genome.interaction(tped_nta, correlations_nta)
pseudo_manhattan(tped_nta, correlations_nta, values="c")
#pseudo_manhattan(tped_nta, correlations_nta, values="c")


##### trying .traw method
## Install datsets
traw_nta <- fread("tfiltered12_0.1_nta.traw", data.table=F)
tped_nta <- fread("tfiltered12_0.1_nta.tped", data.table=F)
ped_nta <- fread("filtered12_0.1_nta.ped", data.table=F)
pheno_sub <- read.table("pheno_sub.txt")
pheno <- pheno_sub[-1, -c(1:6)] # Column V13 for neonate mortality
#pheno_matrix <- as.matrix(pheno_sub)

# Remove non-genotype columns
geno_only_nta <- traw_nta[, 7:ncol(traw_nta)]

## Transpose so rows are individuals and columns are SNPs
genotype_nta <- t(geno_only_nta)

## Helpful to make sure dimensions are the same for both files.
dim(geno_only_nta)
dim(pheno)
dim(genotype_nta)

## Get estimates for the model run time
#epistatic.correlation(pheno$V13, genotype_nta, threads = 4, test = T, glm = T)
epistatic.correlation(pheno$V13, genotype, threads = 4, test = T)

## Run the model
pheno_test <- pheno$V13
as.numeric(pheno_test)
correlations_nta <- epistatic.correlation(pheno_test, genotype, threads=4, test=F)

## Correcting
correlations_nta$Coefficients[(is.na(correlations_nta$Coefficients))]<-0
correlations_nta$Pvalues[(is.na(correlations_nta$Pvalues))]<-1
#cor_mat_nta <- as.matrix(data.frame(correlations_nta))
correlations_nta$Coefficients[correlations_nta$Coefficients > 1000 | correlations_nta$Coefficients < -1000] <- 0 #code that Ricky had


## Generate SNP gene modules
modules_nta <- generate.modules(correlations_nta)
ME <- moduleEigengenes(genotype_nta, colors=modules$modulecolors, softPower=modules$power.estimate)

## Calculate association to cause of death
cor(pheno$V13, ME$eigengenes)


## Visualize ** difficulty
genome.interaction(tped_nta, correlations_nta)
pairwise.chr.map(15, 20, tped_nta, correlations_nta, 25) # from Ricky code
pseudo_manhattan(tped_nta, correlations_nta, values="p")
pseudo_manhattan(tped_nta, correlations_nta, values="c")


