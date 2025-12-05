# Data analysis

## Association test (logistic regression)

To perform an association analysis in PLINK, use the code below. `--bfile` calls all three .fam, .bim, and .bed files with the same prefix. `--pheno` tells PLINK to read the phenotype from a specific file, in this 
case, *pheno_sub.phe*. `--pheno-name` tells PLINK which column of phenotype values to pull from within the phenotype file. `--glm` is the flag that will run the association test. It stands for Generalized
Linear Model, and will perform a logistic regression when the phenotype data is in binary format (1/2/0). Following this flag with `allow-no-covars` means it can be run without adding a
covariate file. This command will produce a file labelled: "neonate_test.neonate_test.glm.logistic.hybrid". It contains the test results.

```
./plink2 --bfile filtered_data --pheno beluga_pheno_sub.phe --pheno-name neonate_test --glm allow-no-covars --allow-extra-chr --out neonate_test
```


When opening up the .glm.logistic.hybrid file, some of the important columns are as follows:
- 1 = Chromosome number
- 2 = Position
- 3 = ID
- 4 = Ref allele
- 5 = Alt allele
- 16 = p-value


After saving as a .txt file, rename columns 1, 2 and 3, to *CHR*, *BP*, and *SNP*, respectively. This is so that the important information is formatted correctly to be used in Rstudio with the qqman package. 


To make files containing significant SNPS, you can use the code below and specify the p-value threshold. The command `awk` filters/processes columns in text files. `'NR==1` tells it to keep the header upon creating
the new file. Here, `$16` specifies which column to process the values from, and the new file will be saved as sig.txt. This may become useful for running the WISH-R package and only including significant SNPs. In the example below, the p-value
threshold is p < 0.1. 

```
awk 'NR==1 || $16 < 0.1' neoplasia_test_genotypic.txt > sig.txt
```


## Epistatic interactions - WISH R package

To identify potential clusters of alleles based on epistatic interactions in R using the WISH-R package, we first need to produce a .tped file from the filtered bfiles. `--recode` creates a new fileset. `--recode12` converts alphabetical allele values to 1/2, where 
1 and 2 are the first and second alleles in the .bim, respectively. `--recode transpose` causes the columns and rows to interchange, since we need each individual to be a singular row for the WISH-R package. 

```
./plink --bfile filtered_data --recode --out WISH_data --allow-extra-chr
./plink --file WISH_data.ped --recode12 --allow-extra-chr 
./plink --file plink --recode -transpose --out transposed_12_data --allow-extra-chr
``` 

**Add in: would then need to follow R code (phenotype formatting and genotype formatting)


