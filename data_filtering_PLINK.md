# Data filtering

## Quality control

Before running commands to filter/quality check, run the below lines before and after to count the number of variants and number of individuals/samples. This should tell you 
if alterations remove data or not (using `.bim` will give you the number of variants and `.fam` will give you the number of samples). 

```
wc -l 0.01_ff.bim 
wc -l 0.01_fully_filtered.fam
```


Next, run the code below. `--geno 0.05` filters out variants with high rates of missing genotypes greater than 5%. `--hwe 1e-6` filters out all variants which have Hardy-Weinberg equilibrium exact test p-values below the provided threshold. 
Lastly, `--maf 0.01` filters out all variants with a minor allele frequency below the provided threshold (1%). This can be completed in one command, along with `--make-bed` and `--out`, and the output will give you how many variants and samples were removed, 
in addition to what remains:

```
./plink -bfile 0.01_ff --geno 0.05 --hwe 1e-6 --maf 0.01 --make-bed --out qc_data --allow-extra-chr
```


Next, we can check for relatedness using the `king-cutoff` command and a threshold value. The threshold value will be important for what degree of relatedness you want to remove. For our case, we will use first
degree relatives, equivalent to a threshold value of **0.177**. This command does not automatically remove samples, and therefore requires two commands to first identify samples that should be removed (which are
saved in a `.king.cutoff.out.id` file), and a subsequent command to actually remove them and create a new set of bfiles. The flag `king-cutoff` requires plink2.

```
./plink2 --bfile qc_data --make-king-table --king-cutoff 0.177 --out king_cutoff
./plink2 --bfile qc_data --remove king_cutoff.king.cutoff.out.id --make-bed --out filtered_data
```


After quality checks, the data is now stored in a set of bfiles with the prefix *filtered_data*.


## Aligning phenotype mortality data
After filtering has been completed and sites/samples have been removed, the code below can be run in Rstudio to align/reduce the phenotype data file to match genetic data. Phenotype data will now be stored as *pheno_data.phe*. 

```
# To make subet of filtered .fam file to match the .phe file
fam <- read.table("filtered_data.fam")
phe <- read.table("pheno_data.phe", header=TRUE)

## Keep only rows in .fam
phe_sub <- phe[phe$FID %in% fam$V1 & phe$IID %in% fam$V2, ]

write.table(phe_sub, "pheno_sub.phe", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
```


*Note*: you can use this code to count how many 0's, 1's and 2's for each phenotype column, where the number after `$` is the column number. 
```
awk 'NR>1{print $6}' beluga_pheno.phe | sort | uniq -c
```



