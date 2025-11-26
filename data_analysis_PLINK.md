# Data analysis

## Association test (logistic regression)

To do an association analysis in PLINK, use the code below. `--bfile` calls all three .fam, .bim, and .bed files with the same prefix. `--pheno` tells plink to read the phenotype from a specific file, in this 
case: beluga_pheno_sub.phe. `--pheno-name` tells plink which column of phenotype values to pull from within the phenotype file. `--glm` is the flag that will run the association test. It stands for Generalized
Linear Model, and when the phenotype data is in binary format (1/2/0), will perform a logistic regression (what we want). Following this flag with `allow-no-covars` means it can be run without adding a
covariate file. This command combined will produce a file labelled: "neonate_test.neonate_test.glm.logistic.hybrid". It contains the test results.

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


*Note*: use this code to count how many 0's, 1's and 2's for each phenotype column, where the number after `$` is the column number. 
```
awk 'NR>1{print $6}' beluga_pheno.phe | sort | uniq -c
```
