# Data analysis

## Association test (logistic regression)

To perform an association analysis in PLINK, use the code below. `--bfile` calls all three .fam, .bim, and .bed files with the same prefix. `--pheno` tells PLINK to read the phenotype from a specific file, in this 
case, *pheno_sub.phe*. `--pheno-name` tells PLINK which column of phenotype values to pull from within the phenotype file. `--glm` is the flag that will run the association test. It stands for Generalized
Linear Model, and will perform a logistic regression when the phenotype data is in binary format (1/2/0). Following this flag with `allow-no-covars` means it can be run without adding a
covariate file. Including the modifier 'genotypic' after '--glm' adds a dominance-deviation column within the results file to compare genotypic associations. 

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


After saving as a .txt file, rename columns 1, 2 and 3 as *CHR*, *BP* and *SNP*, respectively. This is so that the important information is appropriately formatted to be used in Rstudio with the qqman package. 


To make files containing significant SNPS, you can use the code below and specify the p-value threshold. The command `awk` filters/processes columns in text files. `'NR==1` tells it to keep the header upon creating the new file. Here, `$16` specifies which column to process the values from. 

```
awk 'NR>1 && $16 < 0.1 {print $3}' neonate_test_allelic.txt > snps_0.1_nta.txt
```

Saving significant SNPs for genotypic association tests requires two steps: first, the dominance-deviation column is created into a new file, and then similar steps as above are run.

```
grep "GENO_2DF" dystocia_test_genotypic.txt > dystocia_test_geno2df.txt
awk 'NR>1 && $16 < 0.1 {print $3}' dystocia_test_geno2df.txt > snps_0.1_dg.txt
```

## Epistatic interactions - data prep

To identify potential clusters of alleles based on epistatic interactions in Rstudio, we first need to produce a .tped file from the filtered bfiles. `--recode` creates a new fileset. `--recode12` converts alphabetical allele values to 1/2, where 
1 and 2 are the first and second alleles in the .bim, respectively. `--recode transpose` causes the columns and rows to interchange/swap. This step is necessary since we need each SNP to be a singular row. We can also do additional linkage disequilibrium pruning using `--indep-pairwise 50 5 0.2`, where the arguments 50, 5, and 0.2 signify the window size, step size, and r<sup>2</sup> threshold, respectively. 
These steps will have to be run for each cause of death (here nta represents neonate (nt) mortality, with data pulled from the allelic GWAS test results (a)).

```
./plink --bfile filtered_data --extract snps_0.1_nta.txt --recode12 --out filtered12_0.1_nta --allow-extra-chr
./plink --file filtered12_0.1_nta --recode transpose --out tfiltered12_0.1_nta --allow-extra-chr
```


## The code below can be used to determine if any of the sequenced SNPs are within protein-coding regions.

This requires a reference genome, which can be downloaded from DNA Zoo. The one used here is titled *ASM228892v2_HiC.fasta_v2.functional.gff3*. First convert the .map file (epistatic input) to a .bed file.

```
awk '{OFS="\t"; print $1, $4-1, $4, $2}' filtered12_0.1_npg.map > bedtools_input_npg.bed
```

Make a .bed file from the reference .gtf file containing SNPs within gene regions. 

```
awk '$3 == "gene"' ASM228892v2_HiC.fasta_v2.functional.gff3 > beluga_reference.gff3.bed 
```

Intersect the files. Here `wc -l` counts how many SNPs overlap. To save these SNPs to a new file use `>` followed by the name of the new file. 

```
bedtools intersect -a bedtools_input_npg.bed -b beluga_reference.gff3.bed -wa | wc -l
bedtools intersect -a bedtools_input_npg.bed -b beluga_reference.gff3.bed -wa > intersect_npg.bed
```



