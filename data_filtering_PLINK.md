# Data filtering

## Quality control

Before running commands to filter/quality check, run the below line before and after to count the number of variants and number of individuals/samples. This should tell you 
if alterations remove data or not (using `.bim` will give you the number of variants and `.fam` will give you the number of samples). 

```
wc -l 0.01_ff.bim 
wc -l 0.01_fully_filtered.fam
```


Next, we will run the code below. `--geno 0.05` filters out variants with high rates of missing genotypes, while `--mind 0.05` does the same but for samples. Using 0.05 specifies to remove variants or samples
with missing genotype call rates > 5%. `--hwe 1e-6` filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold. Lastly, `--maf 0.01` filters out all variants 
with minor allele frequency below the provided threshold. Doing these one at a time, along with `--make-bed` and `--out`, will allow you to observe the effect of these commands on the data set. Example below:

```
./plink --bfile 0.01_ff --geno 0.05 --make-bed --out qc_data --allow-extra-chr
```


This can also be completed in one command, and the output will give you how many variants and samples were removed, in addition to what remains:

```
./plink -bfile 0.01_ff --geno 0.05 --mind 0.05 --hwe 1e-6 --maf 0.01 --make-bed --out qc_data --allow-extra-chr
```


Next, we can check for relatedness using the `king-cutoff` command and a threshold value. The threshold value will be important for what degree of relatedness you want to remove. For our case, we will use first
degree relatives, equivalent to a threshold value of **0.177**. This command does not automatically remove samples, and therefore requires two commands to first identify samples that should be removed (which are
saved in a `.king.cutoff.out.id` file), and a subsequent command to actually remove them and create a new set of bfiles. Also, the flag `king-cutoff` requires plink2.

```
./plink2 --bfile qc_data --make-king-table --king-cutoff 0.177 --out king_cutoff
./plink2 --bfile qc_data --remove king_cutoff.king.cutoff.out.id --make-bed --out filtered_data
```


So after quality checks, our data is now stored in a set of bfiles with the prefix *filtered_data*.