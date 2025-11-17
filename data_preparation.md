# Data preparation

# Practice
##  Styling text

This is **bold**. This is *italics*. This is <ins>underlined</ins>. 

Below is a line of code - a very useful one. 
```
ls
```

***This is bold and italics, so obviously super important***. 

If I wanted to show a square meter I could use a superscript like this: m<sup>2</sup>. For carbon dioxide, I could instead use subscript 
notation like this: CO<sub>2</sub>.

##  Useful Git commands

Some useful Git commands include:
- `git status` to list recent files that have not yet been committed.
- ` git add` to add new files.
- `git commit` to commit these files to Git hub.

# Beginning to look at data

My data is located within **"0.01_fully_filtered.vcf"**. We are using 0.01, which means we have
99% confidence that each base is correct. 

To show the first 30 lines within the data frame, we will use `gunzip` or `gzip` and `-c`, which means it will not truly unzip the file, 
but show the unzipped file on the screen. Along with the data set name, we can then pipe `|` that into the `head` command and specify `-n 30` 
(or any number of lines). 

```
gunzip -c 0.01_fully_filtered.vcf.gz | head -n 30
```


This shows the first 30 rows of info, where each row begins with ##. To count the number of lines that start with ##, we can use `zgrep` and 
`-c`, which means to count. For example:

```
zgrep -c "^##" 0.01_fully_filtered.vcf.gz
```


To show how many lines **DO NOT** start with ##, we can use `-v`. Below will show the first row (hence 1) that does not start with ##.

```
zgrep -v "^##" 0.01_fully_filtered.vcf.gz | head -n 1
```


Use `cut` to look at just some of the important column headers. For example:

```
zgrep -v "^##" 0.01_fully_filtered.vcf.gz | cut -f 1,2,4,5,9,10 | head -n 2
```


These columns show:
  - 1 -> # chromosome
  - 2 -> position
  - 3 -> ID
  - 4 -> reference allele
  - 5 -> alternate allele
  - 7 -> filter
  - 9 -> format (GT:DP:AD:GQ:GL:VAF:VAF1)
  - 10 -> sample

Column 8 is also very useful as it provides a lot of information, such as frequency of the alternate allele, etc.

  - *Reminder*: hemizgyous = only one allele at a site

Use the below code to count the total number of columns in the header, where `^#CHROM` will tell it to only choose
this row and `cw -w` asks it to count the words (number of columns in this case). Given that there are 9 rows of metadata,
the total number of rows minus 9 will give you the number of samples.

```
zgrep "^#CHROM" 0.01_fully_filtered.vcf.gz | cw -w
```


To count the number of rows starting with #, use the command below, where `cw -1` asks it to count the number of lines:

```
zgrep -v "^#" 0.01_fully_filtered.vcf.gz | wc -l
```

Pressing `Ctrl+C` cancels the command that it is processing. Use it when command line gets stuck.


## Mortality data

Example of how to count the number of rows where the category for mortality is infectious disease (*note that row 13 is "Cateogry"*):

```
cut -f 13 beluga_mortality_2021.txt | grep "Infectious" | wc -l
```


This counts the total number of rows, minus the header (for this data set):

```
tail -n +2 beluga_mortality_2021.txt | wc -l
```

# PLINK-related

Note: this is a useful site for PLINK flags/commands: https://www.cog-genomics.org/plink/1.9/input#vcf

To learn more about plink flags/commands, use the following code (for some reason on my computer plink commands need to include "./"):

```
./plink --help
```

I started by trying to convert the vcf file into plink format. I was at first having issues, but the code below seems to have worked.`Gunzip -c`
is required for plink since we are giving it a zipped file. `--vcf` specifies that the following commands are to be pulled from the vcf file,
while `/dev/stdin` tells it to pull from the unzipped file (standard input). `--make-bed` tells it to generate a new PLINK binary fileset (.bed,
.bim and .fam). `--out` is used to tell it what prefix to name these new files. `--keep-allele-order` is used so that PLINK does not reorganize
the reference and alternate alleles (very important). Command line output told me to use `--allow-extra-chr` (but not entirely sure about this). 

```
gunzip 0.01_fully_filtered.vcf.gz | ./plink2 --vcf /dev/stdin \
  --make-bed \
  --out converted_data \
  --keep-allele-order \
  --allow-extra-chr
```
***^ using this code results in messed-up .bed data

The code below works the same, except it uses plink2 and does not require `--keep-allele-order` or `--allow-extra-chr` (can choose different prefix name if you like).
```
./plink2 --vcf 0.01_fully_filtered.vcf.gz \
--make-bed \
--out converted_data
```

Each file in this PLINK fileset stores different genetic data:
- **.bim** contains variant (SNP) data (important info here is base-pair coordinate and alt/ref allele)
- **.fam** contains sample information, with 6 columns as follows:
  - 1. Family ID
  - 2. Within-family ID
  - 3. Parental ID
  - 4. Maternal ID
  - 5. Sex code ('1' = male, '2' = female, '0' = unknown)
  - 6. Phenotype value (-9 or 0 if unknown)
- **.bed** contains genotype data 
