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
  1 -> # chromosome
  2 -> position
  3 -> ID
  4 -> reference allele
  5 -> alternate allele
  7 -> filter
  9 -> format (GT:DP:AD:GQ:GL:VAF:VAF1)
  10 -> sample

Column 8 is also very useful as it provides a lot of information, such as frequency of the alternate allele, etc.

  *Reminder*: hemizgyous = only one allele at a site

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

Pressing `Ctrl+C` cancels the command that it is processing, for when it gets stuck.

## Mortality data

Example of how to count the number of rows where the category for mortality is infectious disease (*note that row 13 is "Cateogry"*):

```
cut -f 13 beluga_mortality_2021.txt | grep "Infectious" | wc -l
```

This counts the total number of rows, minus the header (for this data set):

```
tail -n +2 beluga_mortality_2021.txt | wc -l
```
