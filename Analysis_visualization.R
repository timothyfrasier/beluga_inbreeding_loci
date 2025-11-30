# GWAS analysis and visualization in R

## Download required packages
install.packages("tidyverse")
install.packages("qqman")
library(qqman)
library(tidyverse)
library(readr)
library(dplyr)

# Visualizing results of --glm logisitc regression test from PLINK

# Case 1: Neonate mortality, for allelic comparisons
neonate_test_allelic <- read_table("neonate_test_allelic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in the correct format (scaffold values given as 1,2,3, etc.):
neonate_test_allelic <- read_delim("neonate_test_allelic.txt") %>%
  filter(TEST == "ADD") %>%  # only includes values where the "additive" test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

neonate_test_allelic = neonate_test_allelic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
neonate_test_allelic$CHR <- as.numeric(neonate_test_allelic$CHR)


## Produce Manhattan plot
manhattan(neonate_test_allelic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)

# Case 1: Neonate mortality, for genotypic comparisons
neonate_test_genotypic <- read_table("neonate_test_genotypic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in correct format (scaffold values given as 1,2,3, etc.):
neonate_test_genotypic <- read_delim("neonate_test_genotypic.txt") %>%
  filter(TEST == "GENO_2DF") %>%  # only includes values where the genotypic test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

neonate_test_genotypic = neonate_test_genotypic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
neonate_test_genotypic$CHR <- as.numeric(neonate_test_genotypic$CHR)


## Produce Manhattan plot
manhattan(neonate_test_genotypic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)



# Case 2: Dystocia, for allelic comparisons
dystocia_test_allelic <- read_table("dystocia_test_allelic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in correct format (scaffold values given as 1,2,3, etc.):
dystocia_test_allelic <- read_delim("dystocia_test_allelic.txt") %>%
  filter(TEST == "ADD") %>%  # only includes values where the "additive" test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

dystocia_test_allelic = dystocia_test_allelic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
dystocia_test_allelic$CHR <- as.numeric(dystocia_test_allelic$CHR)


## Produce Manhattan plot
manhattan(dystocia_test_allelic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)

# Case 2: Dystocia, for allelic comparisons
dystocia_test_genotypic <- read_table("dystocia_test_genotypic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in correct format (scaffold values given as 1,2,3, etc.):
dystocia_test_genotypic <- read_delim("dystocia_test_genotypic.txt") %>%
  filter(TEST == "GENO_2DF") %>%  # only includes values where the genotypic test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

dystocia_test_genotypic = dystocia_test_genotypic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
dystocia_test_genotypic$CHR <- as.numeric(dystocia_test_genotypic$CHR)


## Produce Manhattan plot
manhattan(dystocia_test_genotypic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)




# Case 3: Neoplasia, for allelic comparisons
neoplasia_test_allelic <- read_table("neoplasia_test_allelic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in correct format (scaffold values given as 1,2,3, etc.):
neoplasia_test_allelic <- read_delim("neoplasia_test_allelic.txt") %>%
  filter(TEST == "ADD") %>%  # only includes values where the "additive" test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

neoplasia_test_allelic = neoplasia_test_allelic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
neoplasia_test_allelic$CHR <- as.numeric(neoplasia_test_allelic$CHR)


## Produce Manhattan plot
manhattan(neoplasia_test_allelic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)

# Case 3: Neoplasia, for genotypic comparisons
neoplasia_test_genotypic <- read_table("neoplasia_test_genotypic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in correct format (scaffold values given as 1,2,3, etc.):
neoplasia_test_genotypic <- read_delim("neoplasia_test_genotypic.txt") %>%
  filter(TEST == "GENO_2DF") %>%  # only includes values where the genotypic test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

neoplasia_test_genotypic = neoplasia_test_genotypic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
neoplasia_test_genotypic$CHR <- as.numeric(neoplasia_test_genotypic$CHR)


## Produce Manhattan plot
manhattan(neoplasia_test_genotypic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)





# Case 4: Infection, for allelic comparisons
neoplasia_test_allelic <- read_table("neoplasia_test_allelic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in correct format (scaffold values given as 1,2,3, etc.):
infection_test_allelic <- read_delim("infection_test_allelic.txt") %>%
  filter(TEST == "ADD") %>%  # only includes values where the "additive" test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

infection_test_allelic = infection_test_allelic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
infection_test_allelic$CHR <- as.numeric(infection_test_allelic$CHR)


## Produce Manhattan plot
manhattan(infection_test_allelic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)

# Case 4: Infection, for genotypic comparisons
infection_test_genotypic <- read_table("infection_test_genotypic.txt") %>% select(CHR, BP, SNP, P) # Read in the data (after renaming and saving as a .txt file)

## Modify data so that it is presented in correct format (scaffold values given as 1,2,3, etc.):
infection_test_genotypic <- read_delim("infection_test_genotypic.txt") %>%
  filter(TEST == "GENO_2DF") %>%  # only includes values where the genotypic test was run
  filter(is.finite(P), P > 0) %>% 
  select(CHR, BP, SNP, P)  # just chooses the four columns which the qqman package uses for Manhattan plots

infection_test_genotypic = infection_test_genotypic %>% mutate(CHR = case_when(
  CHR == "HiC_scaffold_1" ~ "1",
  CHR == "HiC_scaffold_2" ~ "2",
  CHR == "HiC_scaffold_3" ~ "3",
  CHR == "HiC_scaffold_4" ~ "4",
  CHR == "HiC_scaffold_5" ~ "5",
  CHR == "HiC_scaffold_6" ~ "6",
  CHR == "HiC_scaffold_7" ~ "7",
  CHR == "HiC_scaffold_8" ~ "8",
  CHR == "HiC_scaffold_9" ~ "9",
  CHR == "HiC_scaffold_10" ~ "10",
  CHR == "HiC_scaffold_11" ~ "11",
  CHR == "HiC_scaffold_12" ~ "12",
  CHR == "HiC_scaffold_13" ~ "13",
  CHR == "HiC_scaffold_14" ~ "14",
  CHR == "HiC_scaffold_15" ~ "15",
  CHR == "HiC_scaffold_16" ~ "16",
  CHR == "HiC_scaffold_17" ~ "17",
  CHR == "HiC_scaffold_18" ~ "18",
  CHR == "HiC_scaffold_19" ~ "19",
  CHR == "HiC_scaffold_20" ~ "20",
  CHR == "HiC_scaffold_21" ~ "21",
  CHR == "HiC_scaffold_22" ~ "22",))  # the package/code needs CHR column to be just numeric values
infection_test_genotypic$CHR <- as.numeric(infection_test_genotypic$CHR)


## Produce Manhattan plot
manhattan(infection_test_genotypic,
          xlab = "HiC Scaffold",
          ylim = c(0, 10), 
          cex = 0.4, 
          cex.axis = 0.9)
