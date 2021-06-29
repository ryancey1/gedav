---
title: "Homework 2"
author: "Ryan Yancey"
date: "29 June 2021"
output: 
    html_document:
        highlight: tango
        theme: cosmo
        keep_md: yes
colorlinks: yes
---



------------------------------------------------------------------------

#### For this assignment, we will be evaluating different normalization methods on 2-channel arrays in which 4 biological samples were run. The study is from GEO and the description of the experiment is provided as follows.

> **Series GSE12050**: Subcutaneous adipose tissue from lean and obese subjects ([source][1]).

#### Obtaining adipose tissue samples are paramount to the understanding of human obesity. We have examined the impact of needle-aspirated and surgical biopsy techniques on the study of subcutaneous adipose tissue (scAT) gene expression in both obese and lean subjects. Biopsy sampling methods have a significant impact on data interpretation and revealed that gene expression profiles derived from surgical tissue biopsies better capture the significant changes in molecular pathways associated with obesity. We hypothesize that this is because needle biopsies do not aspirate the fibrotic fraction of scAT; which subsequently results in an under-representation of the inflammatory and metabolic changes that coincide with obesity. This analysis revealed that the biopsy technique influences the gene expression underlying the biological themes commonly discussed in obesity (e.g. inflammation, extracellular matrix, metabolism, etc), and is therefore a caveat to consider when designing microarray experiments. These results have crucial implications for the clinical and physiopathological understanding of human obesity and therapeutic approaches. We will be working with 4 lean subjects from which a needle biopsy was taken.

------------------------------------------------------------------------

#### **1.) First load the marray library, then load the 4 GenePix files, making sure to extract the foreground and background median values from the Cy5 and Cy3 channels.**


```r
suppressPackageStartupMessages(library(marray))
```

First we need to decompress the `GSE12050_amend.zip` dataset downloaded into the "data" directory.


```r
dir(path = "data")
```

```
## [1] "gp"                 "GSE12050_amend.zip" "qRT-PCR.csv"
```

```r
# unzip and gunzip the GSE files
system(command = "unzip -o data/GSE12050_amend.zip -d data/gp; gunzip -q data/gp/*")

# load data
data <- read.GenePix(path = "data/gp")
```

```
## Reading ...  data/gp/GSM304445.gpr 
## Reading ...  data/gp/GSM304446.gpr 
## Reading ...  data/gp/GSM304447.gpr 
## Reading ...  data/gp/GSM304448.gpr
```

#### **2.) Normalize each array using median global, loess, and print-tip-group loess methods. Then plot MvA plots of all 4 arrays comparing no normalization to the other 3 normalization approaches.**

#### **3.) Plot density plots of the log ratio values for each normalization (and pre normalization) for only array #4. Put them all on the same plot. Make sure to label the axes and provide a legend.**

#### **4.) Based on the plots generated so far, which normalization do you think is most preferred for this dataset?**

#### **5.) Research has demonstrated that often a single channel, background subtracted provides as good a normalization as using both channels. To test this, we will be utilizing the fact that these 4 samples are replicates and calculate the correlation between them. So, first extract the Cy5 foreground and background values for each of the 4 arrays and subtract the background from the foreground values, then log2 transform these values. Then calculate global median normalization on these 4 arrays using these background subtracted Cy5 values. Hint, you need to use the median of each array to scale, such that after normalization, all arrays will have a median of 1.**

#### **6.) Next calculate a Spearman’s rank correlation between all 4 arrays that you normalized in #5 and do the same with the M values from loess normalized data that you generated in #2. Plot a scatter plot matrix for each of the two normalizations (pairs() function), and be sure to label the arrays and title the plot. Print the correlation coefficients to the screen.**

#### **7.) Now we want to compare these normalizations to quantile normalized data to see if we gain anything by leveraging the distributions across all 4 arrays. Carry out the steps in the lecture or use the paper from Bolstad et al. entitled: “A comparison of normalization methods for high density oligonucleotide array data based on variance and bias” (on the course website), but we are only going to conduct this on the Cy5 channel. The basic steps are as follows (these 6 steps are calculated on non-logged data; the data is logged after these steps are carried out):**

1. Subtract the foreground – background for each of the 4 chips for only the Cy5 channel. This should all be on the linear or raw scale (no logging yet).

2. Sort each column independently in this new matrix

3. Calculate row means for the sorted matrix

4. Create a new matrix with each row having the same values as the sorted row mean vectors from step #3 (you should have a new R matrix)

5. Rank the columns independently on the original background subtracted matrix (from step #1) Hint: use the `rank()` function with the argument `ties=”first”` or `order()`

6. Reorder the columns in the new matrix from step #4 using the ranks from step #5

**To verify that each array has the same distribution, use the `hist()` function to look at various arrays (e.g., `hist(c5.norm[,1])`; `hist(c5.norm[,2])`; etc.). Slight differences in distributions are a result of the ties in the ranking.**

#### **8.) Now log (base 2) the new R matrix you created from step 6 (question #7) and calculate a Spearman’s rank correlation between the 4 arrays and plot a scatter plot matrix as you did before. Print the correlation coefficients to the screen.**

#### **9.) Of the 4 normalization methods, which do you suggest as optimal and why?**

#### **10.) Now we want to work with a qRT-PCR dataset from patients with an inflammatory disease. The genes measured for this experiment included a set of proinflammatory chemokines and cytokines that are related to the disease. Download the raw qRT-PCR file called Inflammation_qRT-PCR.csv. Then change the normalization script from the lecture notes to include the housekeeping genes beta actin, GAPDH, and 18S. Look at the file to make sure the housekeepers are spelled correctly. Run the normalization script and output a data matrix of fold change values.**

#### **11.) Read the normalized qRT-PCR data matrix into R, using a Spearman’s rank correlation, which two patients are most correlated? Plot these two patients against each other in a scatter plot.**

[1]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12050
