# A rareGWAMA R package

**Table of Contents**

- [Introduction](#introduction)
- [Installation](#Installing-the-rareGWAMA-R-package)
- [Quick tutorial](#quick-tutorial)
    - [Single variant tests](#single-variant-tests)
    - [Conditional single variant tests](#conditional-single-variant-tests)
- [Input files and arguments](#input-files)
    - [Score statistic files files (Summary statistics)](#Score-statistics-files)
    - [Imputation quality files](#Imputation-quality-files)
    - [Tabix range](#Tabix-range)
    - [Alternative](#Alternative)
- [Feedback/Contact](#Feedback/Contact)


## Introduction

rareGWAMA is a flexible software package for imputation based GWAS meta-analysis. 
Right now it is available as a beta version.


## Installing the rareGWAMA R package <a name="Installing-the-rareGWAMA-R-package"></a>

The package is hosted on github, which allows installation and update to be very easy. First, make sure you have the `mvtnorm` and `data.table` packages installed:

    install.packages("devtools")

Then you could use:

    library(devtools)
    install_github("dajiangliu/rareGWAMA", subdir ="rareGWAMA")
    
With `library(rareGWAMA)`, your are ready to go!


## Quick tutorial <a name="quick-tutorial"></a>

### Single variant tests <a name="conditional-single-variant-tests"></a>
`rareGWAMA.single(study, study.imp, "1:11000-58000", alternative="two.sided", col.impqual=5, impQual.lb=0, impQualWeight=FALSE, weight="Npq+impQ",gc=FALSE, rmMultiAllelicSite=TRUE);`

### Conditional single variant tests <a name="single-variant-tests"></a>



## Feedback/Contact <a name="Feedback/Contact"></a>

Questions and requests can be sent to
Github issue page ([link](https://github.com/dajiangliu/rareGWAMA/issues))
or
Dajiang Liu ([dajiang.liu@outlook.com](mailto:dajiang.liu@outlook.com "mailto:dajiang.liu@outlook.com"))
