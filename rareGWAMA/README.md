# A rareGWAMA R package

**Table of Contents**

- [Introduction](#introduction)
- [Installation](#Installing-the-rareGWAMA-R-package)
- [Quick tutorial](#quick-tutorial)
    - [Single variant tests](#single-variant-tests)
    - [Conditional single variant tests](#Conditional-single-variant-tests)
- [Input files and arguments](#input-files)
    - [Score statistic files files (Summary statistics)](#Score-statistics-files)
    - [Imputation quality files](#Imputation-quality-files)
    - [Tabix range](#Tabix-range)
    - [Alternative](#Alternative)
- [Feedback/Contact](#Feedback/Contact)


## Introduction

rareGWAMA is a flexible software package for imputation based GWAS meta-analysis. 
Right now it is available as a beta version.


## Installing the rareGWAMA R package <a name="#Installing-the-rareGWAMA-R-package"></a>

The package is hosted on github, which allows installation and update to be very easy. First, make sure you have the `mvtnorm` and `data.table` packages installed:

    install.packages("devtools")

Then you could use:

    library(devtools)
    install_github("dajiangliu/rareGWAMA", subdir ="rareGWAMA")
    
With `library(rareGWAMA)`, your are ready to go!



## Feedback/Contact <a name="Feedback/Contact"></a>

Questions and requests can be sent to
Github issue page ([link](https://github.com/dajiangliu/rareGWAMA/issues))
or
Dajiang Liu ([dajiang.liu@outlook.com](mailto:dajiang.liu@outlook.com "mailto:dajiang.liu@outlook.com"))
