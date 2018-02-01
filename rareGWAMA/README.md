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

1.The very basic test is using:  
`res <- rareGWAMA.single(study, imp.qual, "1:11000-58000", alternative="two.sided", col.impqual=5, impQual.lb=0, impQualWeight=FALSE, weight="Npq+impQ",gc=FALSE, rmMultiAllelicSite=TRUE);`   
please find more detail in the [input and arguments part](#Input-files-and-arguments) for the arguments:
* study: The file names of score statistic files, which could be a **vector object**;
* imp.qual: The file names of imputation quality, which could be a **vector object**;
* "1:11000-58000": The tabix range, which must be in quote and provided as a string like this;
* alternative: The alternative hypothesis. Default is two.sided;
* col.impqual: The column number for the imputation quality score;
* impQual.lb: The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
* impQualWeight: Using imputation quality as weight;
* rmMultiAllelicSite: Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function;  

2.The out put should be as follows:  
`head(res$res.formatted))`
```
     POS           REF ALT  AF         STAT       PVALUE     BETA        SD         N        DIRECTION                          EFFECTIVE_N
[1,] "9:100000040" "G" "A"  "2.28e-05" "1.01e+00" "3.16e-01" " 7.29e-01" " 0.72624" "41634"  "XXXXXXXXXXXXXXXXXXXXXXXXXXXX++XX" "7925"
[2,] "9:100000056" "T" "TA" "7.07e-01" "1.62e+00" "2.04e-01" "-9.09e-03" " 0.00715" "47219"  "XXXXXXXXXXXXXXXX-XXXXXXXXXXX--XX" "46951"
[3,] "9:100000102" "G" "C"  "6.91e-04" "8.00e-01" "3.71e-01" "-7.78e-02" " 0.08691" "95869"  "+++---+-------++XX-++--X++----XX" "36378"
[4,] "9:100000172" "C" "T"  "2.26e-05" "2.25e-01" "6.35e-01" "-1.99e-01" " 0.41982" "125656" "XXXXXXXXXXXXXXXX-XXXXXXXXXXX--X+" "20434"
[5,] "9:100000177" "C" "G"  "2.16e-04" "8.38e-01" "3.60e-01" "-1.04e-01" " 0.11335" "179891" "-+---X+XX-+--+-+-X-+-+-X++-+--X+" "89179"
[6,] "9:100000187" "A" "T"  "1.13e-04" "6.78e-01" "4.10e-01" "-2.35e-01" " 0.28514" "54235"  "+-++-XX+--++--+-XX+++++X----XXXX" "21285"
```

### Conditional single variant tests <a name="single-variant-tests"></a>


## Input files and arguments <a name="input-files"></a> 

Both score statistics file and imputation quality file should be tabixed.

### Score statistics files:  
if you use [RVTEST](https://github.com/zhanxw/rvtests), the output is ready to go:
```
CHROM   POS     REF     ALT     N_INFORMATIVE   AF      INFORMATIVE_ALT_AC      CALL_RATE       HWE_PVALUE      N_REF   N_HET   N_ALT   U_STAT  SQRT_V_STAT     ALT_EFFSIZE     PVALUE
1       10177   A       AC      2352    0.5     2352    1       0       0       2352    0       1.67496 2.51553 0.264695        0.505508
1       10235   T       TA      2352    0       0       1       1       2352    0       0       -0.108472       0.207841        -2.51104        0.601742
1       10352   T       TA      2352    0.5     2352    1       0       0       2352    0       0.665562        2.61389 0.0974122       0.799013
1       10539   C       A       2352    0       0       1       1       2352    0       0       -0.00020902     0.0626305       -0.0532862      0.997337
```

###Imputation quality files:

```
CHROM   POS     REF     ALT     Rsq
1       10177   A       AC      0.00581
1       10235   T       TA      0.00396
1       10352   T       TA      0.00608
1       10539   C       A       0.00154
1       10616   CCGCCGTTGCAAAGGCGCGCCG  C       0.02085
1       10642   G       A       0.00013
1       11008   C       G       0.01251
1       11012   C       G       0.01252
```

## Feedback/Contact <a name="Feedback/Contact"></a>

Questions and requests can be sent to
Github issue page ([link](https://github.com/dajiangliu/rareGWAMA/issues))
or
Dajiang Liu ([dajiang.liu@outlook.com](mailto:dajiang.liu@outlook.com "mailto:dajiang.liu@outlook.com"))
