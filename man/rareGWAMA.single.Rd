% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareGWAMA.single.R
\name{rareGWAMA.single}
\alias{rareGWAMA.single}
\title{single variant meta-analysis integrating imputation quality;}
\usage{
rareGWAMA.single(score.stat.file, imp.qual.file = NULL, tabix.range,
  alternative = "two.sided", col.impqual = 5, impQual.lb = 0.7,
  impQualWeight = FALSE, rmMultiAllelicSite = FALSE, gc = FALSE, ...)
}
\arguments{
\item{score.stat.file}{the file names of score statistic files;}

\item{imp.qual.file}{the file names of imputation quality;}

\item{tabix.range}{the tabix range. IT must be in quote and provided as a string;}

\item{alternative}{The alternative hypothesis. Default is two.sided;}

\item{col.impqual}{The column number for the imputation quality score;}

\item{impQual.lb}{The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;}

\item{impQualWeight}{Using imputation quality as weight}

\item{rmMultiAllelicSite}{Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function;}
}
\value{
A list of analysis results;
}
\description{
single variant meta-analysis integrating imputation quality;
}
\examples{
study.vec <- c("rareGWAMA/inst/extdata/study1.gz", "rareGWAMA/inst/extdata/study1.R2.gz", "rareGWAMA/inst/extdata/study3.gz");
r2.vec <- c("rareGWAMA/inst/extdata/study1.R2.gz", "rareGWAMA/inst/extdata/study2.R2.gz", "rareGWAMA/inst/extdata/study3.R2.gz");
res <- rareGWAMA.single(score.stat.file = study.vec, imp.qual.file = r2.vec ,"1:10177-57200", alternative="two.sided", 
                        col.impqual=5, impQual.lb=0, impQualWeight=FALSE, weight="Npq+impQ", 
                        gc=FALSE, rmMultiAllelicSite=TRUE);
head(res$res.formatted);
}
