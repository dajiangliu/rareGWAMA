#' This is the function for flipping alleles
#'
#' @param raw.data The input datasets to be considered flipped
#' @param raw.data.ori The input datasets to be considered flipped
#' @param refaltList The list consists of ref, alt, pos, af and af.diff.max, as well as the option of whether throw away sites with large af.differences checkAF;
#' @param ix.pop index of the population
#' @param ix.var index of the variant;
#' @param log.mat.var The log for QC procedure;
#' @param correctFlip Correct for score and covariance matrices for flipped alleles;
#' @return A list consist of modified raw.data, ix.include and log.mat.var
#' @export
flipAllele <- function(raw.data,raw.data.ori,refaltList,ix.pop,ix.var,log.mat.var,correctFlip=TRUE,analyzeRefAltListOnly=TRUE)
    {
        ##########################print(c("REF","ALT",refaltList$ref[ix.var],refaltList$alt[ix.var]))
        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & analyzeRefAltListOnly )
            {
                ii <- ix.pop;
                
                ix.include <- rep(0,length(raw.data$ustat));
                
                if(length(raw.data$cov)>0)
                    {
                        raw.data$cov[[ii]][ix.var,] <- NA;
                        raw.data$cov[[ii]][,ix.var] <- NA;
                    }
                log.mat.var[ii] <- "NotInRefAltList";
                raw.data$nSample[[ii]][ix.var] <- NA;
                raw.data$af[[ii]][ix.var] <- NA;
                raw.data$ac[[ii]][ix.var] <- NA;
                raw.data$ustat[[ii]][ix.var] <- NA;
                raw.data$vstat[[ii]][ix.var] <- NA;
                raw.data$nref[[ii]][ix.var] <- NA;
                raw.data$nhet[[ii]][ix.var] <- NA;
                raw.data$nalt[[ii]][ix.var] <- NA;                
                return(list(raw.data=raw.data,
                            log.mat.var=log.mat.var,
                            ix.include=ix.include));
            }

        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & !analyzeRefAltListOnly )
            {
                ii <- ix.pop;                
                ix.include <- rep(0,length(raw.data$ustat));
                ##########################print("NotInRefAltList");
                log.mat.var[ii] <- "NotInRefAltList";
                return(list(raw.data=raw.data,
                            log.mat.var=log.mat.var,
                            ix.include=ix.include));
            }


        
        ii <- ix.pop;
        ref.gold <- refaltList$ref;alt.gold <- refaltList$alt;af.gold <- refaltList$af;af.diff.max <- refaltList$af.diff.max;checkAF <- refaltList$checkAF;
        if(length(checkAF)==0) checkAF <- FALSE;
        
        af.diff <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
        ix.include <- rep(0,length(raw.data$ustat));
        if(is.na(af.diff)) af.diff <- 0;
        if(!is.na(raw.data.ori$ref[[ii]][ix.var]) | !is.na(raw.data.ori$alt[[ii]][ix.var]))
            {
                if(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1)
                    {

                        if(af.diff<=af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$vstat[[ii]][ix.var] <- 0;
                                log.mat.var[ii] <- "Monomorphic";
                            }
                        if(af.diff>=af.diff.max)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                                    }
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$vstat[[ii]][ix.var] <- 0;
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";

                            }
                    }           
            }
        
        if(!is.na(raw.data$ref[[ii]][ix.var]) & !is.na(raw.data$alt[[ii]][ix.var]) & !(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1))
            {
                strandAmbiguous <- (((ref.gold[ix.var]=="A") & (alt.gold[ix.var]=="T")) | ((ref.gold[ix.var]=="T") & (alt.gold[ix.var]=="A")) | ((ref.gold[ix.var]=="C") & (alt.gold[ix.var]=="G")) | ((ref.gold[ix.var]=="G") & (alt.gold[ix.var]=="C")));
                flip.ref.alt <- (raw.data$ref[[ii]][ix.var]==refaltList$alt[ix.var] & raw.data$alt[[ii]][ix.var]==refaltList$ref[ix.var]);
                match.ref.alt <- (refaltList$ref[ix.var]==(raw.data$ref[[ii]][ix.var]) & refaltList$alt[ix.var]==raw.data$alt[[ii]][ix.var]);
                ##mono <- (((raw.data$ref[[ii]][ix.var]==".") | (raw.data$ref[[ii]][ix.var]==0) | (raw.data$alt[[ii]][ix.var]==".") | (raw.data$alt[[ii]][ix.var]==0) | (rm.na(raw.data$af[[ii]][ix.var])==0)| rm.na(raw.data$af[[ii]][ix.var])==1) & !is.na(raw.data.ori$nSample[[ii]][ix.var]))
                mono <- (((raw.data$ref[[ii]][ix.var]==".") | (raw.data$ref[[ii]][ix.var]==0) | (raw.data$alt[[ii]][ix.var]==".") | (raw.data$alt[[ii]][ix.var]==0) ) & (rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1) & !is.na(raw.data.ori$nSample[[ii]][ix.var]))
                if(!match.ref.alt & !correctFlip &!mono)
                    {
                        
                        ix.include[ii] <- 1;
                        if(length(raw.data$cov)>0)
                            {
                                raw.data$cov[[ii]][ix.var,] <- NA;
                                raw.data$cov[[ii]][,ix.var] <- NA;
                            }
                        
                        log.mat.var[ii] <- "MismatchRemove";
                        
                        raw.data$nSample[[ii]][ix.var] <- NA;
                        raw.data$af[[ii]][ix.var] <- NA;
                        raw.data$ac[[ii]][ix.var] <- NA;
                        raw.data$ustat[[ii]][ix.var] <- NA;
                        raw.data$vstat[[ii]][ix.var] <- NA;
                        raw.data$nref[[ii]][ix.var] <- NA;
                        raw.data$nhet[[ii]][ix.var] <- NA;
                        raw.data$nalt[[ii]][ix.var] <- NA;
                    }
                
                if(flip.ref.alt & (!strandAmbiguous))
                    {
                        raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                        raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                        if(length(raw.data$cov)>0)
                            {
                                raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                            }
                        
                        nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                        nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                        nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                        raw.data$nref[[ii]][ix.var] <- nref.tmp;
                        raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                        raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                        ix.include[ii] <- 1;
                     
                        log.mat.var[ii] <- "FlipRefAlt";

                    }
                if(match.ref.alt & (!strandAmbiguous))
                    {
                        ix.include[ii] <- 1;
                        log.mat.var[ii] <- "Match";
                    }
                af.diff.min <- 0.05;
                if(flip.ref.alt & strandAmbiguous)
                    {                        
                        if(af.diff<=af.diff.min)
                            {
                                ix.include[ii] <- 1;
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- NA;
                                        raw.data$cov[[ii]][,ix.var] <- NA;
                                    }
                                
                                log.mat.var[ii] <- "FlipStrand";
                                
                                raw.data$nSample[[ii]][ix.var] <- NA;
                                raw.data$af[[ii]][ix.var] <- NA;
                                raw.data$ac[[ii]][ix.var] <- NA;
                                raw.data$ustat[[ii]][ix.var] <- NA;
                                raw.data$vstat[[ii]][ix.var] <- NA;
                                raw.data$nref[[ii]][ix.var] <- NA;
                                raw.data$nhet[[ii]][ix.var] <- NA;
                                raw.data$nalt[[ii]][ix.var] <- NA;
                                
                                
                            }
                        if(af.diff>af.diff.min)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                                if(length(raw.data$cov)>0)
                                    {
                                        
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]

                                    }
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                                ix.include[ii] <- 1;
                                
                                log.mat.var[ii] <- "FlipRefAlt";
                            }
                    }

                if(match.ref.alt & strandAmbiguous)
                    {
                        if(af.diff<=af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                
                                log.mat.var[ii] <- "Match";
                            }
                        if(af.diff>af.diff.max)
                            {
                                log.mat.var[ii] <- "FlipStrand";
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- NA;
                                        raw.data$cov[[ii]][,ix.var] <- NA;
                                    }
                                
                                raw.data$nSample[[ii]][ix.var] <- NA;
                                raw.data$af[[ii]][ix.var] <- NA;
                                raw.data$ac[[ii]][ix.var] <- NA;
                                raw.data$ustat[[ii]][ix.var] <- NA;
                                raw.data$vstat[[ii]][ix.var] <- NA;
                                raw.data$nref[[ii]][ix.var] <- NA;
                                raw.data$nhet[[ii]][ix.var] <- NA;
                                raw.data$nalt[[ii]][ix.var] <- NA;


                            }
                    }
                
                if(mono)
                    {
                        if(af.diff<=af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";
                                raw.data$ustat[[ii]][ix.var] <- 0;
                            }
                        if(af.diff>=af.diff.max)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                                    }
                                
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";

                            }                        
                    }
                if(!match.ref.alt & !flip.ref.alt & !mono)
                    {
                        ix.include[ii] <- 1;
                        if(length(raw.data$cov)>0)
                            {
                                raw.data$cov[[ii]][ix.var,] <- NA;
                                raw.data$cov[[ii]][,ix.var] <- NA;
                            }
                        
                        log.mat.var[ii] <- "Unmatched";
                        raw.data$nSample[[ii]][ix.var] <- NA;
                        raw.data$af[[ii]][ix.var] <- NA;
                        raw.data$ac[[ii]][ix.var] <- NA;
                        raw.data$ustat[[ii]][ix.var] <- NA;
                        raw.data$vstat[[ii]][ix.var] <- NA;
                        raw.data$nref[[ii]][ix.var] <- NA;
                        raw.data$nhet[[ii]][ix.var] <- NA;
                        raw.data$nalt[[ii]][ix.var] <- NA;
                    }
                if(checkAF==TRUE)
                    {
                        ########################print(c("raw.data$af[[ii]][ix.var]",raw.data$af[[ii]][ix.var]));
                        ########################print(c("af.gold[ix.var]",af.gold[ix.var]));
                        af.diff.new <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
                        if(is.na(af.diff.new)) af.diff.new <- 0;
                        if(af.diff.new>af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "DiffAF";
                                raw.data$nSample[[ii]][ix.var] <- NA;
                                raw.data$af[[ii]][ix.var] <- NA;
                                raw.data$ac[[ii]][ix.var] <- NA;
                                raw.data$ustat[[ii]][ix.var] <- NA;
                                raw.data$vstat[[ii]][ix.var] <- NA;
                                
                                raw.data$nref[[ii]][ix.var] <- NA;
                                raw.data$nhet[[ii]][ix.var] <- NA;
                                raw.data$nalt[[ii]][ix.var] <- NA;                                
                            }
                    }
     
            }
        return(list(raw.data=raw.data,
                    log.mat.var=log.mat.var,
                    ix.include=ix.include));
    }
