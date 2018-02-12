#' Single variant meta-analysis
#'
#' @param score.stat.file files of score statistics
#' @param cov.file covariance matrix files
#' @param range tabix range of variants to be analyzed
#' @param refaltList A list of reference alternative and positions of variants to be analyzed; Each variant in the dataset will be match against ref/alt alleles specified in refaltList; Only variants with matched ref and alt alleles can be included; we also need AF and af.diff.max to determine if the flips are due to strand issues or due to ref/alt alleles flips;
#' @param alternative alternative hypothesis to be specified
#' @param callrate.cutoff Cutoffs of call rate, lower than which will NOT be analyzed (labelled as missing)
#' @param hwe.cutoff Cutoffs of HWE p-values
#' @param refaltList A list that contains pos, ref and alt for calibrating variant sites to have the right ref and alt alleles;
#' @param correctFlip Correcting for flipped alleles; Default is TRUE; If FALSE, studies with incorrect REF/ALT alleles will be labelled as missing, and dropped from meta-analyses
#' @param analyzeRefAltListOnly Only analyze variants that are included in the refaltList; Default is TRUE; If FALSE, variant sites in the dataset but not specified in the refaltList will be labelled as missing and dropped from studies;

#' @return a list consisting of results
#' @export
rareMETALS.single.group <- function(score.stat.file,cov.file,range,refaltList,alternative=c('two.sided','greater','less'),callrate.cutoff=0,hwe.cutoff=0,correctFlip=TRUE,analyzeRefAltListOnly=TRUE)
  {
    ix.gold <- 1;
    extra.par <- list(ix.gold=ix.gold,QC.par=list(callrate.cutoff=callrate.cutoff,hwe.cutoff=hwe.cutoff));
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, cov.file, range));
    if(length(raw.data.all)==0)
      return(list(list(p.value=NA,
                       skip=1,
                       statistic=NA,
                       no.var=0,
                       no.sample=NA)));
    if(length(alternative)>1) alternative <- alternative[1];
    ix.gold <- extra.par$ix.gold;
    if(length(extra.par$ix.gold)==0) ix.gold <- 1;
    raw.data <- raw.data.all[[1]];
    raw.data.ori <- raw.data;
    raw.data$cov <- NULL;
    ix.match <- match(raw.data$pos,refaltList$pos);
    refaltList <- list(pos=refaltList$pos[ix.match],ref=refaltList$ref[ix.match],alt=refaltList$alt[ix.match],af=refaltList$af[ix.match],anno=refaltList$anno[ix.match],af.diff.max=refaltList$af.diff.max,checkAF=refaltList$checkAF);
    
    if(length(extra.par$QC.par)>0)
        {
            raw.data <- QC(raw.data,extra.par$QC.par,cov=0);
        }
    log.mat <- raw.data$log.mat;
    ix.pop <- 1:length(raw.data$ref);
    U.stat <- 0;V.stat.sq <- 0;V.stat.sq <- 0;
    p.value <- 0;statistic <- 0;
    direction.by.study <- "+";
    beta1.est <- 0;
    beta1.sd <- 0;
    hsq.est <- 0;
    no.sample <- 0;
    maf.vec <- 0;
    no.sample.pop <- 0;
    no.sample.var <- 0;
    ustat <- 0;
    vstat <- 0;
    nref <- 0;nalt <- 0;nhet <- 0;
    no.sample.mat <- matrix(0,nrow=length(ix.pop),ncol=length(raw.data$nSample[[1]]));
    ix.var <- 1;nearby <- character(0);
    maf.byStudy <- NA;
    if(length(raw.data$ref[[ix.gold]])==0)
      return(list(p.value=NA,
                  ref=NA,
                  alt=NA,
                  raw.data=raw.data.ori,
                  raw.data.clean=raw.data,
                  statistic=NA,
                  direction.by.study=NA,
                  anno=NA,
                  maf=NA,
                  no.sample=NA,
                  beta1.est=NA,
                  beta1.sd=NA,
                  hsq.est=NA,
                  nearby=NA,
                  pos=NA));
    
    ref.gold <- refaltList$ref;
    alt.gold <- refaltList$alt;
    af.gold <- refaltList$af;
    anno.gold <- refaltList$anno;
    af.diff.max <- refaltList$af.diff.max;
    if(length(af.diff.max)==0) af.diff.max <- 0.7;
    maf.sd.vec <- 0;maf.maxdiff.vec <- 0;ix.maf.maxdiff.vec <- 0;
    maf.pop.ori <- 0;
    QC.by.study <- "";
    cochranQ.stat <- NA;cochranQ.pVal <- NA;cochranQ.df <- NA;I2 <- NA;
    for(ix.var in 1:length(raw.data$ref[[ix.gold]]))
      {
          cat("Analyzing ",raw.data$pos[ix.var],"\n");
          direction.by.study.var <- rep("X",length(ix.pop));
          U.stat <- 0;V.stat.sq <- 0;maf.pop <- 0;
          nref.var <- 0;nhet.var <- 0;nalt.var <- 0;
          ix.include <- rep(0,length(ix.pop));
          Z.byStudy <- 0;
          for(ii in 1:length(ix.pop))
              {
                  
                  if(length(raw.data$covXZ[[ii]])>0) {

                      warning(paste0("Study  ",ii," is analyzed as binary trait. It is advised to use rareMETALS2 for meta-analysis"))
                  }
                  res.flipAllele <- flipAllele(raw.data,raw.data.ori,refaltList,ii,ix.var,log.mat[ix.var,],correctFlip,analyzeRefAltListOnly);
                  raw.data <- res.flipAllele$raw.data;
                  ix.include <- res.flipAllele$ix.include;
                  log.mat[ix.var,] <- res.flipAllele$log.mat.var;
                  U.stat <- U.stat+rm.na(raw.data$ustat[[ii]][ix.var]);
                  V.stat.sq <- V.stat.sq+(rm.na(raw.data$vstat[[ii]][ix.var]))^2;
                  nref.var <- nref.var+rm.na(raw.data$nref[[ii]][ix.var]);
                  nalt.var <- nalt.var+rm.na(raw.data$nalt[[ii]][ix.var]);
                  nhet.var <- nhet.var+rm.na(raw.data$nhet[[ii]][ix.var]);
                  no.sample.mat[ii,ix.var] <- (raw.data$nSample[[ii]][ix.var]);
                  if(is.na(raw.data$ustat[[ii]][ix.var]))
                      {
                          direction.by.study.var[ii] <- "X";
                          ix.include[ii] <- 1;
                          log.mat[ix.var,ii] <- "ustatNA";
                          raw.data$nSample[[ii]][ix.var] <- NA;
                          raw.data$af[[ii]][ix.var] <- NA;
                          raw.data$ac[[ii]][ix.var] <- NA;
                          raw.data$ustat[[ii]][ix.var] <- NA;
                          raw.data$vstat[[ii]][ix.var] <- NA;
                          raw.data$nref[[ii]][ix.var] <- NA;
                          raw.data$nhet[[ii]][ix.var] <- NA;
                          raw.data$nalt[[ii]][ix.var] <- NA;
                      }
                  if(!is.na(raw.data$ustat[[ii]][ix.var]))
                      {
                          if(raw.data$ustat[[ii]][ix.var]>0) direction.by.study.var[ii] <- "+";
                          if(raw.data$ustat[[ii]][ix.var]<0) direction.by.study.var[ii] <- "-";
                          if(raw.data$ustat[[ii]][ix.var]==0) direction.by.study.var[ii] <- "=";
                      }
                  Z.byStudy[ii] <- (raw.data$ustat[[ii]][ix.var])/(raw.data$vstat[[ii]][ix.var]);
                  maf.pop[ii] <- ((raw.data$af[[ii]])[ix.var]);
                  maf.pop.ori[ii] <- (raw.data$af[[ii]])[ix.var];
              }
          
          maf.byStudy[ix.var] <- paste(maf.pop,collapse=",",sep=",");
          maf.vec[ix.var] <- sum(maf.pop*no.sample.mat[,ix.var],na.rm=TRUE)/sum(no.sample.mat[,ix.var],na.rm=TRUE);
          maf.sd.vec[ix.var] <- sqrt(sum((no.sample.mat[,ix.var])*((maf.pop.ori-maf.vec[ix.var])^2),na.rm=TRUE)/sum(no.sample.mat[,ix.var],na.rm=TRUE));
          maf.maxdiff.tmp <- max(abs(maf.pop.ori-maf.vec[ix.var]),na.rm=TRUE);
          maf.maxdiff.vec[ix.var] <- 0;ix.maf.maxdiff.vec[ix.var] <- sample(1:length(maf.vec),1);
          if(maf.maxdiff.tmp>0)
              {
                  maf.maxdiff.vec[ix.var] <- maf.maxdiff.tmp;
                  ix.maf.maxdiff.vec[ix.var] <- which.max(abs(maf.pop.ori-maf.vec[ix.var]));
              }
          ix.0 <- which(maf.pop==0 | maf.pop==1);
          if(length(ix.0)>0)
              {
                  maf.sd.vec[ix.var] <- sqrt(sum((no.sample.mat[-ix.0,ix.var])*((maf.pop[-ix.0]-maf.vec[ix.var])^2),na.rm=TRUE)/sum(no.sample.mat[,ix.var],na.rm=TRUE));
                  maf.maxdiff.vec[ix.var] <- 0;ix.maf.maxdiff.vec[ix.var] <- ix.0[1];
                  if(max(abs(maf.pop[-ix.0]-maf.vec[ix.var]),na.rm=TRUE)>0)
                      {
                          maf.maxdiff.vec[ix.var] <- max(abs(maf.pop[-ix.0]-maf.vec[ix.var]),na.rm=TRUE);
                          ix.rm0 <- (1:length(maf.pop))[-ix.0];
                          ix.maf.maxdiff.vec[ix.var] <- ix.rm0[which.max(abs(maf.pop[-ix.0]-maf.vec[ix.var]))];
                      }                  
              }
          direction.by.study[ix.var] <- paste(direction.by.study.var,sep='',collapse='');
          QC.by.study[ix.var] <- paste(log.mat[ix.var,],sep=",",collapse=",");

          beta1.est[ix.var] <- U.stat/V.stat.sq;
          beta1.sd[ix.var] <- sqrt(1/V.stat.sq);
          Z.meta <- U.stat/sqrt(V.stat.sq);
          cochranQ.stat[ix.var] <- sum((Z.byStudy-Z.meta)^2,na.rm=T);
          cochranQ.df[ix.var] <- sum(which(!is.na(Z.byStudy-Z.meta)^2))-1;
          if(cochranQ.df[ix.var]>0)
              cochranQ.pVal[ix.var] <- pchisq(cochranQ.stat[ix.var],df=cochranQ.df[ix.var],lower.tail=FALSE);
          if(cochranQ.df[ix.var]<=0)
              cochranQ.pVal[ix.var] <- NA;
          I2[ix.var] <- (cochranQ.stat[ix.var]-cochranQ.df[ix.var])/cochranQ.stat[ix.var];
          hsq.est[ix.var] <- (beta1.est[ix.var])*(beta1.est[ix.var])*V.stat.sq/sum(no.sample.mat[,ix.var],na.rm=TRUE);
          no.sample.var[ix.var] <- sum(no.sample.mat[,ix.var],na.rm=TRUE);
          nearby[ix.var] <- paste(get.gene.inWindow(raw.data$pos[ix.var],1e6),sep=",",collapse=",");
          if(alternative=='two.sided'){
              statistic[ix.var] <- U.stat^2/V.stat.sq;
              ustat[ix.var] <- U.stat;
              vstat[ix.var] <- sqrt(V.stat.sq);
          nref[ix.var] <- nref.var;
          nalt[ix.var] <- nalt.var;
          nhet[ix.var] <- nhet.var;
          ix.inf <- which(V.stat.sq==0);
          if(length(ix.inf)>0) statistic[ix.var][ix.inf] <- NA;
          p.value[ix.var] <- pchisq(statistic[ix.var],df=1,lower.tail=FALSE)
        }
        if(alternative=='greater'){
            statistic[ix.var] <- (U.stat/sqrt(V.stat.sq));
            p.value[ix.var] <- pnorm(statistic[ix.var],lower.tail=FALSE)
        }
        if(alternative=='less'){
            statistic[ix.var] <- (U.stat/sqrt(V.stat.sq));
            p.value[ix.var] <- pnorm(statistic[ix.var],lower.tail=TRUE)
        }
      }

    ref.out <- list(ref.gold);
    alt.out <- list(alt.gold);
    nSample.out <- list(no.sample.var);
    af.out <- list(maf.vec);
    ac.out <- list(as.integer(2*maf.vec*no.sample.var));
    callrate <- list(rep(1,length(ref.gold)));
    nref.out <- list(nref);
    nalt.out <- list(nalt);
    nhet.out <- list(nhet);
    ustat.out <- list(ustat);
    vstat.out <- list(vstat);
    callrate.out <- list(rep(1,length(maf.vec)));
    hwe.out <- list(rep(1,length(maf.vec)));
    effect.out <- list(beta1.est);
    pVal.out <- list(p.value);
    pos.out <- raw.data$pos;
    cov.out <- list(matrix(nrow=0,ncol=0));
    integratedData <- list();
    integratedData[[1]] <- list(ref=ref.out,
                                alt=alt.out,
                                nSample=nSample.out,
                                af=af.out,
                                ac=ac.out,
                                callrate=callrate.out,
                                hwe=hwe.out,
                                nref=nref.out,
                                nalt=nalt.out,
                                ustat=ustat.out,
                                vstat=vstat.out,
                                nhet=nhet.out,
                                effect=effect.out,
                                pVal=pVal.out,
                                cov=cov.out,
                                pos=pos.out,
                                covZZ=list(matrix(nrow=0,ncol=0)),
                                covXZ=list(matrix(ncol=0,nrow=length(ref.out[[1]]))),
                                hweCase=list(rep(NA,length(ref.out[[1]]))),
                                hweCtrl=list(rep(NA,length(ref.out[[1]]))),
                                afCtrl=list(rep(NA,length(ref.out[[1]]))),
                                afCase=list(rep(NA,length(ref.out[[1]]))));
    
    
    return(list(p.value=p.value,
                ref=ref.gold,
                alt=alt.gold,
                integratedData=integratedData,
                raw.data=raw.data.ori,
                clean.data=raw.data,
                statistic=statistic,
                direction.by.study=direction.by.study,
                anno=raw.data$anno,
                maf=maf.vec,
                maf.byStudy=maf.byStudy,
                maf.maxdiff.vec=maf.maxdiff.vec,
                ix.maf.maxdiff.vec=ix.maf.maxdiff.vec,
                maf.sd.vec=maf.sd.vec,
                no.sample.mat=no.sample.mat,
                no.sample=no.sample.var,
                beta1.est=beta1.est,
                beta1.sd=beta1.sd,
                QC.by.study=QC.by.study,
                hsq.est=hsq.est,
                nearby=nearby,
                cochranQ.stat=cochranQ.stat,
                cochranQ.df=cochranQ.df,
                cochranQ.pVal=cochranQ.pVal,
                I2=I2,
                log.mat=log.mat,
                pos=raw.data$pos));
  }
