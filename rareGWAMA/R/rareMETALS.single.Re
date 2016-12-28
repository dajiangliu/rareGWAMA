#' Single variant meta-analysis
#'
#' @param score.stat.file files of score statistics
#' @param cov.file covariance matrix files
#' @param range tabix range of variants to be analyzed
#' @param alternative alternative hypothesis to be specified
#' @param ix.gold Gold standard population to align reference allele to
#' @param callrate.cutoff Cutoffs of call rate, lower than which will NOT be analyzed (labelled as missing)
#' @param hwe.cutoff Cutoffs of HWE p-values
#' @return a list consisting of results
#' @export
rareMETALS.single <- function(score.stat.file,cov.file,range,alternative=c('two.sided','greater','less'),ix.gold=1,callrate.cutoff=0,hwe.cutoff=0)
  {
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
    if(length(extra.par$QC.par)>0) 
      raw.data <- QC(raw.data,extra.par$QC.par,cov=0);
    log.mat <- raw.data$log.mat;
    ix.pop <- 1:length(raw.data$ref);
    ref.gold <- raw.data$ref[[ix.gold]];
    alt.gold <- raw.data$alt[[ix.gold]];
    U.stat <- 0;V.stat.sq <- 0;V.stat.sq <- 0;
    p.value <- 0;statistic <- 0;
    direction.by.study <- "+";
    QC.by.study <- "";
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
    if(length(raw.data$ref[[ix.gold]])==0)
      return(list(p.value=NA,
                  ref=NA,
                  alt=NA,
                  raw.data=raw.data,
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
    for(ix.var in 1:length(raw.data$ref[[ix.gold]]))
      {
        direction.by.study.var <- rep("?",length(ix.pop));
        U.stat <- 0;V.stat.sq <- 0;maf.pop <- 0;
        nref.var <- 0;nhet.var <- 0;nalt.var <- 0;
        for(ii in 1:length(ix.pop)) {
            if(length(raw.data$covXZ[[ii]])>0) {
                warning(paste0("Study  ",ii," is analyzed as binary trait. It is advised to use rareMETALS2 for meta-analysis"))
            }
              
            if(is.na(ref.gold[ix.var]) & !is.na(raw.data$ref[[ii]][ix.var]))
              {
                ref.gold[ix.var] <- raw.data$ref[[ii]][ix.var];
                alt.gold[ix.var] <- raw.data$alt[[ii]][ix.var];
              }

            if((is.na(alt.gold[ix.var]) | alt.gold[ix.var]=='0'| alt.gold[ix.var]=='.')  & (!is.na(raw.data$alt[[ii]][ix.var]) & raw.data$alt[[ii]][ix.var]!="0" & raw.data$alt[[ii]][ix.var]!="." ))
              {
                alt.gold[ix.var] <- raw.data$alt[[ii]][ix.var];
              }
            if(is.na(raw.data$ref[[ii]][ix.var]))
              {
                log.mat[ix.var,ii] <- paste(log.mat[ix.var,ii],"O",sep="",collapse="");
                direction.by.study.var[ii] <- "X";
              }
            if(!is.na(ref.gold[ix.var]) & !is.na(raw.data$ref[[ii]][ix.var]) )
              {
                match.ref.alt <- (ref.gold[ix.var]==(raw.data$ref[[ii]][ix.var]) & alt.gold[ix.var]==(raw.data$alt[[ii]][ix.var]));                
                flip.ref.alt <- (alt.gold[ix.var]==(raw.data$ref[[ii]][ix.var]) & ref.gold[ix.var]==(raw.data$alt[[ii]][ix.var]));
                mono <- ((raw.data$af[[ii]][ix.var]==0) | (raw.data$af[[ii]][ix.var]==1) | (raw.data$ref[[ii]][ix.var]==".") | (raw.data$ref[[ii]][ix.var]==0) | (raw.data$alt[[ii]][ix.var]==".") | (raw.data$alt[[ii]][ix.var]==0));
                if(is.na(match.ref.alt)) match.ref.alt <- 0;
                if(is.na(flip.ref.alt)) flip.ref.alt <- 0;
                if(is.na(mono)) mono <- 0;
                
                if(match.ref.alt)
                  {
                    U.stat <- U.stat+rm.na(raw.data$ustat[[ii]][ix.var]);
                    V.stat.sq <- V.stat.sq+(rm.na(raw.data$vstat[[ii]][ix.var]))^2;
                    nref.var <- nref.var+rm.na(raw.data$nref[[ii]][ix.var]);
                    nalt.var <- nalt.var+rm.na(raw.data$nalt[[ii]][ix.var]);
                    nhet.var <- nhet.var+rm.na(raw.data$nhet[[ii]][ix.var]);
                    no.sample.mat[ii,ix.var] <- raw.data$nSample[[ii]][ix.var];                   
                    log.mat[ix.var,ii] <- paste(log.mat[ix.var,ii],"M",sep="",collapse="");
                  }
                if(flip.ref.alt)
                  {
                    raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                    raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                    U.stat <- U.stat+rm.na(raw.data$ustat[[ii]][ix.var]);
                    nref.var <- nref.var+rm.na(raw.data$nalt[[ii]][ix.var]);
                    nalt.var <- nalt.var+rm.na(raw.data$nref[[ii]][ix.var]);
                    nhet.var <- nhet.var+rm.na(raw.data$nhet[[ii]][ix.var]);
                    V.stat.sq <- V.stat.sq+(rm.na(raw.data$vstat[[ii]][ix.var]))^2;                    
                    log.mat[ix.var,ii] <- paste(log.mat[ix.var,ii],"F",sep="",collapse="");
                    no.sample.mat[ii,ix.var] <- raw.data$nSample[[ii]][ix.var];
                   
                  }
                if(mono)
                  {
                    nref.var <- nref.var+raw.data$nref[[ii]][ix.var];
                    nalt.var <- nalt.var+raw.data$nalt[[ii]][ix.var];
                    nhet.var <- nhet.var+raw.data$nhet[[ii]][ix.var];
                    log.mat[ix.var,ii] <- paste(log.mat[ix.var,ii],"0",sep="",collapse="");
                    no.sample.mat[ii,ix.var] <- raw.data$nSample[[ii]][ix.var];
                   
                  }
                if(!match.ref.alt & !flip.ref.alt & !mono)
                  {
                    raw.data$nSample[[ii]][ix.var] <- NA;
                    raw.data$af[[ii]][ix.var] <- NA;
                    raw.data$ac[[ii]][ix.var] <- NA;
                    raw.data$nref[[ii]][ix.var] <- NA;
                    raw.data$nhet[[ii]][ix.var] <- NA;
                    raw.data$nalt[[ii]][ix.var] <- NA;
                    raw.data$ustat[[ii]][ix.var] <- NA;
                    raw.data$vstat[[ii]][ix.var] <- NA;
                    log.mat[ix.var,ii] <- paste(log.mat[ix.var,ii],"O",sep="",collapse="");
                  }
                if(is.na(raw.data$ustat[[ii]][ix.var]))
                  {
                    direction.by.study.var[ii] <- "X";
                  }
                if(!is.na(raw.data$ustat[[ii]][ix.var]))
                  {
                    if(raw.data$ustat[[ii]][ix.var]>0) direction.by.study.var[ii] <- "+";
                    if(raw.data$ustat[[ii]][ix.var]<0) direction.by.study.var[ii] <- "-";
                    if(raw.data$ustat[[ii]][ix.var]==0) direction.by.study.var[ii] <- "=";
                  }
              }
            maf.pop[ii] <- rm.na((raw.data$af[[ii]])[ix.var]);
          }
        maf.vec[ix.var] <- Sum(maf.pop*no.sample.mat[,ix.var])/Sum(no.sample.mat[,ix.var]);
        direction.by.study[ix.var] <- paste(direction.by.study.var,sep='',collapse='');
        QC.by.study[ix.var] <- paste(log.mat[ix.var,],sep=",",collapse=",");
        beta1.est[ix.var] <- U.stat/V.stat.sq;
        beta1.sd[ix.var] <- sqrt(1/V.stat.sq);
        hsq.est[ix.var] <- (beta1.est[ix.var])*(beta1.est[ix.var])*V.stat.sq/Sum(no.sample.mat[,ix.var]);
        no.sample.var[ix.var] <- Sum(no.sample.mat[,ix.var]);
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
    anno.out <- raw.data$anno;
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
                                 afCase=list(rep(NA,length(ref.out[[1]]))),
                                 anno=anno.out);
    
    
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
                QC.by.study=QC.by.study,
                no.sample=no.sample.var,
                beta1.est=beta1.est,
                beta1.sd=beta1.sd,
                hsq.est=hsq.est,
                nearby=nearby,
                pos=raw.data$pos));
  }
