#' Perform conditional analysis for gene-level tests
#'
#' @param range.name name of the range to be analyzed (for example, it can be a gene name e.g. APOE)
#' @param score.stat.file files of score statistics
#' @param cov.file covariance matrix files
#' @param candidate.variant.vec Vectors of candidate variants: e.g. c("1:123","1:1234"). The quotation is necessary!! 
#' @param known.variant.vec  Vectors of known variants: e.g. c("1:123","1:1234"). The quotation is necessary!! 
#' @param test test of rare variant tests
#' @param maf.cutoff Cutoffs of MAF used for determining rare variants
#' @param alternative Alternative hypothesis to be tested
#' @param ix.gold Index of the gold standard population: Used for flipping alleles
#' @param out.digits The number of digits used in the output
#' @param callrate.cutoff Cutoff of call rates. Sites with callrates lower than the cutoff will be labeled as missing
#' @param hwe.cutoff Cutoff of HWE p-values. Sites with HWE pvalues lower than the cutoff will be labeled as missing
#' @param max.VT The maximum number of thresholds used in VT; Setting max.VT to 10 can improve the speed for calculation without affecting the power too much. The default parameter is NULL, which does not set upper limit on the number of variable frequency threhsold. 
#' @export
conditional.rareMETALS.range <- function(range.name=NULL,score.stat.file,cov.file,candidate.variant.vec,known.variant.vec,test='GRANVIL',maf.cutoff,alternative=c('two.sided','greater','less'),ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,max.VT=NULL)
  {
    ANNO <- "gene";
    no.boot <- 0;
    alpha <- 0.05;
    extra.pars <- list(ix.gold=ix.gold,
                       out.digits=out.digits,
                       QC.par=list(callrate.cutoff=callrate.cutoff,hwe.cutoff=hwe.cutoff));
    tabix.known.variant <- get.tabix.range(known.variant.vec);
    tabix.candidate.variant <- get.tabix.range(candidate.variant.vec);
    tabix.all <- paste(tabix.known.variant,tabix.candidate.variant,collapse=",",sep=",");
    raw.data.all <- list();
    capture.output(raw.data.all[[1]] <- rvmeta.readDataByRange( score.stat.file, cov.file,tabix.all)[[1]]);
    res.null <- list(gene.name=NA,
                     p.value=NA,
                     statistic=NA,
                     no.var=NA,
                     no.sample=NA,
                     anno=NA,
                     ref=NA,
                     alt=NA,
                     maf.cutoff=NA,
                     direction.burden.by.study=NA,
                     direction.meta.single.var=NA,
                     beta1.est=NA,
                     no.site=NA,
                     beta1.sd=NA,
                     hsq.est=NA,
                     pos=NA);
    if(length(raw.data.all)==0)
      {
        return(list(res.null));        
      }
    res <- list(list(p.value=NA,statistic=NA,no.var=NA,no.sample=NA));
    p.value <- double(0);
    ref <- character(0);
    alt <- ref;
    direction.by.study <- character(0);
    ix.gold <- extra.pars$ix.gold;
    if(length(extra.pars$ix.gold)==0) {
      ix.gold <- 1;
    }
    out.digits <- extra.pars$out.digits;
    if(length(out.digits)==0) out.digits <- 4;
    statistic <- double(0);pos <- integer(0);anno <- character(0);direction <- integer(0);res.maf.vec <- double(0);beta1.est.vec <- double(0);beta1.sd.vec <- double(0);
    gene.name.out <- NA;p.value.out <- NA;statistic.out <- NA;no.site.out <- NA;beta1.est.out <- NA;beta1.sd.out <- NA;maf.cutoff.out <- NA;direction.burden.by.study.out <- NA;direction.meta.single.var.out <- NA;pos.ref.alt.out <- NA;pos.ref.alt.known.out <- NA;top.singlevar.pos <- NA;
    top.singlevar.refalt <- NA;top.singlevar.pval <- NA;top.singlevar.af <- NA;
    pos.single.out <- character(0);p.value.single.out <- double(0);ref.single.out <- character(0);alt.single.out <- character(0);anno.single.out <- character(0);maf.single.out <- double(0);beta1.est.single.out <- double(0);beta1.sd.single.out <- double(0);pos.ref.alt.known.single.out <- character(0);direction.single.out <- character(0);
    for(kk in 1:length(raw.data.all))
      {
        raw.data <- raw.data.all[[kk]];
        if(length(extra.pars$QC.par)>0) 
          raw.data <- QC(raw.data,extra.pars$QC.par,cov=1);        
        ix.var <- integer(0);
        if(ANNO=='gene')
          {
            if(length(raw.data$ustat[[1]])>0)
              ix.var <- 1:length(raw.data$ustat[[1]]);
          }
        ix.var <- c(ix.var,match(known.variant.vec,raw.data$pos));
        ix.var <- sort(unique(ix.var));
        ix.tmp <- match(known.variant.vec,(raw.data$pos)[ix.var]);
        if(length(ix.var)==length(ix.tmp))
          {
            res.null$gene.name <- range.name[kk];
            res[[kk]] <- res.null;
          }
        if(length(ix.var)>length(ix.tmp)) {
          ix.pop <- 1:length(raw.data$nSample);
          score.stat.vec.list <- list();mac.vec.list <- list();maf.vec.list <- list();cov.mat.list <- list();var.Y.list <- list();N.list <- list();mean.Y.list <- list();pos.list <- list();anno.list <- list();ac.vec.list <- list();af.vec.list <- list();
          ref.list <- list();
          alt.list <- list();
          no.sample <- 0;
          for(ii in 1:length(ix.pop))
            {
                if(length(raw.data$covXZ[[ii]])>0) {
                    warning(paste0("Study  ",ii," is analyzed as binary trait. It is advised to use rareMETALS2 for meta-analysis"))
                }
                
              N.list[[ii]] <- rm.na(as.integer(mean(raw.data$nSample[[ii]],na.rm=TRUE)));
              no.sample <- no.sample+N.list[[ii]];
              U.stat <- rm.na(raw.data$ustat[[ii]][ix.var]);
              V.stat <- rm.na(raw.data$vstat[[ii]][ix.var]);
              
              score.stat.vec.list[[ii]] <- rm.na(U.stat/V.stat);
              
              cov.mat.list[[ii]] <- as.matrix(rm.na(as.matrix(raw.data$cov[[ii]])[ix.var,ix.var]));
              var.Y.list[[ii]] <- 1;
              mean.Y.list[[ii]] <- 0;
              af.vec.list[[ii]] <- rm.na((raw.data$af[[ii]])[ix.var]);
              ac.vec.list[[ii]] <- rm.na((raw.data$ac[[ii]])[ix.var]);
              pos.list[[ii]] <- (raw.data$pos)[ix.var];
              ref.list[[ii]] <- (raw.data$ref)[[ii]][ix.var];
              alt.list[[ii]] <- (raw.data$alt)[[ii]][ix.var];
              anno.list[[ii]] <- (raw.data$anno)[ix.var];
            }
          if(length(ix.pop)>1)
            {
              for(ii in 1:length(ix.var))
                {
                  for(jj in (1:length(ix.pop))[-ix.gold])
                    {
                      if(is.na(ref.list[[ix.gold]][ii]) & !is.na(ref.list[[jj]][ii]))
                        {
                          ref.list[[ix.gold]][ii] <- ref.list[[jj]][ii];
                          alt.list[[ix.gold]][ii] <- alt.list[[jj]][ii];
                        }
                      if((is.na(alt.list[[ix.gold]][ii]) | (alt.list[[ix.gold]][ii])=='0') & !is.na(alt.list[[jj]][ii]))
                        {
                          alt.list[[ix.gold]][ii] <- alt.list[[jj]][ii];
                        }
                      if(!is.na(ref.list[[jj]][ii]) & !is.na(ref.list[[ix.gold]][ii]))
                        {
                          if(ref.list[[jj]][ii]==alt.list[[ix.gold]][ii] & (ref.list[[jj]][ii])!=(ref.list[[ix.gold]][ii]))
                            {
                              tmp <- ref.list[[jj]][ii];
                              ref.list[[jj]][ii] <- alt.list[[jj]][ii];
                              alt.list[[jj]][ii] <- tmp;                             
                              score.stat.vec.list[[jj]][ii] <- (-1)*(score.stat.vec.list[[jj]][ii]);
                              cov.mat.list[[jj]][ii,] <- (-1)*(cov.mat.list[[jj]][ii,]);
                              cov.mat.list[[jj]][,ii] <- (-1)*(cov.mat.list[[jj]][,ii]);
                              af.vec.list[[jj]][ii] <- 1-af.vec.list[[jj]][ii];
                              ac.vec.list[[jj]][ii] <- 2*N.list[[jj]]-ac.vec.list[[jj]][ii];
                            }
                        }
                    }
                }
            }

          res.check.mono <- check.mono(af.vec.list,ac.vec.list,N.list)
          af.vec.list <- res.check.mono$af.vec.list;
          ac.vec.list <- res.check.mono$ac.vec.list;
          
          maf.vec <- rep(0,length(af.vec.list[[1]]));
          af.vec <- maf.vec;
          for(ii in 1:length(ix.pop))
            {
              af.vec <- af.vec+(af.vec.list[[ii]])*(2*N.list[[ii]]);
              ac.vec <- ac.vec+ac.vec.list[[ii]];
            }
          af.vec <- af.vec/sum(2*unlist(N.list));
          maf.vec <- af.vec;
          mac.vec <- ac.vec;
          ix.major <- which(af.vec>0.5);
          if(length(ix.major)>0)
            {
              maf.vec[ix.major] <- 1-maf.vec[ix.major];
              mac.vec[ix.major] <- no.sample-ac.vec[ix.major];
              tmp.major <- flip.score.cov(score.stat.vec.list,cov.mat.list,ix.major);
              score.stat.vec.list <- tmp.major$score.stat.vec.list;
              cov.mat.list <- tmp.major$cov.mat.list;
            }
          ix.rare <- which(maf.vec<maf.cutoff);
          ix.tmp <- match(known.variant.vec,pos.list[[ix.gold]]);
          ix.rare <- unique(c(ix.rare,ix.tmp));
          ix.tmp <- match(known.variant.vec,pos.list[[ix.gold]][ix.rare]);
          ix.rare <- c(ix.rare[-ix.tmp],ix.rare[ix.tmp]);
          maf.vec.rare <- maf.vec[ix.rare];
          mac.vec.rare <- mac.vec[ix.rare];
          
          if(length(ix.rare)==length(ix.tmp))
            {
              res.null$gene.name <- range.name[kk];
              res[[kk]] <- res.null;
            }
          if(length(ix.rare)>length(ix.tmp))
            {
              for(ii in 1:length(ix.pop))
                {
                  score.stat.vec.list[[ii]] <- score.stat.vec.list[[ii]][ix.rare];
                  cov.mat.list[[ii]] <- as.matrix(cov.mat.list[[ii]][ix.rare,ix.rare]);
                  af.vec.list[[ii]] <- af.vec.list[[ii]][ix.rare];
                  ac.vec.list[[ii]] <- ac.vec.list[[ii]][ix.rare];
                  anno.list[[ii]] <- anno.list[[ix.gold]][ix.rare];
                  pos.list[[ii]] <- pos.list[[ix.gold]][ix.rare];
                  ref.list[[ii]] <- ref.list[[ix.gold]][ix.rare];
                  alt.list[[ii]] <- alt.list[[ix.gold]][ix.rare];
                }
              res.extra <- list(anno=anno.list[[ix.gold]],
                                pos=pos.list[[ix.gold]],
                                ref=ref.list[[ix.gold]],
                                alt=alt.list[[ix.gold]],
                                af.vec.list=af.vec.list,
                                score.stat.vec.list=score.stat.vec.list,
                                cov.mat.list=cov.mat.list,
                                pos.list=pos.list,
                                alt.list=alt.list,
                                ref.list=ref.list,
                                raw.data=raw.data,
                                gene.name=range.name[kk]);
              ix.X1 <- 1:(length(ix.rare)-length(ix.tmp));
              ix.X2 <- (length(ix.rare)-length(ix.tmp)+1):length(ix.rare);
              
              if(test=='WSS')
                {
                  res.kk <- (c(cond.rvmeta(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='WSS',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2,maf.vec=maf.vec.rare,mac.vec=mac.vec.rare))));
                  res[[kk]] <- c(res.kk,res.extra);
                  res[[kk]]$maf.cutoff <- maf.cutoff;
                }
              if(test=='GRANVIL')
                {
                  weight <- "MZ";
                  res.kk <- c(cond.rvmeta(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='WSS',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2,weight=weight,mac.vec.list=ac.vec.list,maf.vec=maf.vec.rare,mac.vec=mac.vec.rare)));
                  res[[kk]] <- c(res.kk,res.extra);
                  res[[kk]]$maf.cutoff <- maf.cutoff;
                }
              if(test=='SKAT')
                {
                  res.kk <- (c(cond.rvmeta(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='SKAT',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2,mac.vec.list=ac.vec.list,maf.vec=maf.vec.rare,mac.vec=mac.vec.rare))));
                  res[[kk]] <- c(res.kk,res.extra);
                  res[[kk]]$maf.cutoff <- maf.cutoff;
                }
              if(test=='VT')
                {
                  res.kk <- (c(cond.rvmeta(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='VT',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2,mac.vec.list=ac.vec.list,maf.vec=maf.vec.rare,mac.vec=mac.vec.rare,max.TH=max.VT))));
                  res[[kk]] <- c(res.kk,res.extra);
                }
              gene.name.out[kk] <- res[[kk]]$gene.name;
              p.value.out[kk] <- format(res[[kk]]$p.value,digits=out.digits);
              statistic.out[kk] <- format(res[[kk]]$statistic,digits=out.digits);
              no.site.out[kk] <- res[[kk]]$no.site;
              beta1.est.out[kk] <- format(res[[kk]]$beta1.est,digits=out.digits);
              beta1.sd.out[kk] <- format(res[[kk]]$beta1.sd,digits=out.digits);
              maf.cutoff.out[kk] <- format(res[[kk]]$maf.cutoff,digits=out.digits);

              direction.burden.by.study.out[kk] <- res[[kk]]$direction.burden.by.study;
              direction.meta.single.var.out[kk] <- res[[kk]]$direction.meta.single.var;          
              pos.ref.alt.out[kk] <- paste(res[[kk]]$pos[ix.X1],res[[kk]]$ref[ix.X1],res[[kk]]$alt[ix.X1],sep='/',collapse=',');
              
              pos.ref.alt.known.out[kk] <- paste(res[[kk]]$pos[ix.X2],res[[kk]]$ref[ix.X2],res[[kk]]$alt[ix.X2],sep='/',collapse=',');
              p.value.single.out <- c(p.value.single.out,format(res[[kk]]$p.value.single,digit=out.digits));

              pos.single.out <- c(pos.single.out,paste('chr',res[[kk]]$pos[ix.X1],sep=''));
              direction.single.out <- c(direction.single.out,res[[kk]]$direction.single.vec);
              ref.single.out <- c(ref.single.out,res[[kk]]$ref[ix.X1]);
              alt.single.out <- c(alt.single.out,res[[kk]]$alt[ix.X1]);
              beta1.est.single.out <- c(beta1.est.single.out,res[[kk]]$beta1.est.single[ix.X1]);
              beta1.sd.single.out <- c(beta1.sd.single.out,res[[kk]]$beta1.sd.single[ix.X1]);
              maf.single.out <- c(maf.single.out,res[[kk]]$maf.vec[ix.X1]);
              anno.single.out <- c(anno.single.out,res[[kk]]$anno[ix.X1]);
              pos.ref.alt.known.single.out <- c(pos.ref.alt.known.single.out,
                                                rep(pos.ref.alt.known.out[kk],length(ix.X1)));
              ix.best <- res[[kk]]$ix.best;
              top.singlevar.pos[kk] <- res[[kk]]$pos[ix.X1][ix.best];
              top.singlevar.refalt[kk] <- paste(c(res[[kk]]$ref[ix.X1][ix.best],res[[kk]]$alt[ix.X1][ix.best]),sep="/",collapse="/");
              top.singlevar.pval[kk] <- format(res[[kk]]$singlevar.pval.vec[ix.best],digits=out.digits);
              top.singlevar.af[kk] <- format(res[[kk]]$singlevar.af.vec[ix.best],digits=out.digits);

            }
        }
    }
    res.out <- cbind(gene.name.out,p.value.out,statistic.out,no.site.out,beta1.est.out,beta1.sd.out,maf.cutoff.out,direction.burden.by.study.out,direction.meta.single.var.out,
                     top.singlevar.pos,top.singlevar.refalt,top.singlevar.pval,top.singlevar.af,
                     pos.ref.alt.out,pos.ref.alt.known.out);
    res.single.out <- cbind(pos.single.out,ref.single.out,alt.single.out,p.value.single.out,maf.single.out,beta1.est.single.out,beta1.sd.single.out,direction.single.out,anno.single.out,pos.ref.alt.known.single.out);
    return(list(res.list=res,
                res.single.out=res.single.out,
                res.out=res.out));
  }
r2cov.mat <- function(r2.mat,maf.vec)
  {
    var.vec <- sqrt(maf.vec*(1-maf.vec)*2);
    var.mat <- (var.vec%*%t(var.vec));
    cov.mat <- r2.mat*var.mat;
    return(cov.mat);
  }
