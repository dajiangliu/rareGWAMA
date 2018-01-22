rareMETALS.gene.pool <- function(ANNO,score.stat.file,cov.file,gene,test='GRANVIL',maf.cutoff,no.boot=10000,alternative=c('two.sided','greater','less'),alpha=0.05,extra.pars=list(),gene.file="refFlat_hg19.txt.gz")
  {
    capture.output(raw.data.all <- rvmeta.readDataByGene( score.stat.file, cov.file, gene.file,gene));
    if(length(raw.data.all)==0)
      {
        return(list(list(p.value=NA,
                         skip=1,
                         statistic=NA,
                         no.var=0,
                         no.site=0,
                         no.sample=NA)));
        
      }
    res <- list(list(p.value=NA,statistic=NA,no.var=NA,no.sample=NA));
    gene.name <- names(raw.data.all);
    p.value <- double(0);
    ref <- character(0);
    alt <- ref;
    direction.by.study <- character(0);
    ix.gold <- extra.pars$ix.gold;
    if(length(extra.pars$ix.gold)==0) {
      ix.gold <- 1;
    }
    statistic <- double(0);pos <- integer(0);anno <- character(0);direction <- integer(0);res.maf.vec <- double(0);beta1.est.vec <- double(0);beta1.sd.vec <- double(0);
    for(kk in 1:length(raw.data.all))
      {
        raw.data <- raw.data.all[[kk]];
        ix.var <- integer(0);
        if(ANNO=='gene') ix.var <- 1:length(raw.data$ustat[[1]]);
        if(ANNO!='gene')
          {
            if(length(raw.data$ustat[[1]])>0)
              {
                ix.var <- c(ix.var,grep(ANNO,raw.data$anno));
              }
          }
        ix.var <- sort(unique(ix.var));
        if(length(ix.var)==0)
          {
            res[[kk]] <- list(gene.name=gene.name[kk],
                              p.value=NA,
                              statistic=NA,
                              no.var=NA,
                              no.sample=NA,
                              anno=NA,
                              ref=NA,
                              alt=NA,
                              beta1.est=NA,
                              no.site=NA,
                              beta1.sd=NA,
                              hsq.est=NA,
                              pos=NA);
          }

        if(length(ix.var)>=1) {
          ix.pop <- 1:length(raw.data$nSample);
          score.stat.vec.list <- list();maf.vec.list <- list();cov.mat.list <- list();var.Y.list <- list();N.list <- list();mean.Y.list <- list();pos.list <- list();anno.list <- list();
          ref.list <- list();
          alt.list <- list();
          for(ii in 1:length(ix.pop))
            {
              N.list[[ii]] <- rm.na(as.integer(mean(raw.data$nSample[[ii]],na.rm=TRUE)));
              U.stat <- rm.na(raw.data$ustat[[ii]][ix.var]);
              V.stat <- rm.na(raw.data$vstat[[ii]][ix.var]);
              score.stat.vec.list[[ii]] <- rm.na(U.stat/V.stat);              
              cov.mat.list[[ii]] <- as.matrix(rm.na(as.matrix(raw.data$cov[[ii]])[ix.var,ix.var]));
              var.Y.list[[ii]] <- 1;
              mean.Y.list[[ii]] <- 0;
              maf.vec.list[[ii]] <- rm.na((raw.data$af[[ii]])[ix.var]);
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
                              maf.vec.list[[jj]][ii] <- 1-maf.vec.list[[jj]][ii];
                            }
                        }
                    }
                }
            }          
          maf.vec <- rep(0,length(maf.vec.list[[1]]));
          for(ii in 1:length(ix.pop))
            {
              ix.1 <- which(maf.vec.list[[ii]]==1);
              if(length(ix.1)>0) maf.vec.list[[ii]][ix.1] <- 0;
              maf.vec <- maf.vec+(maf.vec.list[[ii]])*(2*N.list[[ii]]);
            }
          maf.vec <- maf.vec/sum(2*unlist(N.list));
          ix.rare <- which(maf.vec<maf.cutoff);
          maf.vec.rare <- maf.vec[ix.rare];
          if(length(ix.rare)==0)
            {
              res[[kk]] <- list(
                p.value=NA,
                statistic=NA,
                no.site=0,
                skip=1,
                gene.name=gene.name[[kk]],
                anno=NA,
                beta1.est=NA,
                beta1.sd=NA,
                hsq.est=NA,
                ref=NA,
                alt=NA,
                pos=NA);
            }
          if(length(ix.rare)>=1)
            {
              for(ii in 1:length(ix.pop))
                {
                  score.stat.vec.list[[ii]] <- score.stat.vec.list[[ii]][ix.rare];
                  cov.mat.list[[ii]] <- as.matrix(cov.mat.list[[ii]][ix.rare,ix.rare]);
                  maf.vec.list[[ii]] <- maf.vec.list[[ii]][ix.rare];
                  anno.list[[ii]] <- anno.list[[ix.gold]][ix.rare];
                  pos.list[[ii]] <- pos.list[[ix.gold]][ix.rare];
                  ref.list[[ii]] <- ref.list[[ix.gold]][ix.rare];
                  alt.list[[ii]] <- alt.list[[ix.gold]][ix.rare];
                }
              if(test=='WSS')
                {
                  res.kk <- (c(rvmeta.CMH(score.stat.vec.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='WSS',extra.pars=list(weight='MB'))));
                  res[[kk]] <- c(res.kk,
                                 list(anno=anno.list[[ix.gold]],
                                      pos=pos.list[[ix.gold]],
                                      maf.vec.list=maf.vec.list,
                                      ref=ref.list[[ix.gold]],
                                      alt=alt.list[[ix.gold]],
                                      gene.name=gene.name[kk]));
                }
              if(test=='GRANVIL')
                {
                  res.kk <- (c(rvmeta.CMH(score.stat.vec.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='WSS',extra.pars=list(weight='MZ'))));
                  res[[kk]] <- c(res.kk,
                                 list(anno=anno.list[[ix.gold]],
                                      pos=pos.list[[ix.gold]],
                                      ref=ref.list[[ix.gold]],
                                      alt=alt.list[[ix.gold]],
                                      maf.vec.list=maf.vec.list,
                                      gene.name=gene.name[kk]));
                }
              if(test=='SKAT')
                {
                  res.kk <- (c(rvmeta.CMH(score.stat.vec.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='SKAT',extra.pars=list(kernel='beta'))));
                  res[[kk]] <- c(res.kk,
                                 list(anno=anno.list[[ix.gold]],
                                      pos=pos.list[[ix.gold]],
                                      ref=ref.list[[ix.gold]],
                                      alt=alt.list[[ix.gold]],
                                      maf.vec.list=maf.vec.list,
                                      gene.name=gene.name[kk]));
                }
              if(test=='VT')
                {
                  res.kk <- (c(rvmeta.approx.mega(score.stat.vec.list,maf.vec.list,cov.mat.list,mean.Y.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='VT',extra.pars=list())));
                  res[[kk]] <- c(res.kk,
                                 list(anno=anno.list[[ix.gold]],
                                      pos=pos.list[[ix.gold]],
                                      ref=ref.list[[ix.gold]],
                                      alt=alt.list[[ix.gold]],
                                      maf.vec=maf.vec.rare,
                                      maf.vec.list=maf.vec.list,
                                      gene.name=gene.name[kk]));
                }
            }
        }
    }
    return(res);
  }

r2cov.mat <- function(r2.mat,maf.vec)
  {
    var.vec <- sqrt(maf.vec*(1-maf.vec)*2);
    var.mat <- (var.vec%*%t(var.vec));
    cov.mat <- r2.mat*var.mat;
    return(cov.mat);
  }
