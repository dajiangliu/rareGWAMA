rvmeta.inVAR.single <- function(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot=0,alpha=0.05,weight=c('MB','MZ'))
  {
    if(length(alternative)>1) alternative <- "two.sided";
    if(length(weight)>1) weight <- "MB";
    res.list <- list();
    maf.vec <- rep(0,length(maf.vec.list[[1]]));
    score.stat.numer <- maf.vec;
    score.stat.denom.sq <- maf.vec;

    for(ii in 1:length(X.T.times.X.list))
      {
        maf.vec <- maf.vec.list[[ii]]*N.list[[ii]]+maf.vec;
      }
    N <- sum(unlist(N.list));
    maf.vec <- maf.vec/N;
    beta.est.vec <- 0;
    var.est.vec <- 0;
    w.vec <- 0; const.vec <- 0;
    wss.stat <- 0;
    beta.est.meta <- 0;
    var.est.meta <- 0;
    ix.Inf <- 0;
    for(ii in 1:length(X.T.times.X.list))
    {
      if(length(X.T.times.Y.centered.list[[1]])>1)
        {
          beta.est.vec <- X.T.times.Y.centered.list[[ii]]/diag(X.T.times.X.list[[ii]]);
          var.est.vec <- var.Y.list[[ii]]/diag(X.T.times.X.list[[ii]]);
          ix.Inf <- which(diag(X.T.times.X.list[[ii]])==0);
        }
      if(length(X.T.times.Y.centered.list[[1]])==1)
        {
          beta.est.vec <- X.T.times.Y.centered.list[[ii]]/as.numeric(X.T.times.X.list[[ii]]);
          var.est.vec <- var.Y.list[[ii]]/as.numeric(X.T.times.X.list[[ii]]);
          ix.Inf <- which(as.numeric(X.T.times.X.list[[ii]])==0);
        }

      beta.est.vec <- rm.na(beta.est.vec);
      var.est.vec <- rm.na(var.est.vec);
      w.vec <- 1/var.est.vec;
      if(length(ix.Inf)>0)
        {
          w.vec[ix.Inf] <- 0;
          beta.est.vec[ix.Inf] <- 0;
          var.est.vec[ix.Inf] <- 0;
        }
      const.vec <- const.vec+w.vec;
      beta.est.meta <- beta.est.meta+beta.est.vec*w.vec;
      var.est.meta <- var.est.meta+var.est.vec*w.vec^2
    }
    wss.stat <- beta.est.meta/sqrt(var.est.meta);
    beta.est.meta <- beta.est.meta/const.vec;
    beta1.sd.vec <- sqrt(var.est.meta)/const.vec
    if(alternative=='two.sided')
    {
        statistic.vec <- wss.stat^2;
        p.value.vec <- pchisq(statistic.vec,df=1,lower.tail=FALSE);
    }
    if(alternative=='greater')
    {
        statistic.vec <- wss.stat;
        p.value.vec <- pnorm((-1)*statistic.vec);
    }
    if(alternative=='less')
    {
        statistic.vec <- wss.stat;
        p.value.vec <- pnorm(statistic.vec);
    }

    return(list(p.value=p.value.vec,
                statistic=statistic.vec,
                direction=as.numeric(wss.stat>0)*2-1,
                maf.vec=maf.vec,
                beta1.est.vec=beta.est.meta,
                beta1.sd.vec=beta1.sd.vec));
}
