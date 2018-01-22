rvmeta.CMH.single <- function(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot=0,alpha=0.05,weight=c('MB','MZ'))
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
    
    for(ii in 1:length(X.T.times.X.list))
      {
        score.stat.numer <- score.stat.numer+X.T.times.Y.centered.list[[ii]];
        if(length(score.stat.numer)>1)
            score.stat.denom.sq <- score.stat.denom.sq+as.numeric((N.list[[ii]])*diag(cov.mat.list[[ii]])*(var.Y.list[[ii]]));
        if(length(score.stat.numer)==1)
            score.stat.denom.sq <- score.stat.denom.sq+as.numeric((N.list[[ii]])*as.numeric(cov.mat.list[[ii]])*(var.Y.list[[ii]]));
      }
    wss.stat <- score.stat.numer/sqrt(score.stat.denom.sq);
    if(alternative=='two.sided')
    {
        statistic.vec <- wss.stat^2;
        p.value.vec <- 2*pnorm(-sqrt(statistic.vec));
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
                maf.vec=maf.vec));
  }
