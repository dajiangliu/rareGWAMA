fisher.p <- function(score.stat.vec.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot,alpha=0.05,rv.test,extra.pars=list(weight.by.N="N"))
  {
    res <- list();
    p.value.vec <- 0;
    statistic.vec <- 0;
    for(ii in 1:length(score.stat.vec.list))
      {
        score.stat.vec.list.new <- list();maf.vec.list.new <- list();cov.mat.list.new <- list();var.Y.list.new <- list();N.list.new <- list();mac.vec.list.new <- list();
        maf.vec.list.new[[1]] <- maf.vec.list[[ii]];
        score.stat.vec.list.new[[1]] <- score.stat.vec.list[[ii]];
        maf.vec.list.new[[1]] <- maf.vec.list[[ii]];
        cov.mat.list.new[[1]] <- cov.mat.list[[ii]];
        var.Y.list.new[[1]] <- var.Y.list[[ii]];
        N.list.new[[1]] <- N.list[[ii]];
        mac.vec.list.new[[1]] <- (maf.vec.list[[ii]])*(N.list[[ii]])*2;
        mac.vec <- mac.vec.list.new[[1]];
        maf.vec <- maf.vec.list.new[[1]];
        res[[ii]] <- rvmeta.CMH(score.stat.vec.list.new,maf.vec.list.new,cov.mat.list.new,var.Y.list.new,N.list.new,alternative,no.boot,alpha,rv.test,c(extra.pars,list(mac.vec.list=mac.vec.list.new,maf.vec=maf.vec,mac.vec=mac.vec,af.vec.list=maf.vec.list.new,ac.vec.list=mac.vec.list.new,max.TH=extra.pars$max.TH)));
        p.value.vec[ii] <- res[[ii]]$p.value;
      }
    if(length(extra.pars$weight.by.N)==0) extra.pars$weight.by.N <- "N";
    if(extra.pars$weight.by.N=="N")
      {
        statistic <- sum(-2*log(p.value.vec),na.rm=TRUE);
        p.value <- pchisq(statistic,df=2*length(p.value.vec),lower.tail=FALSE);
      }
    if(extra.pars$weight.by.N=="Y")
      {
        N.vec <- unlist(N.list);
        statistic <- sum(-2*N.vec*log(p.value.vec),na.rm=TRUE);
        p.value <- try(davies(statistic,lambda=rep(N.vec,2))$Qq,silent=TRUE);
        if(p.value==0)
          {
            p.value <- try(liu(statistic,lambda=rep(N.vec,2)),silent=TRUE);
          }
      }    
    return(list(p.value=p.value,
                statistic=statistic));    
  }
