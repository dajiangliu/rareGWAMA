min.p <- function(score.stat.vec.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot,alpha=0.05,rv.test,extra.pars=list())
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
        ##mac.vec.list.new[[1]] <- mac.vec;
        maf.vec <- maf.vec.list.new[[1]];
        res[[ii]] <- rvmeta.CMH(score.stat.vec.list.new,maf.vec.list.new,cov.mat.list.new,var.Y.list.new,N.list.new,alternative,no.boot,alpha,rv.test,c(extra.pars,list(mac.vec.list=mac.vec.list.new,maf.vec=maf.vec,mac.vec=mac.vec,af.vec.list=maf.vec.list.new,ac.vec.list=mac.vec.list.new,max.TH=extra.pars$max.TH)));
        p.value.vec[ii] <- res[[ii]]$p.value;
      }
    statistic <- min(p.value.vec);
    p.value <- 1-(1-statistic)^length(p.value.vec);
    return(list(p.value=p.value,
                statistic=statistic));
    
  }
