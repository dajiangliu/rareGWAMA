meta.p <- function(p.value.vec,N.vec,alternative=c('two.sided','greater','less'))
  {
    if(length(alternative)>1)
      alternative <- 'two.sided';
    
    z.score.vec <- qnorm(1-p.value.vec);
    statistic <- sum(sqrt(N.vec)*z.score.vec)/sqrt(sum(N.vec));
    if(alternative=='greater')
      p.value <- 1-pnorm(statistic);
    if(alternative=='less')
      p.value <- pnorm(statistic);
    if(alternative=='two.sided')
      {
        p.value <- 1-pchisq(statistic^2,df=1);
        statistic <- statistic^2;
      }
    return(list(p.value=p.value,
                statistic=statistic));
  }

rvmeta.p <- function(score.stat.vec.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot=0,method.name=c('WSS','VT','SKAT','RARECOVER'))
  {
    if(method.name=="WSS") method <- rvmeta.wss;
    if(method.name=="VT") method <- rvmeta.vt;
    if(method.name=="SKAT") method <- rvmeta.skat;
    if(method.name=="RARECOVER") method <- rvmeta.rarecover;

    p.value.vec <- 0;
    N.vec <- unlist(N.list);
    
    for(ii in 1:length(score.stat.vec.list))
      {
        score.stat.vec.list.new <- list();score.stat.vec.list.new[[1]] <- score.stat.vec.list[[ii]];
        maf.vec.list.new <- list();maf.vec.list.new[[1]] <- maf.vec.list[[ii]];
        cov.mat.list.new <- list();cov.mat.list.new[[1]] <- cov.mat.list[[ii]];
        var.Y.list.new <- list();var.Y.list.new[[1]] <- var.Y.list[[ii]];
        N.list.new <- list();N.list.new[[1]] <- N.list[[ii]];

        if(method.name!="SKAT"){
          
          res.ii <- method(score.stat.vec.list.new,
                           maf.vec.list.new,
                           cov.mat.list.new,
                           var.Y.list.new,
                           N.list.new,
                           alternative,
                           no.boot);
          
        }
        if(method.name=="SKAT"){  
          res.ii <- method(score.stat.vec.list.new,
                           maf.vec.list.new,
                           cov.mat.list.new,
                           var.Y.list.new,
                           'beta',
                           N.list.new,
                           alternative,
                           no.boot);         
        }
        
        p.value.vec[ii] <- res.ii$p.value;
        
      }
    res.meta.p <- meta.p(p.value.vec,N.vec,'greater');

    return(list(p.value=res.meta.p$p.value));
  }

