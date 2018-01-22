pvt.core <- function(x,mu,sigma,alternative=c('two.sided','greater'))
  {
    sum.all <- 0;
    if(alternative=='greater')
      {
        dimn <- length(mu);
        for(ii in 1:dimn)
          {
            mat.ii <- combn(1:dimn,ii);
            sum.ii <- 0;
            for(jj in 1:ncol(mat.ii))
              {
                sum.ii <- sum.ii+as.numeric(pmvnorm(upper=rep(-x,ii),lower=rep(-Inf,ii),mean=mu[mat.ii[,jj]],sigma=sigma.jj));
              }
            if(ii%%2==1)
              sum.all <- sum.all+sum.ii;
            if(ii%%2==0)
              sum.all <- sum.all-sum.ii;
          }
        return(sum.all);
      }
    
    if(alternative=='two.sided')
      {
        dimn <- length(mu);
        for(ii in 1:dimn)
          {
            mat.ii <- combn(1:dimn,ii);
            sum.ii <- 0;
            for(jj in 1:ncol(mat.ii))
              {
                
                sum.ii <- sum.ii+2*as.numeric(pmvnorm(upper=rep(-sqrt(x),ii),lower=rep(-Inf,ii),mean=mu[mat.ii[,jj]],sigma=sigma.jj));
            }
            if(ii%%2==1)
              sum.all <- sum.all+sum.ii;
            if(ii%%2==0)
              sum.all <- sum.all-sum.ii;
          }
        return(sum.all);
      }
    
  }
pvt <- function(statistic,mu,sigma,alternative=c('greater','less','two.sided'))
  {
    p.value <- NA;
    if(length(alternative)>1) alternative <- "two.sided";
    if(alternative=='two.sided')
      {
        p.tmp <- 1-pmaxnormsq(statistic,mu,sigma);
        if(!is.na(p.tmp)){
          if(p.tmp>0) p.value <- p.tmp;
          
          if(p.tmp==0)
            {              
              p.value <- pvt.core(statistic,mu,sigma,'two.sided')
            }
        }
      }
    if(alternative=='less')
      {
        p.tmp <- 1-pmaxnormsq(-statistic,mu,sigma);
        if(p.tmp>0) p.value <- p.tmp;
        if(p.tmp==0)
          {
            p.value <- pvt.core(statistic,mu,sigma,'greater')
          }        
      }
    if(alternative=='greater')
      {
        p.tmp <- 1-pmaxnorm(statistic,mu,sigma);
        if(p.tmp>0) p.value <- p.tmp;
        if(p.tmp==0)
          {
            p.value <- pvt.core(statistic,mu,sigma,'greater')
          }
        
      }
    return(p.value);
  }
