rvmeta.approx.mega <- function(score.stat.vec.list,maf.vec.list,cov.mat.list,mean.Y.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot,alpha=0.05,rv.test,extra.pars=list())
  {
    if(length(alternative)>1) alternative <- "two.sided";
    res.list <- list();
    X.T.times.Y.centered <- rep(0,length(maf.vec.list[[1]]));
    X.T.times.X <- 0*diag(maf.vec.list[[1]]);
    maf.vec <- rep(0,length(score.stat.vec.list[[1]]));
    if(length(maf.vec)<=1 & rv.test!='SINGLE')
      {
        return(list(p.value=NA,
                    statistic=NA,
                    X.T.times.X=X.T.times.X,
                    no.site=0,
                    beta1.est=NA,
                    beta1.sd=NA,
                    hsq.est=NA,
                    beta1.conf.lower=NA,
                    beta1.conf.upper=NA));
      }
    direction.by.study.matrix <- matrix(rep("+",length(score.stat.vec.list[[1]])*length(score.stat.vec.list)),nrow=length(score.stat.vec.list));
    direction.by.study <- '';    
    for(ii in 1:length(score.stat.vec.list))
      {
        score.stat.vec.list[[ii]] <- rm.na(score.stat.vec.list[[ii]]);        
        maf.vec.list[[ii]] <- rm.na(maf.vec.list[[ii]]);       
        cov.mat.list[[ii]] <- rm.na(cov.mat.list[[ii]]);
        X.T.times.Y.centered <- X.T.times.Y.centered+(sqrt(N.list[[ii]]))*
          (score.stat.vec.list[[ii]])*
            sqrt(diag(cov.mat.list[[ii]]))*
              sqrt(var.Y.list[[ii]]);
        ix.minus <- which(score.stat.vec.list[[ii]]<0);
        if(length(ix.minus)>0)
          direction.by.study.matrix[ii,ix.minus] <- "-";
                
        X.T.times.Y.centered <- X.T.times.Y.centered+2*(maf.vec.list[[ii]])*(mean.Y.list[[ii]])*(N.list[[ii]]);        
        X.T.times.X <- X.T.times.X+N.list[[ii]]*(cov.mat.list[[ii]]+4*((maf.vec.list[[ii]])%*%t(maf.vec.list[[ii]])));
        ix.q <- which(diag(cov.mat.list[[ii]])==0);
        if(length(ix.q)>0) direction.by.study.matrix[ii,ix.q] <- "?";
        maf.vec <- maf.vec.list[[ii]]*N.list[[ii]]+maf.vec;
      }
    direction.by.study <- apply(direction.by.study.matrix,2,paste,collapse='',sep='');
    N.vec <- unlist(N.list);
    N <- sum(N.vec);
    maf.vec <- maf.vec/N;
    var.Y.vec <- unlist(var.Y.list);
    mean.Y.vec <- unlist(mean.Y.list);
    mean.Y.pool <- sum(mean.Y.vec*N.vec)/N;
    var.Y.pool <- sum((var.Y.vec+mean.Y.vec^2)*N.vec)/N-(mean.Y.pool)^2;
    X.T.times.Y.centered <- X.T.times.Y.centered-2*maf.vec*mean.Y.pool*N;
    N.list <- list();X.T.times.Y.centered.list <- list();X.T.times.X.list <- list();maf.vec.list <- list();var.Y.list <- list();cov.mat.list <- list();
    N.list[[1]] <- N;
    cov.mat <- (1/N)*X.T.times.X-4*maf.vec%*%(t(maf.vec));
    if(length(maf.vec)<=1)
      {
        return(list(p.value=NA,
                    statistic=NA,
                    X.T.times.X=X.T.times.X,
                    no.site=0,
                    beta1.est=NA,
                    beta1.sd=NA,
                    hsq.est=NA,
                    beta1.conf.lower=NA,
                    beta1.conf.upper=NA));
      }
    
    X.T.times.Y.centered.list[[1]] <- X.T.times.Y.centered;
    X.T.times.X.list[[1]] <- X.T.times.X;
    maf.vec.list[[1]] <- maf.vec;
    var.Y.list[[1]] <- var.Y.pool;
    cov.mat.list[[1]] <- cov.mat;
    ix.rm <- integer(0);
    res <- list();
    if(rv.test=='WSS')
      {
        weight <- extra.pars$weight;
        if(length(weight)!=1) weight <- 'MB';
        res <- rvmeta.CMH.wss(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,weight);
      }
    if(rv.test=='VT')
      {
        res <- rvmeta.CMH.vt(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha);
      }
    if(rv.test=='SKAT')
      {
        kernel <- extra.pars$kernel;
        if(length(kernel)!=1) kernel <- "beta";
        res <- rvmeta.CMH.skat(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,kernel);
      }
    if(rv.test=='STEPUP')
      {
        Q <- 0;
        res <- rvmeta.CMH.stepup(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha);
      }
    if(rv.test=='SINGLE')
      {
        res <- rvmeta.CMH.single(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha);  
      }
    
    ix.var <- 1:length(maf.vec);
    if(rv.test=='VT')
      {
        ix.var <- which(maf.vec<=res$maf.cutoff);
      }
    X.T.times.X.centered <- as.matrix((N.list[[1]]*cov.mat.list[[1]])[ix.var,ix.var]);

    w <- rep(1,length(X.T.times.Y.centered));

    if(rv.test=='WSS')
      {
        if(extra.pars$weight=='MB')
          {
            q <- ((maf.vec)+1/(2*N))*(2*N)/(2*N+2);
            w <- 1/sqrt(N*q*(1-q));
          }
      }
    if(rv.test!='SINGLE')
      {
        beta1.est <- sum((w[ix.var])*(X.T.times.Y.centered[ix.var]))/as.numeric(t(w[ix.var])%*%X.T.times.X.centered%*%(w[ix.var]));
        beta1.sd <- sqrt(1/as.numeric(t(w[ix.var])%*%X.T.times.X.centered%*%(w[ix.var]))*var.Y.pool);
        macf.vec <- 2*(maf.vec[ix.var])*(1-maf.vec[ix.var])+(maf.vec[ix.var])*(maf.vec[ix.var]);
        hsq.est <- beta1.est*beta1.est*(1-prod(1-macf.vec))*prod(1-macf.vec);
        beta1.conf.lower <- beta1.est-1.96*beta1.sd;
        beta1.conf.upper <- beta1.est+1.96*beta1.sd;
      }
    if(rv.test=='SINGLE')
      {
        beta1.est <- X.T.times.Y.centered/diag(X.T.times.X.centered);
        beta1.sd <- sqrt(1/as.numeric(diag(X.T.times.X.centered)*var.Y.pool));
        beta1.conf.lower <- beta1.est-1.96*beta1.sd;
        beta1.conf.upper <- beta1.est+1.96*beta1.sd;
        hsq.est <- beta1.est^2*diag(X.T.times.X.centered);
      }
    res.return <- c(res,                    
                    list(X.T.times.X=X.T.times.X,
                         direction.by.study=direction.by.study,
                         ix.rm=ix.rm,
                         no.site=length(X.T.times.Y.centered.list[[1]]),
                         beta1.est=beta1.est,
                         beta1.sd=beta1.sd,
                         hsq.est=hsq.est,
                         beta1.conf.lower=beta1.conf.lower,
                         beta1.conf.upper=beta1.conf.upper));
    return(res.return);


}
