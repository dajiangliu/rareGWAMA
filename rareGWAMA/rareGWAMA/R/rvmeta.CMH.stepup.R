#rvmeta.CMH.stepup(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot);
rvmeta.CMH.stepup <- function(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha=0.05)
  {
    if(length(alternative)>1) alternative <- "two.sided";
    X.T.times.Y <- rep(0,length(maf.vec.list[[1]]));
    X.T.times.X <- 0*diag(length(maf.vec.list[[1]]));
    maf.vec <- rep(0,length(maf.vec.list[[1]]));
    for(ii in 1:length(maf.vec.list))
      {     
        X.T.times.Y <- X.T.times.Y+X.T.times.Y.centered.list[[ii]];
        X.T.times.X <- X.T.times.X+(N.list[[ii]])*(cov.mat.list[[ii]])*(var.Y.list[[ii]]);
        maf.vec <- maf.vec+(maf.vec.list[[ii]])*(N.list[[ii]]);
      }
    N <- sum(unlist(N.list));
    maf.vec <- maf.vec/N;
    cov.mat <- X.T.times.X;
    var.Y <- unlist(var.Y.list);
    err.msg <- vector(length=0);

    if(no.boot==0) {
      err.msg <- c(err.msg,"no.boot reset to 1000");
      no.boot <- 1000;
    }
    #t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    if(alternative=='two.sided'){
      alternativePt <- 0
    }
    if(alternative=='greater'){
      alternativePt <- 1;
    }
    if(alternative=='less'){
      alternativePt <- -1;
    }
    res.data <- .C("genericStat",
                   X_T_times_Y_permPt=as.double(X.T.times.Y),
                   X_T_times_XPt=as.double(X.T.times.X),
                   noMarkerPt=as.integer(nrow(cov.mat)),
                   NPt=as.integer(N),
                   noBootPt=as.integer(1),
                   varYPt=as.double(1),
                   statisticPt=as.double(rep(0,1)),
                   alternativePt=as.integer(alternativePt),
                   testPt=as.integer(1),
                   extraPar=as.double(0));
    statData <- (res.data$statisticPt)[1];
    
    stat.perm <- double(0);
    if(no.boot/1000!=as.integer(no.boot/1000)) no.boot <- 1000*as.integer(no.boot/1000)+1000;
    res.svd <- svd(X.T.times.X);
    U <- res.svd$u;
    V <- res.svd$v;
    D <- res.svd$d;
    #X.T.times.X <- X.T.times.X*N;
    res.perm <- .C("genericBoot",
                   U=as.double(U),
                   V=as.double(V),
                   D=as.double(D),
                   noMarkerPt=as.integer(nrow(cov.mat)),
                   X_T_times_XPt=as.double(X.T.times.X),
                   NPt=as.integer(N),
                   noBootPt=as.integer(no.boot),
                   varYPt=as.double(1),
                   statisticPt=as.double(rep(0,1000)),
                   alternativePt=as.integer(alternativePt),
                   testPt=as.integer(1),
                   extraPar=as.double(c((res.data$statisticPt)[1],alpha,0,rep(0,no.boot))));

    p.value <- res.perm$extraPar[3];
    
    return(list(statistic=res.data$statisticPt[1],
                p.value=p.value));
  }
