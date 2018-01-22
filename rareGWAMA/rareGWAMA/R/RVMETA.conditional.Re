RVMETA.conditional <- function(scenario=c('gene','region'),score.stat.file,cov.file,gene,test=c('GRANVIL','WSS','SKAT-O','SKAT','VT'),maf.cutoff,no.boot=10000,alternative=c('two.sided','greater','less'),alpha=0.05,extra.pars=list(option=character(0),ix.X1=integer(0),ix.X2=integer(0)))
  {
      capture.output(raw.data.all <- rvmeta.readData( score.stat.file, cov.file, gene));
      if(length(raw.data.all)==0)
          return(list(p.value=NA,
                      statistic=NA,
                      no.var=NA,
                      no.sample=NA));
      res <- list();
      gene.name <- names(raw.data.all);
      p.value <- double(0);
      option <- extra.pars$option;
      statistic <- double(0);pos <- integer(0);anno <- character(0);direction <- integer(0);res.maf.vec <- double(0);beta1.est.vec <- double(0);beta1.sd.vec <- double(0);
      for(kk in 1:length(raw.data.all))
      {
          raw.data <- raw.data.all[[kk]];
          if(scenario=='gene')
          {
              ix.var <- integer(0);
              if(length(raw.data$stat[[1]])>0)
              {
                  ix.var <- 1:length(raw.data$stat[[1]])
              }
          }
          if(scenario=='NS')
          {
              ix.var <- integer(0);
              ix.var <- c(ix.var,grep('Nonsynonymous',raw.data$anno));
              ix.var <- sort(unique(ix.var));
          }
          if(scenario=='LOF')
          {
              ix.var <- integer(0);
              ix.var <- c(ix.var,grep('Stop_Gain',raw.data$anno));
              ix.var <- sort(unique(ix.var));
          }
          if(scenario=='NS-LOF')
          {
              ix.var <- integer(0);
              ix.var <- c(ix.var,grep('Stop_Gain|Nonsynonymous',raw.data$anno));
              ix.var <- sort(unique(ix.var));
          }
          if(scenario=='region')
          {
              raw.data <- NULL;
              stop("region based meta analysis has not been implemented in this software");
          }
          if(length(ix.var)<2 & test!='SINGLE') {
            res[[kk]] <- list(gene.name=gene.name[kk],
                              p.value=NA,
                              statistic=NA,
                              no.var=NA,
                              no.sample=NA,
                              anno=NA,
                              pos=NA);
          }
          if(length(ix.var)>=2) {
              ix.pop <- 1:length(raw.data$nSample);
              score.stat.vec.list <- list();maf.vec.list <- list();cov.mat.list <- list();var.Y.list <- list();N.list <- list();mean.Y.list <- list();pos.list <- list();anno.list <- list();
              score.stat.vec.list.all <- list();maf.vec.list.all <- list();cov.mat.list.all <- list();var.Y.list <- list();N.list <- list();mean.Y.list <- list();pos.list.all <- list();anno.list.all <- list();

              for(ii in 1:length(ix.pop))
              {
                  score.stat.vec.list.all[[ii]] <- rm.na((sqrt((raw.data$stat[[ii]]))*((raw.data$direction[[ii]]))));
                  N.list[[ii]] <- as.integer(mean(raw.data$nSample[[ii]],na.rm=TRUE));
                  cov.mat.list.all[[ii]] <- rm.na(as.matrix((raw.data$cov[[ii]])));
                  ix.0 <- which(raw.data$maf[[ii]]==0);
                  if(length(ix.0)>0)
                  {
                      raw.data$maf[[ii]][ix.0] <- 1;
                  }
                  var.Y.list[[ii]] <- 1;
                  mean.Y.list[[ii]] <- 0;
                  pos.list.all[[ii]] <- (raw.data$pos);
                  anno.list.all[[ii]] <- (raw.data$anno);
              }
              maf.vec.all <- rep(0,length(maf.vec.list.all[[1]]));
              for(ii in 1:length(ix.pop))
              {
                  maf.vec.all <- maf.vec.all+(maf.vec.list.all[[ii]])*(2*N.list[[ii]]);
              }
              maf.vec.all <- maf.vec.all/sum(2*unlist(N.list));
              if(option=='TOP')
              {
                  res.single.all <- rvmeta.CMH(score.stat.vec.list.all,
                                               maf.vec.list.all,
                                               cov.mat.list.all,
                                               var.Y.list,
                                               N.list,
                                               alternative,
                                               no.boot,
                                               alpha,
                                               rv.test='SINGLE',
                                               extra.pars=list());
                  ix.top.all <- which.min(res.single.all$p.value);
                  ix.rare <- which(maf.vec.all<maf.cutoff);
                  ix.var.rare <- as.integer(set.intersect(ix.var,ix.rare));
                  if(length(which(ix.var.rare==ix.top.all))>0)
                  {
                      ix.tmp <- which(ix.var.rare==ix.top.all);
                      ix.cond <- c(ix.var.rare[-ix.tmp],ix.var.rare[ix.tmp]);
                  }
                  if(length(which(ix.var.rare==ix.top.all))==0)
                  {
                      ix.cond <- c(ix.var.rare,ix.top.all);
                  }
                  score.stat.vec.list.cond <- list();
                  cov.mat.list.cond <- list();
                  maf.vec.list.cond <- list();
                  anno.list.cond <- list();
                  pos.list.cond <- list();
                  for(ii in 1:length(ix.pop))
                  {
                    score.stat.vec.list.cond[[ii]] <- score.stat.vec.list.all[[ii]][ix.cond];
                    cov.mat.list.cond[[ii]] <- cov.mat.list.all[[ii]][ix.cond,ix.cond];
                    maf.vec.list.cond[[ii]] <- maf.vec.list.all[[ii]][ix.cond];
                    pos.list.cond[[ii]] <- pos.list.all[[ii]][ix.cond];
                    anno.list.cond[[ii]] <- anno.list.all[[ii]][ix.cond];
                  }
                  maf.top <- maf.vec.all[ix.top.all];
                  maf.vec.X2 <- maf.top;
                  p.value.X2 <- (res.single.all$p.value)[ix.top.all];
                  ix.X1 <- 1:(length(ix.cond)-1);
                  ix.X2 <- length(ix.cond);
              }
              if(length(ix.rare)>1)
              {
                  if(test=='WSS')
                  {
                    res.kk <- (c(cond.rvmeta(score.stat.vec.list.cond,
                                             maf.vec.list.cond,
                                             cov.mat.list.cond,
                                             mean.Y.list,
                                             var.Y.list,
                                             N.list,
                                             alternative,
                                             no.boot,alpha,rv.test='WSS',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2))));
                  }
                  if(test=='GRANVIL')
                    {
                      res.kk <- (c(cond.rvmeta(score.stat.vec.list.cond,
                                               maf.vec.list.cond,
                                               cov.mat.list.cond,
                                               mean.Y.list,
                                               var.Y.list,
                                               N.list,
                                               alternative,
                                               no.boot,alpha,rv.test='GRANVIL',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2))));
                    }
                  if(test=='SKAT')
                  {
                    res.kk <- (c(cond.rvmeta(score.stat.vec.list.cond,
                                             maf.vec.list.cond,
                                             cov.mat.list.cond,
                                             mean.Y.list,
                                             var.Y.list,
                                             N.list,
                                             alternative,
                                             no.boot,alpha,rv.test='SKAT',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2))));
                  }
                  if(test=='VT')
                    {
                      res.kk <- (c(cond.rvmeta(score.stat.vec.list.cond,
                                               maf.vec.list.cond,
                                               cov.mat.list.cond,
                                               mean.Y.list,
                                               var.Y.list,
                                               N.list,
                                               alternative,
                                               no.boot,alpha,rv.test='VT',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2))));
                  }
                  res[[kk]] <- c(res.kk,
                                 list(anno.X1=anno.list.cond[[1]][ix.X1],
                                      anno.X2=anno.list.cond[[1]][ix.X2],
                                      pos.X1=pos.list.cond[[1]][ix.X1],
                                      pos.X2=pos.list.cond[[1]][ix.X2],
                                      p.value.X2=p.value.X2,
                                      maf.vec.X2=maf.vec.X2,
                                      gene.name=gene.name[kk]));
              }
          }
      }
      return(res);
  }
