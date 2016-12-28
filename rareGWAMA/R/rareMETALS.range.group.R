#' Meta-analysis of gene-level tests by range;
#'
#' @param score.stat.file files of score statistics
#' @param cov.file Covariance matrix files
#' @param range tabix range for each gene/region; Tabix range needs to be in the format of "1:123-1234". Quotation marks are necessary. 
#' @param range.name The name of the range,e.g. gene names can be used;
#' @param refaltList A list of reference, alternative allele, a vector of alternative allele frequencies; Specifically, the list should consist of 
#' @param test rare variant tests to be used
#' @param maf.cutoff MAF cutoff used to analyze variants; Default value is 1, i.e. no cutoffs are applied. MAFs are based upon the sample MAFs. 
#' @param alternative alternative hypothesis to be specified; Only applicable to GRANVIL test;
#' @param out.digits Number of digits used in the output
#' @param callrate.cutoff Cutoffs of call rate, lower than which will NOT be analyzed (labelled as missing)
#' @param hwe.cutoff Cutoffs of HWE p-values
#' @param max.VT The maximum number of thresholds used in VT; For small p-values, the calculation of VT p-values can be very slow. Setting max.VT to 10 can improve the speed for calculation without affecting the power too much. The default parameter is NULL, which does not set upper limit on the number of variable frequency threhsold.
#' @param correctFlip Correct for flipped sites for score statistics and their covariance matrices
#' @return a list consisting of results;
#' @export
rareMETALS.range.group <- function(score.stat.file,cov.file,range,range.name,test='GRANVIL',refaltList,maf.cutoff=1,alternative=c('two.sided','greater','less'),out.digits=4,callrate.cutoff=0,hwe.cutoff=0,max.VT=NULL,correctFlip=TRUE,analyzeRefAltListOnly=TRUE)
{

    cat("Analyzing ",range.name[1],"\n");

    ii <- 1;

    res <- rareMETALS.range.group.core(score.stat.file,cov.file,range[1],range.name[1],test,refaltList,maf.cutoff,alternative,out.digits,callrate.cutoff,hwe.cutoff,max.VT,correctFlip,analyzeRefAltListOnly)
    if(length(range)>1)
        {
            for(ii in 2:length(range))
                {
                    cat("Analyzing ",range.name[ii],"\n");
                    res0 <- rareMETALS.range.group.core(score.stat.file,cov.file,range[ii],range.name[ii],test,refaltList,maf.cutoff,alternative,out.digits,callrate.cutoff,hwe.cutoff,max.VT,correctFlip,analyzeRefAltListOnly);                     
                    res$res.out <- rbind(res$res.out,res0$res.out);
                    res$res.list <- c(res$res.list,res0$res.list);
                    res$integratedData <- c(res$integratedData,res0$integratedData);
                }
        }
    return(res);
}
#' Meta-analysis of gene-level tests by range;
#'
#' @param score.stat.file files of score statistics
#' @param cov.file Covariance matrix files
#' @param range tabix range for each gene/region; Tabix range needs to be in the format of "1:123-1234". Quotation marks are necessary. 
#' @param range.name The name of the range,e.g. gene names can be used;
#' @param refaltList A list of reference, alternative allele, a vector of alternative allele frequencies; Specifically, the list should consist of 
#' @param test rare variant tests to be used
#' @param maf.cutoff MAF cutoff used to analyze variants; Default value is 1, i.e. no cutoffs are applied. MAFs are based upon the sample MAFs. 
#' @param alternative alternative hypothesis to be specified; Only applicable to GRANVIL test;
#' @param out.digits Number of digits used in the output
#' @param callrate.cutoff Cutoffs of call rate, lower than which will NOT be analyzed (labelled as missing)
#' @param hwe.cutoff Cutoffs of HWE p-values
#' @param max.VT The maximum number of thresholds used in VT; For small p-values, the calculation of VT p-values can be very slow. Setting max.VT to 10 can improve the speed for calculation without affecting the power too much. The default parameter is NULL, which does not set upper limit on the number of variable frequency threhsold. 
#' @return a list consisting of results;
#' @export
rareMETALS.range.group.core <- function(score.stat.file,cov.file,range,range.name,test='GRANVIL',refaltList,maf.cutoff=1,alternative=c('two.sided','greater','less'),out.digits=4,callrate.cutoff=0,hwe.cutoff=0,max.VT=NULL,correctFlip=TRUE,analyzeRefAltListOnly=TRUE)
  {
    ANNO <- "gene";gene.name <- range.name;
    no.boot=0;alpha=0.05;ix.gold=1
    raw.data.all <- list();
    for(ii in 1:length(range))
      {
        capture.output(raw.data.all[[ii]] <- rvmeta.readDataByRange( score.stat.file, cov.file, range[ii])[[1]]);
      }
    res.null <- list(gene.name=NA,
                     p.value=NA,
                     statistic=NA,
                     no.var=NA,
                     no.sample=NA,
                     anno=NA,
                     ref=NA,
                     raw.data=raw.data.all,
                     alt=NA,
                     maf.cutoff=NA,
                     direction.burden.by.study=NA,
                     direction.meta.single.var=NA,
                     beta1.est=NA,
                     no.site=NA,
                     beta1.sd=NA,
                     hsq.est=NA,
                     pos=NA);
    if(length(raw.data.all)==0)
      {
        res <- list();
        for(ii in 1:length(range.name))
          {
            res.null$gene.name <- range.name[ii];
            res[[ii]] <- res.null;
          }
        
        return(res);        
      }
    res <- list(list(p.value=NA,statistic=NA,no.var=NA,no.sample=NA));
    
    p.value <- double(0);
    ref <- character(0);
    alt <- ref;
    direction.by.study <- character(0);
    statistic <- double(0);pos <- integer(0);anno <- character(0);direction <- integer(0);res.maf.vec <- double(0);beta1.est.vec <- double(0);beta1.sd.vec <- double(0);
    gene.name.out <- range.name;p.value.out <- rep(NA,length(range.name));statistic.out <- rep(NA,length(range.name));no.site.out <- rep(NA,length(range.name));beta1.est.out <- rep(NA,length(range.name));beta1.sd.out <- rep(NA,length(range.name));maf.cutoff.out <- rep(NA,length(range.name));direction.burden.by.study.out <- rep(NA,length(range.name));direction.meta.single.var.out <- rep(NA,length(range.name));pos.ref.alt.out <- rep(NA,length(range.name));top.singlevar.pval <- rep(NA,length(range.name));top.singlevar.refalt <- rep(NA,length(range.name));top.singlevar.pos <- rep(NA,length(range.name));top.singlevar.af <- rep(NA,length(range.name));integratedData <- list();N.out <- 0;
    refaltList.ori <- refaltList;
    for(kk in 1:length(raw.data.all))
      {
        raw.data <- raw.data.all[[kk]];
        raw.data.ori <- raw.data;
        refaltList <- refaltList.ori;
        raw.data$gene <- range.name[kk];
        QC.par <- list(callrate.cutoff=callrate.cutoff,hwe.cutoff=hwe.cutoff);
        raw.data <- QC(raw.data,QC.par,cov=1);        
        ix.var <- integer(0);
        log.mat <- raw.data$log.mat;
       
        if(length(raw.data$ustat[[1]])>0)
            ix.var <- 1:length(raw.data$ustat[[1]]);
        
        ix.var <- sort(unique(ix.var));
        if(length(ix.var)==0)
          {
            res.null$gene.name <- gene.name[kk];
            res[[kk]] <- res.null;
          }
        if(length(ix.var)>=1) {
          ix.pop <- 1:length(raw.data$nSample);
          score.stat.vec.list <- list();mac.vec.list <- list();maf.vec.list <- list();cov.mat.list <- list();var.Y.list <- list();N.list <- list();mean.Y.list <- list();pos.list <- list();anno.list <- list();ac.vec.list <- list();af.vec.list <- list();ref.list <- list();alt.list <- list();
          no.sample <- 0;nref.list <- list();nalt.list <- list();nhet.list <- list();
          ix.match <- match(raw.data$pos[ix.var],refaltList$pos);
          ref.gold <- refaltList$ref[ix.match];alt.gold <- refaltList$alt[ix.match];af.gold <- refaltList$af[ix.match];checkAF <- refaltList$checkAF;anno.gold <- refaltList$anno[ix.match];pos.gold <- refaltList$pos[ix.match];
          if(length(checkAF)==0) checkAF <- FALSE;
          refaltList <- list(ref=ref.gold,alt=alt.gold,af=af.gold,checkAF=checkAF,af.diff.max=refaltList$af.diff.max,anno=anno.gold,pos=pos.gold);
          
          for(ii in 1:length(ix.pop))
              {
                  for(jj in 1:length(ix.var))
                      {
                          res.flipAllele <- flipAllele(raw.data,raw.data.ori,refaltList,ii,ix.var[jj],log.mat[ix.var[jj],],correctFlip,analyzeRefAltListOnly);

                          raw.data <- res.flipAllele$raw.data;
                          log.mat[ix.var[jj],] <- res.flipAllele$log.mat.var;

                      }
              }
       
          af.mat <- matrix(NA,ncol=length(raw.data$ustat[[ii]]),nrow=length(raw.data$ustat));
          ac.mat <- af.mat;
          N.mat <- af.mat;
          vstat.list <- list();
          for(ii in 1:length(ix.pop))
              {
                  if(length(raw.data$afCase[[ii]])>0) {
                      warning(paste0("Study  ",ii," is analyzed as binary trait. It is advised to use rareMETALS2 for meta-analysis"))
                  }
                  U.stat <- raw.data$ustat[[ii]][ix.var];
                  V.stat <- raw.data$vstat[[ii]][ix.var];
                  vstat.list[[ii]] <- V.stat;
                  N.list[[ii]] <- rm.na(as.integer(mean(raw.data$nSample[[ii]],na.rm=TRUE)));
                  no.sample <- no.sample+N.list[[ii]];
                  N.mat[ii,] <- raw.data$nSample[[ii]];
                  score.stat.vec.list[[ii]] <- (U.stat/V.stat);              
                  cov.mat.list[[ii]] <- matrix((raw.data$cov[[ii]])[ix.var,ix.var],nrow=length(ix.var),ncol=length(ix.var));
                  var.Y.list[[ii]] <- 1;
                  abc <- cbind(V.stat^2,diag(cov.mat.list[[ii]]*N.list[[ii]]));
                  
                  mean.Y.list[[ii]] <- 0;
                  af.vec.list[[ii]] <- ((raw.data$af[[ii]])[ix.var]);
                  ac.vec.list[[ii]] <- ((raw.data$ac[[ii]])[ix.var]);
                  af.mat[ii,] <- af.vec.list[[ii]][ix.var];
                  ac.mat[ii,] <- ac.vec.list[[ii]][ix.var];
                  pos.list[[ii]] <- (raw.data$pos)[ix.var];
                  ref.list[[ii]] <- (raw.data$ref)[[ii]][ix.var];
                  alt.list[[ii]] <- (raw.data$alt)[[ii]][ix.var];
                  anno.list[[ii]] <- (raw.data$anno)[ix.var];
                  nref.list[[ii]] <- rm.na((raw.data$nref)[[ii]])[ix.var];
                  nalt.list[[ii]] <- rm.na((raw.data$nalt)[[ii]])[ix.var];
                  nhet.list[[ii]] <- rm.na((raw.data$nhet)[[ii]])[ix.var];              
              }

          maf.vec <- rep(0,length(af.vec.list[[1]]));
          af.vec <- maf.vec;
          mac.vec <- 0;ac.vec <- 0;
          af.vec <- colSums(af.mat*N.mat,na.rm=TRUE)/colSums(N.mat,na.rm=TRUE);
          ac.vec <- colSums(ac.mat,na.rm=TRUE);

          maf.vec <- af.vec;
          mac.vec <- ac.vec;
          ix.major <- which(af.vec>0.5);
          if(length(ix.major)>0)
            {
              maf.vec[ix.major] <- 1-maf.vec[ix.major];
              mac.vec[ix.major] <- 2*no.sample-ac.vec[ix.major];
              tmp.major <- flip.score.cov(score.stat.vec.list,cov.mat.list,ix.major);
              score.stat.vec.list <- tmp.major$score.stat.vec.list;
              cov.mat.list <- tmp.major$cov.mat.list;
            }
          ix.rare <- which(maf.vec<maf.cutoff & maf.vec>0);
          maf.vec.rare <- maf.vec[ix.rare];
          mac.vec.rare <- mac.vec[ix.rare];
          refaltList.rare <- list(ref=refaltList$ref[ix.rare],
                                  alt=refaltList$alt[ix.rare],
                                  af=refaltList$af[ix.rare],
                                  anno=refaltList$anno[ix.rare],
                                  af.diff.max=refaltList$af.diff.max,
                                  pos=refaltList$pos[ix.rare]);
          if(length(ix.rare)==0)
            {
              res.null$gene.name <- gene.name[kk];
              res[[kk]] <- res.null;
            }
          if(length(ix.rare)>=1)
            {
              for(ii in 1:length(ix.pop))
                {
                    score.stat.vec.list[[ii]] <- score.stat.vec.list[[ii]][ix.rare];
                    vstat.list[[ii]] <- vstat.list[[ii]][ix.rare];
                    cov.mat.list[[ii]] <- matrix(cov.mat.list[[ii]][ix.rare,ix.rare],nrow=length(ix.rare),ncol=length(ix.rare));
                    af.vec.list[[ii]] <- af.vec.list[[ii]][ix.rare];
                    ac.vec.list[[ii]] <- ac.vec.list[[ii]][ix.rare];
                    anno.list[[ii]] <- anno.list[[ix.gold]][ix.rare];
                    pos.list[[ii]] <- pos.list[[ix.gold]][ix.rare];
                    ref.list[[ii]] <- ref.list[[ix.gold]][ix.rare];
                    alt.list[[ii]] <- alt.list[[ix.gold]][ix.rare];
                    
              }
              
              N.mat <- matrix(N.mat[,ix.rare],nrow=nrow(N.mat),ncol=length(ix.rare));
              covG <- matrix(0,nrow=length(ix.rare),ncol=length(ix.rare));
              nSample.covG <- covG;
              
              for(aa in 1:length(ix.rare)) {
                  for(bb in 1:length(ix.rare)) {
                      for(cc in ix.pop) {
                          
                          covG[aa,bb] <- covG[aa,bb]+rm.na(sqrt(N.mat[cc,aa]*N.mat[cc,bb])*cov.mat.list[[cc]][aa,bb]);
                          nSample.covG[aa,bb] <- nSample.covG[aa,bb]+sqrt(rm.na(N.mat[cc,aa])*rm.na(N.mat[cc,bb]));
                      }
                  }
              }
              r2.approx <- cov2cor(rm.na(covG/nSample.covG));

              res.extra <- list(anno=anno.list[[ix.gold]],
                                pos=refaltList.rare$pos,
                                ref=refaltList.rare$ref,
                                alt=refaltList.rare$alt,
                                af.vec.list=af.vec.list,
                                score.stat.vec.list=score.stat.vec.list,
                                cov.mat.list=cov.mat.list,
                                pos.list=pos.list,
                                alt.list=alt.list,
                                ref.list=ref.list,
                                af.mat=af.mat,
                                ac.mat=ac.mat,
                                r2.approx=r2.approx,
                                log.mat=log.mat,
                                N.mat=N.mat,
                                raw.data=raw.data.ori,
                                clean.data=raw.data,
                                gene.name=gene.name[kk]);
              if(test=='WSS')
                  {
                  res.kk <- (c(rvmeta.CMH(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.mat,alternative,no.boot,alpha,rv.test='WSS',extra.pars=list(vstat.list=vstat.list,weight='MB',ac.vec.list=ac.vec.list,
                                                                                                                                                    maf.vec=maf.vec.rare,mac.vec=mac.vec.rare))));
                  res[[kk]] <- c(res.kk,res.extra);
                  res[[kk]]$maf.cutoff <- maf.cutoff;
                }
              if(test=='GRANVIL')
                {
                  res.kk <- (c(rvmeta.CMH(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.mat,alternative,no.boot,alpha,rv.test='WSS',extra.pars=list(vstat.list=vstat.list,weight='MZ',ac.vec.list=ac.vec.list,
                                                                                                                                                    maf.vec=maf.vec.rare,mac.vec=mac.vec.rare))));
                  res[[kk]] <- c(res.kk,res.extra);
                  res[[kk]]$maf.cutoff <- maf.cutoff;
                }
              if(test=='SKAT')
                {
                  res.kk <- (c(rvmeta.CMH(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.mat,alternative,no.boot,alpha,rv.test='SKAT',extra.pars=list(vstat.list=vstat.list,kernel='beta',ac.vec.list=ac.vec.list,
                                                                                                                                                     maf.vec=maf.vec.rare,mac.vec=mac.vec.rare))));
                  res[[kk]] <- c(res.kk,res.extra);
                  res[[kk]]$maf.cutoff <- maf.cutoff;
                }
              if(test=='VT')
                {
                  res.kk <- (c(rvmeta.CMH(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.mat,alternative,no.boot,alpha,rv.test='VT',extra.pars=list(vstat.list=vstat.list,ac.vec.list=ac.vec.list,maf.vec=maf.vec.rare,mac.vec=mac.vec.rare,max.TH=max.VT))));
                  pos.VT <- res.extra$pos[res.kk$ixVar.VT];               
                  res[[kk]] <- c(res.kk,res.extra,list(pos.VT=pos.VT));
                }
              res[[kk]]$vstat.list <- vstat.list;
              gene.name.out[kk] <- res[[kk]]$gene.name;
              N.out[kk] <- res[[kk]]$N;
              p.value.out[kk] <- format(res[[kk]]$p.value,digits=out.digits);
              statistic.out[kk] <- format(res[[kk]]$statistic,digits=out.digits);
              no.site.out[kk] <- res[[kk]]$no.site;
              beta1.est.out[kk] <- format(res[[kk]]$beta1.est,digits=out.digits);
              beta1.sd.out[kk] <- format(res[[kk]]$beta1.sd,digits=out.digits);
              maf.cutoff.out[kk] <- format(res[[kk]]$maf.cutoff,digits=out.digits);
              direction.burden.by.study.out[kk] <- res[[kk]]$direction.burden.by.study;
              direction.meta.single.var.out[kk] <- res[[kk]]$direction.meta.single.var;          
              pos.ref.alt.out[kk] <- paste(res[[kk]]$pos,res[[kk]]$ref,res[[kk]]$alt,sep='/',collapse=',');
              ix.best <- res[[kk]]$ix.best;
              top.singlevar.pos[kk] <- res[[kk]]$pos[ix.best];
              top.singlevar.refalt[kk] <- paste(c(res[[kk]]$ref[ix.best],res[[kk]]$alt[ix.best]),sep="/",collapse="/");
              top.singlevar.pval[kk] <- format(res[[kk]]$singlevar.pval.vec[ix.best],digits=out.digits);
              top.singlevar.af[kk] <- format(res[[kk]]$singlevar.af.vec[ix.best],digits=out.digits);
              ref.out <- list(res[[kk]]$ref);
              alt.out <- list(res[[kk]]$alt);
              nSample.out <- list(rep(res[[kk]]$nSample,length(ref.out[[1]])));
              af.out <- list(res[[kk]]$singlevar.af.vec);
              ac.out <- list(as.integer(af.out[[1]]*nSample.out[[1]]*2));
              callrate <- list(rep(1,length(ref.out[[1]])));
              nref.out <- 0;
              nalt.out <- 0;
              nhet.out <- 0;
              ustat.out <- list(res[[kk]]$ustat);
              vstat.out <- list(res[[kk]]$vstat);

              callrate.out <- list(rep(1,length(ref.out[[1]])));
              hwe.out <- list(rep(1,length(ref.out[[1]])));
              effect.out <- list(ustat.out[[1]]/vstat.out[[1]]^2);
              pVal.out <- list(res[[kk]]$singlevar.pval.vec);
              cov.out <- list(res[[kk]]$cov);
              pos.out <- res[[kk]]$pos;
              anno.out <- res[[kk]]$anno;
              res[[kk]]$cMAC <- sum(mac.vec.rare);
              res[[kk]]$ix.rare <- ix.rare;
              for(aa in 1:length(nref.list))
                {
                  nref.out <- nref.out+nref.list[[aa]];
                  nalt.out <- nalt.out+nalt.list[[aa]];
                  nhet.out <- nhet.out+nhet.list[[aa]];
                }
              nref.out <- list(nref.out);
              nalt.out <- list(nalt.out);
              nhet.out <- list(nhet.out);
              integratedData[[kk]] <- list(ref=ref.out,
                                           alt=alt.out,
                                           nSample=nSample.out,
                                           af=af.out,
                                           ac=ac.out,
                                           callrate=callrate.out,
                                           hwe=hwe.out,
                                           nref=nref.out,
                                           nalt=nalt.out,
                                           ustat=ustat.out,
                                           vstat=vstat.out,
                                           nhet=nhet.out,
                                           effect=effect.out,
                                           pVal=pVal.out,
                                           cov=cov.out,
                                           pos=pos.out,
                                           covZZ=list(matrix(nrow=0,ncol=0)),
                                           covXZ=list(matrix(ncol=0,nrow=length(ref.out[[1]]))),
                                           hweCase=list(rep(NA,length(ref.out[[1]]))),
                                           hweCtrl=list(rep(NA,length(ref.out[[1]]))),
                                           afCtrl=list(rep(NA,length(ref.out[[1]]))),
                                           afCase=list(rep(NA,length(ref.out[[1]]))),
                                           anno=anno.out);
            }
        }
    }
    res.out <- cbind(gene.name.out,N.out,p.value.out,statistic.out,no.site.out,beta1.est.out,beta1.sd.out,maf.cutoff.out,direction.burden.by.study.out,direction.meta.single.var.out,top.singlevar.pos,top.singlevar.refalt,top.singlevar.pval,top.singlevar.af,pos.ref.alt.out);
    
    return(list(res.list=res,
                integratedData=integratedData,
                res.out=res.out));
  }

r2cov.mat <- function(r2.mat,maf.vec)
  {
    var.vec <- sqrt(maf.vec*(1-maf.vec)*2);
    var.mat <- (var.vec%*%t(var.vec));
    cov.mat <- r2.mat*var.mat;
    return(cov.mat);
  }
flip.score.cov <- function(score.stat.vec.list,cov.mat.list,ix.major)
  {
    for(ii in 1:length(score.stat.vec.list))
      {
        score.stat.vec.list[[ii]][ix.major] <- (-1)*score.stat.vec.list[[ii]][ix.major];
        cov.mat.list[[ii]][ix.major,] <- (-1)*cov.mat.list[[ii]][ix.major,];
        cov.mat.list[[ii]][,ix.major] <- (-1)*cov.mat.list[[ii]][,ix.major];
      }
    return(list(cov.mat.list=cov.mat.list,
                score.stat.vec.list=score.stat.vec.list));
  }
check.mono <- function(af.vec.list,ac.vec.list,N.list)
  {
    count <- rep(0,length(af.vec.list[[1]]));
    for(ii in 1:length(af.vec.list))
      {
        count <- count+(af.vec.list[[ii]]==1 | af.vec.list[[ii]]==0);
      }
    ix.change <- which(count==length(af.vec.list));
    if(length(ix.change)>0)
      {
        for(ii in 1:length(af.vec.list))
          {
            af.vec.list[[ii]][ix.change] <- 0;
            ac.vec.list[[ii]][ix.change] <- N.list[[ii]]
          }
      }
    return(list(af.vec.list=af.vec.list,
                ac.vec.list=ac.vec.list));
  }
