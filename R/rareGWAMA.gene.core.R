#' organize the formatted stat into analyzable format;
#'
#' @param dat formatted data;
#' @param maf.cutoff The minor allele frequency cutoff
#' @return dat which consist of formatted data for gene-level tests
#' @export
rareGWAMA.formatGene <- function(dat,...) {
    extraPar <- list(...);
    maf.cutoff <- extraPar$maf.cutoff;
    method <- extraPar$method;
    if(is.null(method)) method <- "pseudoscore";
    if(is.null(maf.cutoff)) maf.cutoff <- 0.05;
    lambda <- extraPar$lambda;
    if(is.null(lambda)) lambda <- 0.1;
    af.meta <- rowSums(dat$af.mat*dat$nSample.mat,na.rm=TRUE)/rowSums(dat$nSample.mat,na.rm=TRUE)
    maf.meta <- rm.na(af.meta);
    
    maf.meta[which(maf.meta>.5)] <- 1-maf.meta[which(maf.meta>.5)];
    ix.rare <- which(maf.meta<maf.cutoff);
    dat <- NULL;
    if(length(ix.rare)>0) {
        pos.rare <- dat$pos[ix.rare];
        r2.rare <- as.matrix(dat$r2[ix.rare,ix.rare]);
        ustat.meta <- rowSums(matrix(dat$ustat.mat[ix.rare,],nrow=length(ix.rare)),na.rm=TRUE);
        
        N.ustat.meta <- rowSums(matrix(dat$nSample.mat[ix.rare,],nrow=length(ix.rare)),na.rm=TRUE);
        V.meta <- matrix(0,nrow=length(ix.rare),ncol=length(ix.rare));
        V.meta.reg <- V.meta;
        N.V.meta <- matrix(0,nrow=length(ix.rare),ncol=length(ix.rare));
        for(ii in 1:ncol(dat$ustat.mat)) {
            sd.mat <- matrix(0,nrow=length(ix.rare),ncol=length(ix.rare));
            diag(sd.mat) <- rm.na(dat$vstat.mat[ix.rare,ii]);
            V.meta <- V.meta+sd.mat%*%r2.rare%*%sd.mat;
            V.meta.reg <- V.meta.reg+sd.mat%*%(cov2cor(r2.rare+lambda*diag(nrow(r2.rare))))%*%sd.mat;
            sqrt.nSample <- as.matrix(sqrt(rm.na(dat$nSample.mat[ix.rare,ii])));
            N.V.meta <- N.V.meta+sqrt.nSample%*%t(sqrt.nSample);
        }

        beta.est <- ginv(V.meta)%*%ustat.meta;
        beta.var <- ginv(V.meta);

        beta.est.reg <- ginv(V.meta.reg)%*%ustat.meta;
        beta.var.reg <- ginv(V.meta.reg);
        
        V.meta.weighted <- rm.na(V.meta/N.V.meta);
        V.meta.weighted.reg <- rm.na(V.meta.reg/N.V.meta);
        
        ustat.meta.weighted <- ustat.meta/N.ustat.meta;
        
        weight <- matrix(0,nrow=length(ustat.meta.weighted),ncol=length(ustat.meta.weighted));
        diag(weight) <- 1/N.ustat.meta;
        var.ustat.meta.weighted <- weight%*%V.meta%*%weight;
        if(method=='pseudoscore') {
            beta.est <- ginv(V.meta.weighted)%*%ustat.meta.weighted;
            beta.var <- ginv(V.meta.weighted)%*%var.ustat.meta.weighted%*%ginv(V.meta.weighted);
            beta.est.reg <- ginv(V.meta.weighted.reg)%*%ustat.meta.weighted;
            beta.var.reg <- ginv(V.meta.weighted.reg)%*%var.ustat.meta.weighted%*%ginv(V.meta.weighted.reg);

        }
        direction.singlevar.assoc <- rep("X",length(ustat.meta));
        direction.singlevar.assoc[which(ustat.meta>0)] <- "+";
        direction.singlevar.assoc[which(ustat.meta<0)] <- "-";
        
        dat$ustat.meta <- ustat.meta;
        dat$V.meta <- V.meta;
        dat$direction.singlevar.assoc <- direction.singlevar.assoc;
        dat$pos.rare <- pos.rare;
        dat$maf.meta <- maf.meta;
        dat$nSample.meta <- N.ustat.meta;
        dat$beta.est <- beta.est;
        dat$beta.var <- beta.var;
        dat$ustat.meta.weighted <- ustat.meta.weighted;
        dat$N.ustat.meta <- N.ustat.meta;
        dat$N.V.meta <- N.V.meta;
        dat$V.meta.weighted <- V.meta.weighted;
        dat$r2.rare <- r2.rare;
        dat$maf.rare <- maf.meta[ix.rare];
        dat$beta.est.reg <- beta.est.reg;
        dat$beta.var.reg <- beta.var.reg;
        sdG <- matrix(0,nrow=length(dat$maf.rare),ncol=length(dat$maf.rare));
        diag(sdG) <- sqrt(2*(dat$maf.rare)*(1-dat$maf.rare));
        covG <- sdG%*%dat$r2.rare%*%sdG;
        dat$h2 <- as.numeric(t(dat$beta.est.reg)%*%covG%*%(dat$beta.est.reg));
        h2.var <- as.numeric(4*t(beta.est.reg)%*%covG%*%beta.var.reg%*%covG%*%beta.est.reg);
        dat$h2.sd <- sqrt(h2.var);
        
    }
    return(dat);
    
}

#' Conduct approximate gene-level tests and estimate variance explained; 
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param vcf.ref.file the file names of the reference panel file;
#' @param tabix.range Tabix range for variants;
#' @return formatted data and assocation testing results; 
#' @export
rareGWAMA.gene <- function(score.stat.file,imp.qual.file=NULL,vcf.ref.file,tabix.range,...) {
    uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
    extraPar <- list(...);
    maf.cutoff <- extraPar$maf.cutoff;
    if(is.null(maf.cutoff)) maf.cutoff <- 1;
    r2.cutoff <- extraPar$r2.cutoff;
    if(is.null(r2.cutoff)) r2.cutoff <- 0.95;
    sizePerBatch <- extraPar$sizePerBatch;
    if(is.null(sizePerBatch)) sizePerBatch <- 100;
    refGeno <- extraPar$refGeno;
    col.impqual <- extraPar$col.impqual;

    impQual.lb <- extraPar$impQual.lb;
    impQualWeight <- FALSE;
    rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
    if(is.null(col.impqual)) col.impqual <- 5;
    if(is.null(impQual.lb)) impQual.lb <- 0.7;
    if(is.null(rmMultiAllelicSite)) rmMultiAllelicSite <- TRUE;
    if(is.null(refGeno)) refGeno <- "DS";
    pseudoScore <- extraPar$pseudoScore;
    if(is.null(pseudoScore)) pseudoScore <- TRUE;
    beta.est <- 0;beta.se <- 0;statistic <- 0;p.value <- 0;ref.tab <- 0;alt.tab <- 0;pos.all <- 0;marginal.statistic <- 0;marginal.p.value <- 0;
    ii <- 0;batchEnd <- 0;batchStart <- 0;nSample <- 0;af <- 0;numStability <- 0;
    a <- Sys.time();
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range,multiAllelic = TRUE));
    vcfIndv <- refGeno;
    annoType <- "";
    vcfColumn <- c("CHROM","POS","REF","ALT");
    vcfInfo <- NULL;      
    geno.list <- readVCFToListByRange(vcf.ref.file, tabix.range, "", vcfColumn, vcfInfo, vcfIndv)        
        raw.imp.qual <- NULL;
    if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
    time.readData <- Sys.time()-a;
    b <- Sys.time();
    raw.data.all <- raw.data.all[[1]];
    cat('Read in',length(raw.data.all$ref[[1]]),'variants\n',sep=' ');
    dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual);
    if(rmMultiAllelicSite==TRUE) {
        tmp <- GWAMA.rmMulti(dat);
        dat <- tmp$dat;posMulti <- tmp$posMulti;
    }
    pos <- gsub("_.*","",dat$pos);
    if(refGeno=="DS") {
        gt <- geno.list$DS;
        gt <- matrix(as.numeric(gt),nrow=nrow(gt),ncol=ncol(gt));
    }
    if(refGeno=="GT") {
        gt.tmp <- geno.list$GT;
        gt <- matrix(NA,nrow=nrow(gt.tmp),ncol=ncol(gt.tmp));
        gt[which(gt.tmp=="0/0",arr.ind=T)] <- 0;
        gt[which(gt.tmp=="1/0",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="0/1",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="1/1",arr.ind=T)] <- 2
        gt[which(gt.tmp=="0|0",arr.ind=T)] <- 0;
        gt[which(gt.tmp=="1|0",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="0|1",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="1|1",arr.ind=T)] <- 2
        
    }
    
    r2.tmp <- cor(gt,use='pairwise.complete');
    r2.tmp <- rm.na(r2.tmp);
    pos.vcf <- paste(geno.list$CHROM,geno.list$POS,sep=":");
    r2 <- matrix(0,nrow=length(pos),ncol=length(pos));
    r2 <- as.matrix(r2.tmp[match(pos,pos.vcf),match(pos,pos.vcf)]);
    colnames(r2) <- pos;
    diag(r2) <- 1;
    r2 <- rm.na(r2);
    dat$r2 <- r2;
    dat <- rareGWAMA.formatGene(dat,maf.cutoff=maf.cutoff,lambda=extraPar$lambda);
    return(list(dat=dat));
    
    
}


#' Estimate the variance explained by the variants in a locus;
#'
#' @param r2.rare the r2 matrix for rare variants;
#' @param maf.rare the maf vector for rare variants;
#' @return r2 estimates;
#' @export
estimateH2 <- function(...) {
    dat <- list(...);
    r2 <- rm.na(dat$r2.rare)
    maf <- rm.na(dat$maf.rare)
   
    beta.est <- rm.na(as.vector(dat$beta.est))
    beta.var <- rm.na(dat$beta.var);
   
    h2 <- 0;h2.var <- 0;
    
    if(!is.null(r2)) {
        sd.mat <- matrix(0,nrow=length(maf),ncol=length(maf));
        diag(sd.mat) <- sqrt(2*maf*(1-maf));
        v.mat <- sd.mat%*%r2%*%sd.mat
        h2 <- t(beta.est)%*%v.mat%*%beta.est;
        svd.v <- svd(v.mat);
        svd.omega <- svd(beta.var%*%v.mat%*%beta.var);
        h2.var <- sum(2*(svd.v$d)^2*svd.omega$d);
    }
    return(list(h2.est=as.numeric(h2),
                h2.var=h2.var));
    
}

#' gene-level test
#'
#' @param dat
#' @return a list with statistics and p-values;
#' @export
rareGWAMA.burden <- function(dat,...) {
    if(is.null(dat$ustat.meta)) return(list(statistic=NULL, p.value=NULL));

    ustat.burden <- sum(dat$ustat.meta,na.rm=TRUE);
    vstat.sq.burden <- (sum(dat$V.meta,na.rm=TRUE));
    statistic <- ustat.burden^2/vstat.sq.burden;
    p.value <- pchisq(statistic,df=1,lower.tail=FALSE);
    return(list(statistic=statistic,
                p.value=p.value));
}



#' t test;
#'
#' @param dat;
#' @return a list with statistic and p-value;
#' @export
rareGWAMA.t <- function(dat,...) {
    statistic <- t(dat$ustat.meta)%*%ginv(dat$V.meta)%*%dat$ustat.meta;
    p.value <- pchisq(statistic,df=length(dat$ustat.meta),lower.tail=FALSE);
    return(list(statistic=statistic,
                p.value=p.value));

}


#' SKAT test
#' @param dat
#' @return a list with statistics and p-values;
#' @export
rareGWAMA.skat <- function(dat,...) {
    extraPar <- list(...);
    weight <- extraPar$weight;
    if(is.null(weight)) weight <- 'beta';
    
    if(is.null(dat$ustat.meta)) return(list(statistic=NULL, p.value=NULL));
    W <- matrix(0,nrow=length(dat$ustat.meta),ncol=length(dat$ustat.meta));
    
    if(weight=='beta')  diag(W) <- dbeta(dat$maf.meta,1,25)^2;
    if(weight=='linear') diag(W) <- 1;
    Q <- sum((dat$ustat.meta)^2*diag(W));
    svd.V <- svd(dat$V.meta);
    lambda.V <- matrix(0,nrow=length(abs(svd.V$d)),ncol=length(abs(svd.V$d)));
    diag(lambda.V) <- abs(svd.V$d);
    
    L <- (svd.V$u)%*%(sqrt(lambda.V))%*%t(svd.V$v);
    lambda <- try(get.eigen(W,L,t(L)),silent=TRUE);
    if(class(lambda)=='try-error')
    {
        return(list(statistic=NA,
                    p.value=NA));
    }
    p.value.liu <- NA;
    p.value.davies <- NA;
    p.value.imhof <- NA;
    p.value <- NA;
    p.value.davies <- try(davies(Q,lambda=lambda)$Qq,silent=TRUE);
    p.value.liu <- try(liu(Q,lambda=lambda),silent=TRUE);
    p.value.imhof <- try(imhof(Q,lambda=lambda)$Qq,silent=TRUE);
    if(length(attr(p.value.davies,'class'))+length(attr(p.value.liu,'class'))+length(attr(p.value.imhof,'class'))>0)
        return(list(statistic=NA,
                    p.value=NA));
    p.value <- p.value.davies;
    if(p.value<=0 | p.value>=1) p.value <- p.value.liu;
    return(list(statistic=Q,
                p.value=p.value));

}

#' vt test;
#'
#' @param dat
#' @return a list with statistics and p-values;
#' @export
rareGWAMA.vt <- function(dat,...) {
    if(is.null(dat$ustat.meta)) return(list(statistic=NULL, p.value=NULL));
    extraPars <- list(...);
    maf.vec <- dat$maf.meta;
    maf.TH <- sort(unique(maf.vec));
    max.TH <- extraPars$max.TH;
    ustat.meta <- dat$ustat.meta
    alternative <- "two.sided";
    if(length(max.TH)>0)
    {
        ix.tmp <- as.integer(seq(1,length(maf.TH),length=max.TH));
        maf.TH <- maf.TH[ix.tmp];
        maf.TH <- unique(sort(maf.TH));
    }
    maf.TH.old <- maf.TH;
    if(maf.TH[1]==0) {maf.TH <- maf.TH[-1];}
    if(length(maf.TH)==0)
    {
        return(list(p.value=NA,
                    statistic=NA,
                    maf.cutoff.vt=0,
                    no.site.VT=length(which(mac.vec==0))));
    }
    if(length(maf.TH)==1)
    {
        ix.ii <- which(maf.vec<=maf.TH);    
        dat$ustat.meta <- sum(ustat.meta[ix.ii]);
        ind.ii <- rep(0,length(maf.vec));
        ind.ii[ix.ii] <- 1;
        V.stat.sq <- as.numeric(t(ind.ii)%*%dat$V.meta%*%(ind.ii));
 
        statistic <- dat$ustat.meta^2/V.stat.sq;
        p.value <- pchisq(statistic,df=1,lower.tail=FALSE);
        
        return(list(p.value=p.value,
                    statistic=statistic,
                    maf.cutoff.vt=maf.TH,
                    no.site.VT=length(ix.ii)));
    }
    err.msg <- vector(length=0);
    ix.list <- list();
    ustat.meta.VT <- 0;
    VT.mat <- matrix(0,nrow=length(maf.TH),ncol=ncol(dat$V.meta));
    no.TH <- length(maf.TH);
    for(ii in 1:length(maf.TH))
    {
        ix.ii <- which(maf.vec<=maf.TH[ii]);
        ix.list[[ii]] <- ix.ii;
        VT.mat[ii,ix.ii] <- 1;
        ustat.meta.VT[ii] <- sum(ustat.meta[ix.ii]);
    }
    cov.X.VT <- matrix(0,nrow=length(ix.list),ncol=length(ix.list));
    cov.mat <- dat$V.meta;
    a=Sys.time();
    for(ii in 1:length(maf.TH))
    {
        for(jj in 1:length(maf.TH))
        {
            ix.ii <- rep(0,nrow(cov.mat));
            ix.jj <- rep(0,nrow(cov.mat));
            ix.ii[ix.list[[ii]]] <- 1;
            ix.jj[ix.list[[jj]]] <- 1;
            cov.X.VT[ii,jj] <- as.numeric(t(ix.ii)%*%cov.mat%*%(ix.jj));
        }
      }
    
    cor.X.VT <- cov.X.VT/sqrt(diag(cov.X.VT)%*%t(diag(cov.X.VT)));
    
    vt.stat.vec <- (ustat.meta.VT^2/(diag(cov.X.VT)));
    ix.rm <- which(is.na(vt.stat.vec));
    if(length(ix.rm)>0) {
        vt.stat.vec <- vt.stat.vec[-ix.rm];
        cor.X.VT <- as.matrix(cor.X.VT[-ix.rm,-ix.rm]);          
        maf.TH <- maf.TH[-ix.rm];
    }
    ix.max <- which.max(vt.stat.vec);
    vt.max.stat <- vt.stat.vec[ix.max];
    maf.cutoff <- maf.TH[ix.max];
    if(length(ix.max)==0) {
        vt.max.stat <- NA;
        maf.cutoff <- maf.TH.old[length(maf.TH.old)];
    }
    p.value <- NA;
    if(!is.na(vt.max.stat)){
        p.value <- pvt(vt.max.stat,mu=rep(0,nrow(cor.X.VT)),sigma=cor.X.VT,alternative);
    }
    return(list(p.value=p.value,
                statistic=vt.max.stat,
                maf.cutoff.vt=maf.cutoff,
                no.site.VT=length(ix.ii)));
    
}
