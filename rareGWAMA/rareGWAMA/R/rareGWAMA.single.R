#' get imputation quality from tables based upon imputation quality files;
#'
#' @param imp.qual the list with imputation qualities
#' @param pos A vector of positions where the imputatoin quality will be retrieved;
#' @param col.impqual The column where the imputation quality will be retrieved;
#' @return a matrix with imputaton qualities
#' @export
getImpQual <- function(imp.qual,pos,col.impqual) {
    pos.imp <- try(paste(imp.qual[,1],imp.qual[,2],sep=":"),silent=TRUE);
    if(class(pos.imp)=='try-error') return(rep(NA,length(pos)))
    ix <- match(gsub("_.*","",pos),pos.imp);
    return(as.numeric(imp.qual[ix,col.impqual]));
}
#' unique alternative alleles;
#'
#' @param x A vector of allele labels;
#' @return A vector of unique alleles;
#' @export
uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
#' single variant meta-analysis integrating imputation quality;
#'
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param tabix.range the tabix range. IT must be in quote and provided as a string;
#' @param alternative The alternative hypothesis. Default is two.sided;
#' @param col.impqual The column number for the imputation quality score;
#' @param impQual.lb The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
#' @param impQualWeight Using imputation quality as weight
#' @param rmMultiAllelicSite Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function; 
#' @return A list of analysis results;
#' @export 
rareGWAMA.single <- function(score.stat.file,imp.qual.file=NULL,tabix.range,alternative="two.sided",col.impqual=5,impQual.lb=0.7,impQualWeight=FALSE,rmMultiAllelicSite=FALSE,gc=FALSE,...) {
    a <- Sys.time();
    extraPar <- list(...);
    weight <- extraPar$weight;
    if(is.null(weight)) weight <- "Npq+impQ";
    lambda <- rep(1,length(score.stat.file));
    var.y <- extraPar$var.y;
    if(is.null(var.y)) var.y <- 1;
    binaryTrait <- extraPar$binaryTrait;
    if(is.null(binaryTrait)) binaryTrait <- FALSE;
    if(gc==TRUE) {
        lambda <- extraPar$lambda;
        if(is.null(lambda)) lambda <- rep(1,length(score.stat.file));
    }
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range, multiAllelic = TRUE));
    raw.imp.qual <- NULL;
    if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
    time.readData <- Sys.time()-a;
    b <- Sys.time();
    raw.data.all <- raw.data.all[[1]];
    cat('Read in',length(raw.data.all$ref[[1]]),'variants',sep=' ');
    dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual);
    dat.withMulti <- dat;
    tmp <- GWAMA.rmMulti(dat);
    dat <- tmp$dat;posMulti <- tmp$posMulti;
    direction.meta <- apply(dat$direction.mat,1,paste,sep="",collapse="");   
    
    dat$ustat.mat <- scale(dat$ustat.mat,center=FALSE,scale=sqrt(lambda));
    ustat.meta <- rowSums((dat$ustat.mat)*(dat$w.mat),na.rm=TRUE);
    vstat.sq.meta <- rowSums((dat$w.mat)^2*(dat$vstat.mat)^2,na.rm=TRUE);
    beta.meta <- (ustat.meta)/(vstat.sq.meta);
    beta.sd.meta <- sqrt(1/vstat.sq.meta);
    if(weight=="N") {
        w.mat <- sqrt(dat$nSample);
        z.mat <- dat$ustat.mat/dat$vstat.mat;
        statistic.meta <- (rowSums(w.mat*z.mat,na.rm=TRUE))^2/rowSums(w.mat^2,na.rm=TRUE)
    }
    sum.weight <- rowSums((dat$nSample)*(dat$w.mat)^2,na.rm=TRUE);
    if(weight=='Npq+impQ') {
        w.mat <- sqrt(dat$nSample*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat;
        z.mat <- dat$ustat.mat/dat$vstat.mat;
        statistic.meta <- (rowSums(w.mat*z.mat,na.rm=TRUE))^2/rowSums(w.mat^2,na.rm=TRUE)
    }
    
    if(weight=='N+impQ') {
        w.mat <- sqrt(dat$nSample)*dat$w.mat;
        z.mat <- dat$ustat.mat/dat$vstat.mat;
        statistic.meta <- (rowSums(w.mat*z.mat,na.rm=TRUE))^2/rowSums(w.mat^2,na.rm=TRUE)
    }
    I2 <- NULL;cochranQ.pVal.mixChisq <- NULL;cochranQ.pVal <- NULL;cochranQ.df <- NULL;cochranQ.stat <- NULL;
    if(!is.null(extraPar$hetStat)) {
        if(extraPar$hetStat==TRUE) {
            beta.byStudy.mat <- dat$ustat.mat/(dat$vstat.mat)^2;
            for(ix.var in 1:nrow(beta.byStudy.mat)) {
                beta.byStudy <- beta.byStudy.mat[ix.var,];
                beta.var.byStudy <- 1/(dat$vstat.mat[ix.var,])^2
                w.mat <- matrix(0,nrow=length(beta.byStudy),ncol=length(beta.byStudy));
                diag(w.mat) <- 1;
                weight.byStudy <- (dat$vstat.mat[ix.var,])^2;
                weight.byStudy <- weight.byStudy/sum(weight.byStudy,na.rm=TRUE);
                w.mat <- w.mat+rm.na(matrix((-1)*rep(weight.byStudy,length(beta.byStudy)),nrow=length(beta.byStudy),ncol=length(beta.byStudy),byrow=TRUE));
                cochranQ.stat.mixChisq <- t(w.mat%*%rm.na(beta.byStudy))%*%(w.mat%*%rm.na(beta.byStudy))
                v.mat <- matrix(0,nrow=length(beta.byStudy),ncol=length(beta.byStudy));
                diag(v.mat) <- sqrt(rm.na(beta.var.byStudy));
                cochranQ.stat[ix.var] <- sum((beta.byStudy-beta.meta[ix.var])^2/(beta.var.byStudy+(beta.sd.meta[ix.var])^2-2*weight.byStudy*beta.var.byStudy),na.rm=TRUE);
                cochranQ.df[ix.var] <- length(which(!is.na(beta.byStudy-beta.meta[ix.var])^2))-1;
                if(cochranQ.df[ix.var]>0)
                {
                    cochranQ.pVal[ix.var] <- pchisq(cochranQ.stat[ix.var],df=cochranQ.df[ix.var],lower.tail=FALSE);
                    svd.mat <- try(svd(v.mat%*%w.mat%*%w.mat%*%v.mat),silent=TRUE);
                    cochranQ.pVal.mixChisq[ix.var] <- cochranQ.pVal[ix.var];
                    if(class(svd.mat)!="try-error") {
                        lambda <- svd.mat$d;
                        cochranQ.pVal.mixChisq[ix.var] <- try(liu(cochranQ.stat.mixChisq,lambda),silent=TRUE);
                    }
                }
                if(cochranQ.df[ix.var]<=0)
                {
                    cochranQ.pVal[ix.var] <- NA;
                    cochranQ.pVal.mixChisq[ix.var] <- NA;
                }
                
                I2[ix.var] <- (cochranQ.stat[ix.var]-cochranQ.df[ix.var])/cochranQ.stat[ix.var];
            }
            I2[which(I2<0)] <- 0;
            I2 <- paste(myFormat(I2*100,digits=2),"%",sep="");
        }
    }
    if(weight=="optim") statistic.meta <- ustat.meta^2/vstat.sq.meta;
    p.value.meta <- pchisq(statistic.meta,df=1,lower.tail=FALSE);
    nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
    af.meta <- rowSums((dat$nSample.mat)*(dat$af.mat),na.rm=TRUE)/rowSums((dat$nSample.mat),na.rm=TRUE);
    if(binaryTrait==FALSE) {
        beta.meta <- sign(ustat.meta)*sqrt(statistic.meta)*sqrt(var.y)/sqrt(2*nSample.meta*af.meta*(1-af.meta));
        beta.sd.meta <- sqrt(var.y/(2*nSample.meta*af.meta*(1-af.meta)));
    }
    if(binaryTrait==TRUE) {
        beta.meta <- sign(ustat.meta)*sqrt(statistic.meta)/sqrt(2*nSample.meta*af.meta*(1-af.meta))/sqrt(var.y);
        beta.sd.meta <- sqrt(1/sqrt(2*nSample.meta*af.meta*(1-af.meta)*var.y));
    }

    
    time.compAssoc <- Sys.time()-b;
    res.out <- cbind(gsub("_.*","",dat$pos),dat$ref.tab,dat$alt.tab,statistic.meta,p.value.meta,beta.meta,beta.sd.meta,nSample.meta,direction.meta);
    
    res.formatted <- cbind(gsub("_.*","",dat$pos),
                           dat$ref.tab,
                           dat$alt.tab,
                           format(af.meta,digits=3),
                           format(statistic.meta,digits=3),
                           format(p.value.meta,digits=3),
                           format(beta.meta,digits=3),
                           format(beta.sd.meta,digits=3),
                           nSample.meta,
                           direction.meta,
                           as.integer(sum.weight));
    colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION","EFFECTIVE_N");
    pos.only <- gsub("_.*","",dat.withMulti$pos);
    res.formatted.multi <- matrix(nrow=0,ncol=ncol(res.formatted));
    res.collapse <- matrix(nrow=0,ncol=2);
    if(rmMultiAllelicSite==FALSE) {
        for(ii in 1:length(dat$posMulti))  {
            ix.ii <- which(pos.only==dat$posMulti[ii]);
            cor.multi <- list();
            cor.multi$chrpos <- dat$posMulti[ii]
            if(length(unique(apply(matrix(dat.withMulti$ref.mat[ix.ii,],nrow=length(ix.ii)),1,uniq.allele)))==1) {
                cor.multi$ref <- paste(apply(matrix(dat.withMulti$ref.mat[ix.ii,],nrow=length(ix.ii)),1,uniq.allele),sep=',',collapse=',');
                cor.multi$alt <- paste(apply(matrix(dat.withMulti$alt.mat[ix.ii,],nrow=length(ix.ii)),1,uniq.allele),sep=',',collapse=',');
                
                af.meta <- rowSums((dat.withMulti$af.mat[ix.ii,])*(dat.withMulti$nSample.mat[ix.ii,]),na.rm=TRUE)/rowSums(dat.withMulti$nSample.mat[ix.ii,],na.rm=TRUE);
                cor.mat <- (-2)*as.vector(af.meta)%*%t(as.vector(af.meta));
                diag(cor.mat) <- 2*as.vector(af.meta)*(1-as.vector(af.meta));
                cor.mat <- rm.na(cov2cor(cor.mat));                
                diag(cor.mat) <- 1;
                cor.multi$cor <- paste(as.vector(cor.mat),sep=',',collapse=',')
                res.multi.ii <- (multiAlleleAssoc(dat$posMulti[ii],dat.withMulti,cor.multi));
                res.formatted <- rbind(res.formatted,res.multi.ii$res.formatted);
                res.collapse <- rbind(res.collapse,res.multi.ii$res.collapse);
            }
        }
    }
    
    
    return(list(res.out=res.out,
                res.formatted=res.formatted,
                res.collapse=res.collapse,
                I2=I2,
                cochranQ.pVal.mixChisq=cochranQ.pVal.mixChisq,
                cochranQ.pVal=cochranQ.pVal,
                cochranQ.df=cochranQ.df,
                cochranQ.stat=cochranQ.stat,
                formattedData=dat,
                formattedData.withMulti=dat.withMulti,
                posMulti=dat$posMulti,
                raw.data.all=raw.data.all,
                raw.imp.qual=raw.imp.qual,
                time.compAssoc=time.compAssoc,
                time.readData=time.readData));
}

#' Perform multi-allelic association tests;
#'
#' @param pos The position of multi-allelic site;
#' @param dat The formatted data;
#' @param corMultiAllele.mat the correlation matrix generated from multiple alleles;
#' @return a matrix with formatted result;
#' @export 
multiAlleleAssoc <- function(pos,dat,corMultiAllele.mat) {
    res.formatted <- matrix(rep(NA,11),nrow=1);
    res.collapse <- matrix(rep(NA,2),nrow=1);
    ix <- which(corMultiAllele.mat$chrpos==pos);
    if(length(ix)==1) {
        ref.cor <- unlist(strsplit(corMultiAllele.mat$ref[ix],split=','));
        alt.cor <- unlist(strsplit(corMultiAllele.mat$alt[ix],split=','));
        if(length(ref.cor)<length(alt.cor))  ref.cor[(length(ref.cor)+1):length(alt.cor)] <- ref.cor[1];
        if(length(ref.cor)>length(alt.cor))  alt.cor[(length(alt.cor)+1):length(ref.cor)] <- alt.cor[1];
        
        cor.mat <- rm.na(matrix(as.numeric(unlist(strsplit(corMultiAllele.mat$cor[ix],split=','))),nrow=length(ref.cor),ncol=length(ref.cor)));
        ref.alt.cor <- paste(ref.cor,alt.cor,sep=',');
        pos.dat <- gsub("_.*","",dat$pos);
        ix <- which(pos.dat==pos);
        dat$pos <- pos.dat[ix];
        dat$ref.mat <- matrix(dat$ref.mat[ix,],nrow=length(ix))
        dat$alt.mat <- matrix(dat$alt.mat[ix,],nrow=length(ix));
        dat$nref.mat <- matrix(dat$nref.mat[ix,],nrow=length(ix))
        dat$nalt.mat <- matrix(dat$nalt.mat[ix,],nrow=length(ix));
        dat$nhet.mat <- matrix(dat$nhet.mat[ix,],nrow=length(ix));
        dat$ustat.mat <- matrix(dat$ustat.mat[ix,],nrow=length(ix));
        dat$vstat.mat <- matrix(dat$vstat.mat[ix,],nrow=length(ix));
        dat$af.mat <- matrix(dat$af.mat[ix,],nrow=length(ix));
        dat$w.mat <- matrix(dat$w.mat[ix,],nrow=length(ix));
        dat$nSample.mat <- matrix(dat$nSample.mat[ix,],nrow=length(ix));
        dat$ref.alt.mat <- matrix(paste(dat$ref.mat,dat$alt.mat,sep=','),nrow=length(ix));
        ref.alt.dat <- unique(as.character(paste(dat$ref.mat,dat$alt.mat,sep=',')));
        ref.alt.both <- intersect(ref.alt.dat,ref.alt.cor);
        
        if(length(ref.alt.both)>0) {            
            ref.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$ref.mat));
            alt.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$alt.mat));
            ref.alt.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$ref.alt.mat));
            
            nref.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nref.mat));
            nalt.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nalt.mat));
            nhet.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nhet.mat));
            ustat.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$ustat.mat));
            vstat.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$vstat.mat));
            af.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$af.mat));
            w.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$w.mat));
            nSample.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nSample.mat));
            dat$pos <- pos.dat;
            
            for(ii in 1:ncol(ref.mat)) {
                ix.match <- match(ref.alt.both,dat$ref.alt.mat[,ii]);
                ref.mat[,ii] <- dat$ref.mat[ix.match,ii];
                alt.mat[,ii] <- dat$alt.mat[ix.match,ii];
                nref.mat[,ii] <- dat$nref.mat[ix.match,ii];
                nalt.mat[,ii] <- dat$nalt.mat[ix.match,ii];
                nhet.mat[,ii] <- dat$nhet.mat[ix.match,ii];
                ustat.mat[,ii] <- dat$ustat.mat[ix.match,ii];
                vstat.mat[,ii] <- dat$vstat.mat[ix.match,ii];
                af.mat[,ii] <- dat$af.mat[ix.match,ii];
                w.mat[,ii] <- dat$w.mat[ix.match,ii];
                nSample.mat[,ii] <- dat$nSample.mat[ix.match,ii];                
            }
            dat$ref.mat <- ref.mat;
            dat$alt.mat <- alt.mat;
            dat$nref.mat <- nref.mat;
            dat$nalt.mat <- nalt.mat;
            dat$nhet.mat <- nhet.mat;
            dat$ustat.mat <- ustat.mat;
            dat$vstat.mat <- vstat.mat;
            dat$af.mat <- af.mat;
            dat$w.mat <- w.mat;
            dat$nSample.mat <- nSample.mat;
            dat$direction.mat <- dat$ustat.mat;
            dat$direction.mat[which(dat$ustat.mat>0,arr.ind=T)] <- "+";
            dat$direction.mat[which(dat$ustat.mat<0,arr.ind=T)] <- "-";
            dat$direction.mat[which(dat$ustat.mat==0,arr.ind=T)] <- "=";
            dat$direction.mat[which(is.na(dat$ustat.mat),arr.ind=T)] <- "X";
            direction.meta <- apply(dat$direction.mat,1,paste,sep="",collapse="");
            
            af.meta <- rowSums((dat$af.mat)*(dat$nSample.mat),na.rm=T)/rowSums((dat$nSample.mat),na.rm=T);
            
            ix <- match(ref.alt.both,ref.alt.cor);
            if(length(ix)>0) {
                cor.mat <- matrix(as.numeric(cor.mat[ix,ix]),nrow=length(ix),ncol=length(ix));
                ustat.meta.tmp <- rowSums((dat$ustat.mat)*(dat$w.mat),na.rm=TRUE);
                V <- matrix(0,nrow=nrow(cor.mat),ncol=ncol(cor.mat));
                for(ii in 1:ncol(dat$ustat.mat)) {
                    S <- matrix(0,nrow=nrow(cor.mat),ncol=ncol(cor.mat));
                    diag(S) <- rm.na(sqrt((dat$w.mat[,ii])^2*(dat$vstat.mat[,ii])^2))
                    V <- V+S%*%cor.mat%*%S;
                }
                V <- regMat(V,0.1);

                beta.est.meta <- as.vector(ginv(V)%*%ustat.meta.tmp);
                beta.sd.meta <- as.vector(sqrt(diag(ginv(V))));

                statistic.meta <- (beta.est.meta/beta.sd.meta)^2;
                p.value.meta <- pchisq(statistic.meta,df=1,lower.tail=FALSE);
                nSample.meta <- rowSums(dat$nSample.mat,na.rm=T);
                
                ref.tab <- apply(dat$ref.mat,1,uniq.allele);
                alt.tab <- apply(dat$alt.mat,1,uniq.allele);
                res.formatted <- cbind(pos,
                                       ref.tab,
                                       alt.tab,
                                       format(af.meta,digits=3),
                                       format(statistic.meta,digits=3),
                                       format(p.value.meta,digits=3),
                                       format(beta.est.meta,digits=3),
                                       format(beta.sd.meta,digits=3),
                                       nSample.meta,
                                       direction.meta,
                                       nSample.meta);
                stat.collapse <- (sum(ustat.meta.tmp))^2/sum(V);
                p.value.collapse <- pchisq(stat.collapse,df=1,lower.tail=FALSE);
                res.collapse <- matrix(c(unique(pos),p.value.collapse),nrow=1);
            }
        }
    }
    return(list(res.formatted=res.formatted,
                res.collapse=res.collapse)); 
}

#' format data into matrices;
#' @param raw.data.all The read in from summary assoicaiton statistic files;
#' @param raw.imp.qual The raw imputation quality files;
#' @param impQualWeight If using imputation quality as weight
#' @return a list of converted data for downstream meta-analysis;
#' @export
GWAMA.formatData <- function(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual,...) {
    pos <- raw.data.all$pos;
    extraPar <- list(...);
    maf.cutoff <- extraPar$maf.cutoff;
    if(is.null(maf.cutoff)) maf.cutoff <- 1.01;
    knownVar <- extraPar$knownVar;
    ref.mat <- matrix(unlist(raw.data.all$ref),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    alt.mat <- matrix(unlist(raw.data.all$alt),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    
    ref.tab <- apply(ref.mat,1,uniq.allele);
    alt.tab <- apply(alt.mat,1,uniq.allele);
    nSample.mat <- matrix(unlist(raw.data.all$nSample),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    nref.mat <- matrix(unlist(raw.data.all$nref),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    nalt.mat <- matrix(unlist(raw.data.all$nalt),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    nhet.mat <- matrix(unlist(raw.data.all$nhet),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    af.mat <- matrix(unlist(raw.data.all$af),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    ustat.mat <- matrix(unlist(raw.data.all$ustat),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    vstat.mat <- matrix(unlist(raw.data.all$vstat),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    ref.mat <- matrix(unlist(raw.data.all$ref),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    alt.mat <- matrix(unlist(raw.data.all$alt),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));

    if(impQualWeight==FALSE) 
        w.mat <- matrix(1,nrow=nrow(vstat.mat),ncol=ncol(vstat.mat));
    if(!is.null(raw.imp.qual)) {
        imp.qual <- matrix(unlist(lapply(raw.imp.qual,getImpQual,pos=pos,col.impqual=col.impqual)),ncol=ncol(ref.mat),nrow=nrow(ref.mat));
        ix.lowQual <- which(imp.qual<impQual.lb,arr.ind=TRUE);
        ustat.mat[ix.lowQual] <- NA;
        vstat.mat[ix.lowQual] <- NA;
        w.mat <- sqrt(imp.qual);
        w.mat[which(is.na(w.mat),arr.ind=TRUE)] <- 1;
        nSample.mat[ix.lowQual] <- NA;
    }
    direction.mat <- ustat.mat;
    direction.mat[which(ustat.mat>0,arr.ind=T)] <- "+";
    direction.mat[which(ustat.mat<0,arr.ind=T)] <- "-";
    direction.mat[which(ustat.mat==0,arr.ind=T)] <- "=";
    direction.mat[which(is.na(ustat.mat),arr.ind=T)] <- "X";

    af.vec <- rowSums(af.mat*nSample.mat,na.rm=TRUE)/rowSums(nSample.mat,na.rm=TRUE);
    maf.vec <- rm.na(af.vec);
    maf.vec[which(maf.vec>.5)] <- 1-maf.vec[which(maf.vec>.5)];
    pos.norefalt <- gsub("_.*","",pos);
    ix.rare <- unique(c(which(maf.vec<=maf.cutoff & maf.vec>=0),which(pos.norefalt%in%knownVar)));
    nref.mat <- matrix(nref.mat[ix.rare,],nrow=length(ix.rare));
    nalt.mat <- matrix(nalt.mat[ix.rare,],nrow=length(ix.rare));
    ref.tab <- ref.tab[ix.rare];
    alt.tab <- alt.tab[ix.rare];
    af.mat <- matrix(af.mat[ix.rare,],nrow=length(ix.rare));
    ref.mat <- matrix(ref.mat[ix.rare,],nrow=length(ix.rare));
    alt.mat <- matrix(alt.mat[ix.rare,],nrow=length(ix.rare));
    nhet.mat <- matrix(nhet.mat[ix.rare,],nrow=length(ix.rare));
    ustat.mat <- matrix(ustat.mat[ix.rare,],nrow=length(ix.rare));
    w.mat <- matrix(w.mat[ix.rare,],nrow=length(ix.rare));
    pos <- pos[ix.rare]
    nSample.mat <- matrix(nSample.mat[ix.rare,],nrow=length(ix.rare));
    direction.mat <- matrix(direction.mat[ix.rare,],nrow=length(ix.rare));
    vstat.mat <- matrix(vstat.mat[ix.rare,],nrow=length(ix.rare));
    
    return(list(nref.mat=nref.mat,
                nalt.mat=nalt.mat,
                ref.tab=ref.tab,
                alt.tab=alt.tab,
                af.mat=af.mat,
                ref.mat=ref.mat,
                alt.mat=alt.mat,
                nhet.mat=nhet.mat,
                ustat.mat=ustat.mat,
                w.mat=w.mat,
                pos=pos,
                maf.meta=maf.vec[ix.rare],
                nSample.mat=nSample.mat,
                direction.mat=direction.mat,
                vstat.mat=vstat.mat));
}
#' remove multi-allelic sites;
#'
#' @param dat The dataset;
#' @return the data set without multi-allelic sites;
#' @export
GWAMA.rmMulti <- function(dat) {
    pos.tmp <- gsub("_.*","",dat$pos);
    pos.tab <- table(pos.tmp);
    pos.multi <- names(pos.tab)[which(pos.tab>1)];
    ix.multi <- (1:length(pos.tmp))[pos.tmp %in% pos.multi];
    if(length(ix.multi)>0) {
        dat$pos <- dat$pos[-ix.multi];
        dat$ref.mat <- matrix(dat$ref.mat[-ix.multi,],nrow=nrow(dat$ref.mat)-length(ix.multi),ncol=ncol(dat$ref.mat));
        dat$alt.mat <- matrix(dat$alt.mat[-ix.multi,],nrow=nrow(dat$alt.mat)-length(ix.multi),ncol=ncol(dat$alt.mat));
        dat$ref.tab <- dat$ref.tab[-ix.multi];
        dat$alt.tab <- dat$alt.tab[-ix.multi];
        dat$w.mat <- matrix(dat$w.mat[-ix.multi,],nrow=nrow(dat$w.mat)-length(ix.multi),ncol=ncol(dat$w.mat));
        dat$nSample.mat <- matrix(dat$nSample.mat[-ix.multi,],nrow=nrow(dat$nSample.mat)-length(ix.multi),ncol=ncol(dat$nSample.mat));
        dat$nref.mat <- matrix(dat$nref.mat[-ix.multi,],nrow=nrow(dat$nref.mat)-length(ix.multi),ncol=ncol(dat$nref.mat));
        dat$nalt.mat <- matrix(dat$nalt.mat[-ix.multi,],nrow=nrow(dat$nalt.mat)-length(ix.multi),ncol=ncol(dat$nalt.mat));
        dat$nhet.mat <- matrix(dat$nhet.mat[-ix.multi,],nrow=nrow(dat$nhet.mat)-length(ix.multi),ncol=ncol(dat$nhet.mat));
        dat$af.mat <- matrix(dat$af.mat[-ix.multi,],nrow=nrow(dat$af.mat)-length(ix.multi),ncol=ncol(dat$af.mat));
        dat$ustat.mat <- matrix(dat$ustat.mat[-ix.multi,],nrow=nrow(dat$ustat.mat)-length(ix.multi),ncol=ncol(dat$ustat.mat));
        dat$vstat.mat <- matrix(dat$vstat.mat[-ix.multi,],nrow=nrow(dat$vstat.mat)-length(ix.multi),ncol=ncol(dat$vstat.mat));
        dat$direction.mat <- matrix(dat$direction.mat[-ix.multi,],nrow=nrow(dat$direction.mat)-length(ix.multi),ncol=ncol(dat$direction.mat));
        dat$posMulti <- pos.multi;
    }
    return(list(dat=dat,
                posMulti=pos.multi));
}

#' gateway function for multi-allelic analysis in rareGWAMA
#'
#' @param score stat.file The score statistics file names;
#' @param imp.qual.file The imputation quality file names;
#' @param cor.multiallele.file The correlation information for multi-allelic sites;
#' @param tabix.range The tabix range for the variants to be analyzed;
#' @param alternative The alternative hypothesis; only the default choice two.sided is implemented now;
#' @param col.impqual The column number for the imputation quality files; The default choice is 5.
#' @param impQual.lb The lower bound for the imputation quality cutoff; Variants with imputation quality less than the cutoffs are labelled as missing and removed from the meta-analysis
#' @param impQualWeight Whether to apply the imputation quality based optimal weighting;
#' @export
#' @return A list of formatted output; 
rareGWAMA.single.multiAllele <- function(score.stat.file,imp.qual.file=NULL,tabix.range,alternative="two.sided",col.impqual=5,impQual.lb=0.7,impQualWeight=FALSE) {
    a <- Sys.time();
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range,multiAllelic = TRUE));
    raw.imp.qual <- NULL;
    if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
    
    time.readData <- Sys.time()-a;
    b <- Sys.time();
    raw.data.all <- raw.data.all[[1]];
    cat('Read in',length(raw.data.all$ref[[1]]),'variants',sep=' ');
    dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual);
    pos <- unique(gsub("_.*","",dat$pos));
    res.assoc.tmp <- sapply(pos,multiAlleleAssoc,dat=dat,corMultiAllele.mat=corMultiAllele.mat);
    res.formatted <- do.call(rbind,res.assoc.tmp);
    colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION");
    return(list(res.formatted=res.formatted,
                formattedData=dat,
                raw.data.all=raw.data.all,
                raw.imp.qual=raw.imp.qual));
}
#' convert chisq statistic to beta for binary trait assuming the variance for y;
#' @param statistic chisq stat;
#' @param beta.ori original beta; for getting the sign;
#' @param var.y the variance for Y. for binary trait, the variance of y is frac.case*(1-frac.case); for continuous trait, the variance for y is typically set to 1; 
#' @param af allele freq;
#' @param N sample size;
#' @export
convertChisq2Beta <- function(statistic,beta.ori,var.y,af,N,binaryTrait=FALSE) {

    if(binaryTrait==FALSE) {
        beta.out <- sign(beta.ori)*sqrt(statistic)*sqrt(var.y)/sqrt(2*N*af*(1-af));
        beta.sd <- sqrt(var.y/(2*N*af*(1-af)));
    }
    if(binaryTrait==TRUE) {
        beta.out <- sign(beta.ori)*sqrt(statistic)/sqrt(2*N*af*(1-af))/sqrt(var.y);
        beta.sd <- sqrt(1/sqrt(2*N*af*(1-af)*var.y));
    }
    return(list(beta.est=beta.out,
                beta.sd=beta.sd));
}
