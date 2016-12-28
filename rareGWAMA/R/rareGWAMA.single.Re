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
    ix <- match(pos,pos.imp);
    return(as.numeric(imp.qual[ix,col.impqual]));
}
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
rareGWAMA.single <- function(score.stat.file,imp.qual.file=NULL,tabix.range,alternative="two.sided",col.impqual=5,impQual.lb=0.7,impQualWeight=FALSE,rmMultiAllelicSite=TRUE) {
    a <- Sys.time();
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range));
    raw.imp.qual <- NULL;
    if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
    time.readData <- Sys.time()-a;
    b <- Sys.time();
    raw.data.all <- raw.data.all[[1]];
    cat('Read in',length(raw.data.all$ref[[1]]),'variants',sep=' ');
    dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual);
    if(rmMultiAllelicSite==TRUE) {
        tmp <- GWAMA.rmMulti(dat);
        dat <- tmp$dat;posMulti <- tmp$posMulti;
    }
    direction.meta <- apply(dat$direction.mat,1,paste,sep="",collapse="");   
    ustat.meta <- rowSums((dat$ustat.mat)*(dat$w.mat),na.rm=TRUE);
    vstat.sq.meta <- rowSums((dat$w.mat)^2*(dat$vstat.mat)^2,na.rm=TRUE);
    beta.meta <- (ustat.meta)/(vstat.sq.meta);
    beta.sd.meta <- sqrt(1/vstat.sq.meta);
    statistic.meta <- ustat.meta^2/vstat.sq.meta;
    p.value.meta <- pchisq(statistic.meta,df=1,lower.tail=FALSE);
    nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
    af.meta <- rowSums((dat$nSample.mat)*(dat$af.mat),na.rm=TRUE)/rowSums((dat$nSample.mat),na.rm=TRUE);
    time.compAssoc <- Sys.time()-b;
    res.out <- cbind(dat$pos,dat$ref.tab,dat$alt.tab,statistic.meta,p.value.meta,beta.meta,beta.sd.meta,nSample.meta,direction.meta);    
    res.formatted <- cbind(dat$pos,
                           dat$ref.tab,
                           dat$alt.tab,
                           format(af.meta,digits=3),
                           format(statistic.meta,digits=3),
                           format(p.value.meta,digits=3),
                           format(beta.meta,digits=3),
                           format(beta.sd.meta,digits=3),
                           nSample.meta,
                           direction.meta);
    colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION");
    return(list(res.out=res.out,
                res.formatted=res.formatted,
                formattedData=dat,
                posMulti=dat$posMulti,
                raw.data.all=raw.data.all,
                raw.imp.qual=raw.imp.qual,
                time.compAssoc=time.compAssoc,
                time.readData=time.readData));
}
#' format data into matrices;
#' @param raw.data.all The read in from summary assoicaiton statistic files;
#' @param raw.imp.qual The raw imputation quality files;
#' @param impQualWeight If using imputation quality as weight
#' @return a list of converted data for downstream meta-analysis;
#' @export
GWAMA.formatData <- function(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual) {
    pos <- raw.data.all$pos;
    ref.mat <- matrix(unlist(raw.data.all$ref),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    alt.mat <- matrix(unlist(raw.data.all$alt),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
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
        if(impQualWeight==TRUE) 
            w.mat <- sqrt(imp.qual);
        nSample.mat[ix.lowQual] <- NA;
    }
    direction.mat <- ustat.mat;
    direction.mat[which(ustat.mat>0,arr.ind=T)] <- "+";
    direction.mat[which(ustat.mat<0,arr.ind=T)] <- "-";
    direction.mat[which(ustat.mat==0,arr.ind=T)] <- "=";
    direction.mat[which(is.na(ustat.mat),arr.ind=T)] <- "X";
       
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
    pos.multi <- unique(gsub('/[0-9]+','',grep("/",dat$pos,value=TRUE)));
    ix.multi <- unique(c(grep("/",dat$pos),match(pos.multi,dat$pos)));
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
rareGWAMA.single.multiAllele <- function(score.stat.file,imp.qual.file=NULL,cor.multiallele.file,tabix.range,alternative="two.sided",col.impqual=5,impQual.lb=0.7,impQualWeight=FALSE) {
    a <- Sys.time();
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range));
    raw.imp.qual <- NULL;
    if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
    corMultiAllele.mat <- read.table(cor.multiallele.file,header=TRUE,as.is=TRUE);
    corMultiAllele.mat$chrpos <- paste(corMultiAllele.mat$chr,corMultiAllele.mat$pos,sep=":");
    
    time.readData <- Sys.time()-a;
    b <- Sys.time();
    raw.data.all <- raw.data.all[[1]];
    cat('Read in',length(raw.data.all$ref[[1]]),'variants',sep=' ');
    dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual);
    pos <- unique(gsub("/[0-9]+","",dat$pos));
    res.assoc.tmp <- sapply(pos,multiAlleleAssoc,dat=dat,corMultiAllele.mat=corMultiAllele.mat);
    res.formatted <- do.call(rbind,res.assoc.tmp);
    colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION");
    return(list(res.formatted=res.formatted,
                formattedData=dat,
                raw.data.all=raw.data.all,
                raw.imp.qual=raw.imp.qual));
}

#' Perform multi-allelic association tests;
#'
#' @param pos The position of multi-allelic site;
#' @param dat The formatted data;
#' @param corMultiAllele.mat the correlation matrix generated from multiple alleles;
#' @return a matrix with formatted result;
#' @export 
multiAlleleAssoc <- function(pos,dat,corMultiAllele.mat) {
    res.formatted <- matrix(rep(NA,10),nrow=1);
    ix <- which(corMultiAllele.mat$chrpos==pos);
    if(length(ix)==1) {
        ref.cor <- unlist(strsplit(corMultiAllele.mat$ref[ix],split=','));
        alt.cor <- unlist(strsplit(corMultiAllele.mat$alt[ix],split=','));
        if(length(ref.cor)<length(alt.cor))  ref.cor[(length(ref.cor)+1):length(alt.cor)] <- ref.cor[1];
        if(length(ref.cor)>length(alt.cor))  alt.cor[(length(alt.cor)+1):length(ref.cor)] <- alt.cor[1];
        
        cor.mat <- rm.na(matrix(as.numeric(unlist(strsplit(corMultiAllele.mat$cor[ix],split=','))),nrow=length(ref.cor),ncol=length(ref.cor)));
        ref.alt.cor <- paste(ref.cor,alt.cor,sep=',');
        pos.dat <- gsub("/[0-9]+","",dat$pos);
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
        
        matchAllele <- function(ref.alt.dat,ref.alt.both) return(match(ref.alt.dat,ref.alt.both))
        if(length(ref.alt.both)>0) {
            ix.allele.match <- apply(dat$ref.alt.mat,1,matchAllele,ref.alt.both=ref.alt.both);
            
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
            ix.allele.match <- as.matrix(ix.allele.match);
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
                vstat.sq.meta.tmp <- rowSums((dat$w.mat)^2*(dat$vstat.mat)^2,na.rm=TRUE);
                S <- matrix(0,nrow=nrow(cor.mat),ncol=ncol(cor.mat));
                diag(S) <- sqrt(vstat.sq.meta.tmp);
                diag(cor.mat) <- 0;
                cor.mat[which(abs(cor.mat)>0.9,arr.ind=TRUE)] <- 0;
                diag(cor.mat) <- 1;
                cor.mat <- regMat(cor.mat,0.1);
                V <- S%*%cor.mat%*%S
                beta.est.meta <- as.vector(ginv(V)%*%ustat.meta.tmp);
                beta.sd.meta <- as.vector(sqrt(diag(ginv(V))));
                statistic.meta <- (beta.est.meta/beta.sd.meta)^2;
                p.value.meta <- pchisq(statistic.meta,df=1,lower.tail=FALSE);
                nSample.meta <- rowSums(dat$nSample.mat,na.rm=T);
                uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
                
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
                                       direction.meta);
            }
        }
    }
    return(list(res.formatted=res.formatted)); 
}
