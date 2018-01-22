#' perform single tissue twas and multiple tissue twas;
#'
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param vcf.ref.file the file names of the reference panel file;
#' @param transcript The list of transcript for twas;
#' @return formatted data and assocation testing results; 
#' @export
rareGWAMA.TWAS <- function(score.stat.file,imp.qual.file,vcf.ref.file,transcript,...) {
    extraPar <- list(...);
    weightRef <- extraPar$weightRef;
    if(is.null(weightRef)) weightRef <- "1kg";
    load.r2 <- extraPar$load.r2;
    if(is.null(load.r2)) load.r2 <- FALSE;
    fname.r2.prefix <- extraPar$fname.r2.prefix;
    fname.r2.postfix <- extraPar$fname.r2.postfix;
    if(is.null(fname.r2.prefix)) fname.r2.prefix <- "./";
    if(is.null(fname.r2.postfix)) fname.r2.postfix <- "_hrc.RData";
        
    refGeno <- extraPar$refGeno;
    col.impqual <- extraPar$col.impqual;
    impQual.lb <- extraPar$impQual.lb;
    no.per.batch <- extraPar$no.per.batch;
    if(is.null(no.per.batch)) no.per.batch <- 1000;;

    impQualWeight <- FALSE;
    rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
    if(is.null(col.impqual)) col.impqual <- 5;
    if(is.null(impQual.lb)) impQual.lb <- 0.7;
    if(is.null(rmMultiAllelicSite)) rmMultiAllelicSite <- TRUE;
    if(is.null(refGeno)) refGeno <- "GT";
    refGeno <- extraPar$refGeno;
    col.impqual <- extraPar$col.impqual;
    impQual.lb <- extraPar$impQual.lb;
    impQualWeight <- FALSE;
    rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
    if(is.null(col.impqual)) col.impqual <- 5;
    if(is.null(impQual.lb)) impQual.lb <- 0.7;
    if(is.null(rmMultiAllelicSite)) rmMultiAllelicSite <- TRUE;
    if(is.null(refGeno)) refGeno <- "GT";
    transcript.all <- transcript;
    
    counter <- 0;        pval.vec <- 0;statistic.vec <- 0;transcript.out <- 0;pval.t.vec <- 0;
    res.twas <- NULL;dat.twas <- NULL;
    dat <- NULL;
    for(aa in 1:max(1,as.integer(length(transcript.all)/no.per.batch))) {
        
        aa.start <- (aa-1)*no.per.batch+1;
        aa.end <- aa*no.per.batch;

        if(aa>=as.integer(length(transcript.all)/no.per.batch)) aa.end <- length(transcript.all);
        transcript <- transcript.all[aa.start:aa.end];
        ix.match <- match(transcript,transcript.list[,1]);
        
        transcript.sub <- matrix(transcript.list[ix.match,],ncol=2);
        eqtl.chrpos <- character(0);

        chrpos <- character(0);
        transcript.vec <- transcript.sub[,1];
        chr.vec <- transcript.sub[,2];
        
        chrpos <- character(0);
        for(jj in 1:length(twas.weight)) {
            ix <- which(twas.weight[[jj]]$gene%in%transcript.vec);
            chrpos <- c(chrpos,twas.weight[[jj]]$chrpos[ix]);
        }
        chrpos <- unique(chrpos);
        cat('Reading ',length(chrpos),' eQTL SNPs\n',sep=' ');
        
        a <- Sys.time();
        tabix.range <- get.tabix.range(chrpos);
        
        capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range,multiAllelic = TRUE)[[1]]);
        raw.imp.qual <- NULL;
        if(!is.null(imp.qual.file))
            raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);    
        cat('Read in',length(raw.data.all$ref[[1]]),'eQTL SNPs\n',sep=' ');
        dat.all <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual);
        dat.all$pos <- gsub("_.*","",dat.all$pos);
        
        
        for(ii in 1:length(transcript.vec)) {
            transcript.ii <- transcript.vec[ii];
            chrpos <- character(0);
            for(jj in 1:length(twas.weight)) {
                ix <- which(twas.weight[[jj]]$gene==transcript.ii);
                chrpos <- c(chrpos,twas.weight[[jj]]$chrpos[ix]);
            }
            chrpos <- unique(chrpos);
            if(length(chrpos)>0) {
                tabix.range <- get.tabix.range(chrpos);
                a <- Sys.time();
                cat('Analyzing ',transcript.ii,' ');
                chr.ii <- table(gsub(":.*","",chrpos));
                chr.ii <- as.numeric(names(chr.ii[which.max(chr.ii)]));

                if(load.r2==FALSE) 
                    r2 <- calc.r2(vcf.ref.file[as.numeric(chr.vec[ii])],tabix.range,refGeno);
                if(load.r2==TRUE) {
                    fname.r2 <- paste0(fname.r2.prefix,transcript.ii,fname.r2.postfix);
                    load(fname.r2);
                }
                if(is.null(r2)) {
                    warning.msg <- paste0(transcript.ii,' has no r2 information, skip!');
                    warning(warning.msg);
                }
                if(!is.null(r2)) {
                    pos.r2 <- colnames(r2);
                    ix <- dat.all$pos%in%chrpos;
                    dat <- list(pos=dat.all$pos[ix],
                                ustat.mat=as.matrix(dat.all$ustat.mat[ix,],nrow=length(ix)),
                                vstat.mat=as.matrix(dat.all$vstat.mat[ix,],nrow=length(ix)),
                                nSample.mat=dat.all$nSample.mat
                                );
                    
                    ix.match <- match(dat$pos,pos.r2);
                    r2 <- matrix(r2[ix.match,ix.match],nrow=length(ix.match));
                    diag(r2) <- 1;
                    r2 <- rm.na(r2);
                    r2 <- regMat(r2,.1);
                    counter <- counter+1;
                    
                    cov.exp <- matrix(0,nrow=length(twas.weight),ncol=length(twas.weight));
                    exp.vec <- rep(0,length(twas.weight));
                    weight.mat <- matrix(0,ncol=length(twas.weight),nrow=length(dat$pos));
                    for(kk in 1:length(twas.weight)) {
                        pos.eqtl <- intersect(twas.weight[[kk]]$chrpos[which(twas.weight[[kk]]$gene==transcript.ii)],dat$pos);
                        ix.eqtl <- match(pos.eqtl,dat$pos);
                        weight.mat[ix.eqtl,kk] <- twas.weight[[kk]]$weight[match(pos.eqtl,twas.weight[[kk]]$chrpos)];
                        
                    }
                    exp.vec <- rowSums(t(weight.mat)%*%rm.na(dat$ustat.mat),na.rm=TRUE);
                    for(ix.study in 1:ncol(dat$vstat.mat)) {
                        cov.mat <- t(r2*rm.na(dat$vstat.mat[,ix.study]))*rm.na(dat$vstat.mat[,ix.study]);
                        
                        cov.exp <- cov.exp+t(weight.mat)%*%cov.mat%*%weight.mat;
                    }
                    ix.rm <- which(exp.vec==0);
                    if(length(ix.rm)>0) {
                        exp.vec <- exp.vec[-ix.rm];
                        cov.exp <- as.matrix(cov.exp[-ix.rm,-ix.rm]);
                    }
                    dat.twas <- list();
                    dat.twas$ustat.meta <- exp.vec;
                    dat.twas$V.meta <- cov.exp;
                    res.twas <- rareGWAMA.skat(dat.twas,weight='linear');
                    
                    pval.vec[counter] <- res.twas$p.value;
                    res.twas <- rareGWAMA.t(dat.twas);
                    pval.t.vec[counter] <- res.twas$p.value;
                    statistic.vec[counter] <- res.twas$statistic;
                    transcript.out[counter] <- transcript.ii;
                    cat('Time usage ', Sys.time()-a,' \n');
                }
            }
        }
    }
    res.out <- cbind(transcript.out,statistic.vec,pval.vec,pval.t.vec);
    colnames(res.out) <- c("GENE","STATISTIC","PVALUE","PVALUE.t");

    return(list(res.out=res.out,
                dat.twas=dat.twas,
                res.twas=res.twas,
                dat=dat,
                dat.all=dat.all,
                r2=r2,
                weight.mat=weight.mat,
                raw.data.all=raw.data.all
                ));
}


#'
#' calculate R2 from a list of variants;
#' @param vcf.ref.file VCF reference file
#' @param tabix.range tabix range for a list of variants;
#' @param refGeno the FORMAT field for genotypes; could be GT or DS;
#' @return calcualted r2 matrix; with positions being the row and colnames;
#' @export
calc.r2 <- function(vcf.ref.file,tabix.range,refGeno="GT") {
    vcfIndv <- refGeno;
    annoType <- "";
    vcfColumn <- c("CHROM","POS","REF","ALT");
    vcfInfo <- NULL;      
    geno.list <- readVCFToListByRange(vcf.ref.file, tabix.range, "", vcfColumn, vcfInfo, vcfIndv)      
    if(refGeno=="DS") {
        gt <- geno.list$DS;
        if(is.null(gt)) {
            return(NULL);
        }
        gt <- matrix(as.numeric(gt),nrow=nrow(gt),ncol=ncol(gt));
        
        
    }
    if(refGeno=="GT") {
        gt.tmp <- geno.list$GT;
        if(is.null(gt.tmp)) return(NULL);
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
    r2 <- cor(gt,use='pairwise.complete');
    r2 <- rm.na(r2);
    pos.vcf <- paste(geno.list$CHROM,geno.list$POS,sep=":");
    colnames(r2) <- pos.vcf;
    rownames(r2) <- pos.vcf;
    diag(r2) <- 1;
    return(r2);
}
