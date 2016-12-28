#' Estimate the variance explained for a list of variants (possibly in LD);
#'
#' @param score.stat.file files of score statistics
#' @param cov.file Covariance matrix files
#' @param range tabix range for each gene/region; Tabix range needs to be in the format of "1:123-1234". Quotation marks are necessary. 
#' @param range.name The name of the range,e.g. gene names can be used;
#' @param refaltList A list of reference, alternative allele, a vector of alternative allele frequencies; Specifically, the list should consist of 
#' @param maf.cutoff MAF cutoff used to analyze variants; Default value is 1, i.e. no cutoffs are applied. MAFs are based upon the sample MAFs. 
#' @param out.digits Number of digits used in the output
#' @param callrate.cutoff Cutoffs of call rate, lower than which will NOT be analyzed (labelled as missing)
#' @param hwe.cutoff Cutoffs of HWE p-values
#' @param correctFlip Correct for flipped sites for score statistics and their covariance matrices
#' @param lambda Tuning parameter to the variant the variance correlation matrix to avoid near singular LD matrix; Assume the orignal LD matrix is R2, R2.adj=R2+lambda*Id;
#' @return a list consisting of results;
#' @export
rareMETALS.estimatePVE <- function(score.stat.file,cov.file,range,range.name,refaltList,maf.cutoff,out.digits,callrate.cutoff,hwe.cutoff,correctFlip=TRUE,lambda=0)
{

    res.meta <- rareMETALS.range.group(score.stat.file,cov.file,range,range.name,test='GRANVIL',refaltList,maf.cutoff,alternative='two.sided',out.digits,callrate.cutoff,hwe.cutoff,NULL,correctFlip,TRUE);
    beta.joint <- NA;
    locusPVE <- NA;
    if(length(res.meta$integratedData[[1]])>0) {
        datLocus <- (res.meta$integratedData[[1]]);
        corGLocus <- (res.meta$res.list[[1]])$r2.approx;
        covXLocus <- rm.na(datLocus$cov[[1]]);
        
        varLocusSNP <- 0;
        posLocus <- datLocus$pos;
        
        afLocus <- rm.na(datLocus$af[[1]]);
        
        ustatLocus <- rm.na(datLocus$ustat[[1]]);
        vstatLocus <- rm.na(datLocus$vstat[[1]]);
        ix.rm <- which(colSums(is.na(datLocus$cov[[1]]))==nrow(datLocus$cov[[1]]));
            
        if(length(ix.rm)>0) {
            ustatLocus <- ustatLocus[-ix.rm];
            vstatLocus <- vstatLocus[-ix.rm];
            afLocus <- afLocus[-ix.rm];
            covXLocus <- matrix(covXLocus[-ix.rm,-ix.rm],nrow=length(ustatLocus),ncol=length(ustatLocus));
            corGLocus <- matrix(corGLocus[-ix.rm,-ix.rm],nrow=length(ustatLocus),ncol=length(ustatLocus));
            VLocus <- matrix(VLocus[-ix.rm,-ix.rm],nrow=length(ustatLocus),ncol=length(ustatLocus));    
            posLocus <- posLocus[-ix.rm];
        }
        posLocus <- paste(posLocus,sep=",",collapse=",");
        
        Id <- matrix(0,nrow=nrow(covXLocus),ncol=ncol(covXLocus));
        diag(Id) <- 1;
        
        M <- matrix(0,nrow=length(afLocus),ncol=length(afLocus));
        diag(M) <- sqrt(2*afLocus*(1-afLocus));
        Id <- matrix(0,nrow=nrow(corGLocus),ncol=ncol(corGLocus));
        diag(Id) <- 1;
        
        corGLocus <- corGLocus+lambda*Id;
        
        covGLocus <- M%*%corGLocus%*%M;
        v0 <- matrix(0,nrow=length(ustatLocus),ncol=length(ustatLocus));
        diag(v0) <- vstatLocus;
        beta.joint <- ginv(v0%*%corGLocus%*%v0)%*%ustatLocus;
        
        locusPVE <- t(beta.joint)%*%covGLocus%*%beta.joint;
        
    }
    return(list(beta.joint=beta.joint,
                locusPos=posLocus,
                locusPVE=locusPVE));
    
    
}
