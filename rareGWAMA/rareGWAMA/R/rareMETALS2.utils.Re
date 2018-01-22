#' This function quickly paste chrom and pos;
#'
#' @param chrom Chromosome number;
#' @param pos Position
#' @return chrom:pos;
#' @export
#'
paste.pos <- function(chrom,pos) {
    return(paste(chrom,pos,sep=":"));
}


#' This function implement how to split a vector into chunks of certain size
#'
#' @param vecIn The input vector of any type;
#' @param chunkSize An integer for the chunk.size;
#' @return A list with splitted chunks;
#' @export
splitVecByChunkSize <- function(vecIn,chunkSize) {
    vec.split <- split(vecIn, ceiling(seq_along(vecIn)/chunkSize));
    return(vec.split);
}


#' This function implements function to colCenter a matrix so that colSums(matIn.centered) is 0;
#'
#' @param matIn Matrix input;
#' @return colCentered matrix;
#' @export
colCenter <- function(matIn) {
    matIn <- as.matrix(matIn);
    matOut <- matIn-matrix(rep(colMeans(matIn,na.rm=TRUE),nrow(matIn)),nrow=nrow(matIn),ncol=ncol(matIn),byrow=T);
    return(matOut);
}

#' get binomial confidence intervals via normal approx
#'
#' @param n sample size
#' @param p estimated proportions
#' @param alpha alpha level
#' @return conf.upper and conf.lower
#' @export
binom.conf <- function(n,p,alpha=0.95) {
    b <- sqrt(p*(1-p)/n);
    z <- sqrt(qchisq(alpha,df=1));
    conf.lower <- p-z*b;
    conf.upper <- p+z*b;
    return(list(conf.lower=conf.lower,
                conf.upper=conf.upper));
}


#' To minor: change number >.5 to 1-#
#'
#' @param vecIn
#' @return a vector with flipped number;
#' @export
toMinor <- function(vecIn) {
    if(is.vector(vecIn)) {
        ix.flip <- which(vecIn>0.5);
        if(length(ix.flip)>0) {
            vecIn[ix.flip] <- 1-vecIn[ix.flip];
        }
    }
    if(is.matrix(vecIn)) {
        ix.flip <- which(vecIn>0.5,arr.ind=TRUE);
        if(length(ix.flip)>0) {
            vecIn[ix.flip] <- 1-vecIn[ix.flip];
        }
    }
        
    return(vecIn);
}
