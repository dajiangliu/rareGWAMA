#' Safe sum
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
Sum <- function(a)
  {
    return(sum(a,na.rm=TRUE));
  }

#' Safe mean
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
Mean <- function(a)
  {
    return(mean(a,na.rm=TRUE));
  }
    
#' Safe colSums
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
ColSums<- function(a)
  {
    return(colSums(a,na.rm=TRUE));
  }
    
#' Safe rowSums
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
RowSums<- function(a)
  {
    return(rowSums(a,na.rm=TRUE));
  }


#' power calculation
#'
#' @param beta genetic effect sizes;
#' @param maf minor allele frequencies;
#' @param alpha significance threshold
#' @param N sample size
#' @return power
#' @export
pwrCalc <- function(beta,maf,alpha,N=300000)
    {
        z.alpha <- qnorm(1-alpha/2);
        se.beta <- 1/sqrt(2*(maf)*(1-maf)*N);
        pwr <- 1-pnorm(z.alpha-beta/se.beta)+pnorm(-z.alpha-beta/se.beta)
        return(pwr);
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
