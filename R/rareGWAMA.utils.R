#' format according to publication standard;
#'
#' @param n The number to be formatted;
#' @param digits how many digits to be retained;
#' @return formatted result a vector;
#' @export 
myFormat <- function(n,digits) {
    res <- rep(NA,length(n));
    if(length(which(abs(n)>1))>0)
        res[which(abs(n)>1)] <- trimws(as.character(round(n[which(abs(n)>1)],digits=digits)));
    if(length(which(abs(n)>0.0001 & abs(n) <=1)))
        res[which(abs(n)>0.0001 & abs(n) <=1 )] <- trimws(as.character(apply(as.matrix(n[which(abs(n)>0.0001 & abs(n) <=1 )]),1,format,digits=digits,scientific=FALSE)))
    if(length(which(abs(n)<=0.0001))>0)
        res[which(abs(n)<=0.0001)] <- trimws(as.character(apply(as.matrix(n[which(abs(n)<=0.0001)]),1,format,digits=digits,scientific=TRUE)));
    return(res);
}
