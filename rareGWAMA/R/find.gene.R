#' Find the chromosomal position for genes;
#'
#' @param candidate.gene Name of candidate genes;
#' @return a list consisiting of tabix, pos.start and pos.end;
#' @export
find.gene.chrpos <- function(candidate.gene)
  {
    data(refFlat.resource);
    res <- res.refFlat;
    gene <- res$gene;
    chr <- res$chr;
    pos.start <- res$pos.start;
    pos.end <- res$pos.end;
    ix <- match(candidate.gene,gene);
    tabix <- paste(chr[ix],paste(pos.start[ix],pos.end[ix],sep="-"),sep=":");
    return(list(tabix=tabix,
                pos.start=paste(chr[ix],pos.start[ix],sep=":",collapse=":"),
                pos.end=paste(chr[ix],pos.end[ix],sep=":",collapse=":")));
  }
#' change position to tabix range;
#'
#' @param variant.vec The vector of variants of the format "1:12345";
#' @return A vector of tabix ranges;
#' @export
get.tabix.range <- function(variant.vec)
  {
    variant.list <- variant.vec;
    range.vec <- 0;
    variant.list.split <- matrix(unlist(strsplit(variant.list,split=":")),ncol=2,byrow=TRUE);
    variant.list.split[,2] <- paste(variant.list.split[,2],variant.list.split[,2],sep='-');
    range.vec <- paste(variant.list.split[,1],variant.list.split[,2],sep=':',collapse=',');
    
    return(range.vec); 
  }

#' Find nearby genes for a sigificant SNP
#'
#' @param known.variant.vec Positions of known variants, in the format of 1:12345
#' @param window.size The window size. Genes within the window are extracted. Default size is 1e6
#' @return A vector of genes;
#' @export
get.gene.inWindow <- function(known.variant.vec,window.size=1e6)
  {
    data(refFlat.resource);
    res <- res.refFlat;
    gene <- res$gene;
    chr <- res$chr;
    pos.start <- res$pos.start;
    pos.end <- res$pos.end;
    variant.list.split <- matrix(unlist(strsplit(known.variant.vec,split=":")),ncol=2,byrow=TRUE);
    pos.variant.list <- as.integer(variant.list.split[,2]);
    chr.variant.list <- as.integer(variant.list.split[,1]);
    pos.window.start <- max(pos.variant.list-window.size);
    pos.window.end <- min(pos.variant.list+window.size);
    ix.gene <- which(chr==chr.variant.list[1] & pos.start>pos.window.start & pos.end<pos.window.end);
    if(length(ix.gene)>0) {
        pos.start.end <- cbind(pos.start[ix.gene],pos.end[ix.gene]);
        pos.median <- apply(pos.start.end,1,median);
        dist.vec <- abs(pos.variant.list-pos.median);
        ix.sort <- sort(dist.vec,index.return=T)$ix;
        ix.gene <- ix.gene[ix.sort];
    }
    gene.vec <- gene[ix.gene];
    return(gene.vec);
  }

#' Find top variants in a gene region that satisfy a MAF cutoff;
#'
#' @param candidate.vairant.vec The vector of candidate variant;
#' @param window.size The size of window, default is 1e6
#' @param singlevar.result Single varinat association results;
#' @param pval.cutoff The cutoffs of p-values;
#' @export
find.top.variant <- function(candidate.variant.vec,window.size,singlevar.result,pval.cutoff=3.1e-7)
  {
      chr.pos.tmp <- matrix(unlist(strsplit(candidate.variant.vec,split=":")),byrow=TRUE,ncol=2);
    
    chr.candidate <- as.numeric(chr.pos.tmp[,1]);
    pos.candidate <- as.numeric(chr.pos.tmp[,2]);
    
    chrpos.result <- matrix(unlist(strsplit(singlevar.result$POS,split=":")),byrow=TRUE,ncol=2);
    chr.result <- as.numeric(chrpos.result[,1]);
    pos.result <- as.numeric(chrpos.result[,2]);
    
    p.value.vec <- 0;known.variant.collapsed.vec <- 0;candidate.variant.tabix.vec <- 0;anno.known.variant.vec <- 0;anno.candidate.variant.vec <- 0;ix.rm <- integer(0);
    res <- list();
    for(ii in 1:length(candidate.variant.vec))
      {
        
        pos.start.window <- max(pos.candidate[ii]-window.size,1);
        pos.end.window <- pos.candidate[ii]+window.size;
        ix.candidate <- which(singlevar.result$POS==candidate.variant.vec[ii]);
        ix.var <- which(chr.result==chr.candidate[ii] & pos.result>=pos.start.window & pos.result<=pos.end.window);
        res[[ii]] <- list();
        candidate.variant <- candidate.variant.vec[ii];
        if(length(ix.var)==0) ix.rm <- c(ix.rm,ii)
        if(length(ix.var)>0){
          p.value.tmp <- singlevar.result$PVALUE[ix.var];
          ix.top <- which(p.value.tmp<pval.cutoff);
          if(length(ix.top)>0)
            {
              ix.known <- ix.var[ix.top];
              chrpos.known <- singlevar.result$POS[ix.known];
              candidate.in.known <- which(chrpos.known==candidate.variant);
              pval.known <- singlevar.result$PVALUE[ix.known];
              anno.known <- singlevar.result$ANNO[ix.known];              
              if(length(candidate.in.known)>0) {
                chrpos.known <- chrpos.known[-candidate.in.known];
                pval.known <- pval.known[-candidate.in.known];
                anno.known <- anno.known[-candidate.in.known];
                ix.known <- ix.known[-candidate.in.known];
              }
              res[[ii]] <- list(candidate.variant=candidate.variant.vec[ii],
                                ix.candidate=ix.candidate,
                                ix.var=ix.var,
                                ix.known=ix.known);
            }
    
        }
      }
    return(res)
  }

