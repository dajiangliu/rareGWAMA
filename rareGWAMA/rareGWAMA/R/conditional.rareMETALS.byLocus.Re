conditional.rareMETALS.byLocus <- function(ANNO,gene.name,score.stat.file,cov.file,singlevar.result,test='GRANVIL',maf.cutoff,no.boot=10000,alternative=c('two.sided','greater','less'),alpha=0.05,ix.gold=1,out.digits=4,window.size=1e6,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",topN=1)
  {
    #singlevar.result consists of 3 columns: chr, pos, and p-value;
    extra.pars <- list(alpha=alpha,ix.gold=ix.gold,out.digits=out.digits,window.size=window.size,QC.par=list(callrate.cutoff=callrate.cutoff,hwe.cutoff=hwe.cutoff));
    window.size <- extra.pars$window.size;
    if(length(window.size)==0) window.size <- 1e6;
    res.find <- find.topN.variant(gene.name,window.size,singlevar.result,topN);
    raw.data.all <- list();
    known.variant.vec <- res.find$known.variant.collapsed.vec;
    p.value.known.variant.vec <- res.find$p.value.known.variant.vec;
    ##print(p.value.known.variant.vec);
    gene.name <- res.find$gene.name;
    candidate.variant.tabix.vec.tmp <- res.find$candidate.variant.tabix.vec;
    capture.output(raw.data.all <- rvmeta.readDataByRange(score.stat.file,NULL,candidate.variant.tabix.vec.tmp));
    res.list <- list();
    ## Get tabix for variants in gene
    candidate.variant.tabix.vec <- character(0);
    for(kk in 1:length(raw.data.all))
      {
        raw.data <- raw.data.all[[kk]];
        if(length(extra.pars$QC.par)>0) 
          raw.data <- QC(raw.data,extra.pars$QC.par,cov=0);        
        ix.var <- integer(0);
        if(ANNO=='gene') ix.var <- 1:length(raw.data$ustat[[1]]);
        if(ANNO!='gene')
          {
            if(length(raw.data$ustat[[1]])>0)
              {
                ix.var <- c(ix.var,grep(ANNO,raw.data$anno));
              }
          }
        if(length(ix.var)==0)
          {
            print(gene.name[kk],raw.data);
          }
        ##print(c(gene.name[kk],length(ix.var),raw.data$anno));
        ##print(known.variant.vec[kk]);
        candidate.variant.tabix.vec[kk] <- get.tabix.range(raw.data$pos[ix.var]);
      }
    res <- conditional.rareMETALS.basic(candidate.variant.tabix.vec,gene.name,score.stat.file,cov.file,known.variant.vec,test,maf.cutoff,no.boot,alternative,alpha,ix.gold,out.digits,callrate.cutoff,hwe.cutoff,gene.file,p.value.known.variant.vec);    
    return(res);
  }
r2cov.mat <- function(r2.mat,maf.vec)
  {
    var.vec <- sqrt(maf.vec*(1-maf.vec)*2);
    var.mat <- (var.vec%*%t(var.vec));
    cov.mat <- r2.mat*var.mat;
    return(cov.mat);
  }
