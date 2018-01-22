conditional.rareMETALS.gene.byKnownVariant <- function(ANNO,score.stat.file,cov.file,known.variant.vec,singlevar.result,test='GRANVIL',maf.cutoff,no.boot=10000,alternative=c('two.sided','greater','less'),alpha=0.05,ix.gold=1,out.digits=4,window.size=1e6,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz")
  {
    extra.pars <- list(alpha=alpha,ix.gold=ix.gold,out.digits=out.digits,window.size=window.size,QC.par=list(callrate.cutoff=callrate.cutoff,hwe.cutoff=hwe.cutoff));
    window.size <- extra.pars$window.size;
    if(length(window.size)==0) window.size <- 1e6;
    res.find <- find.gene(known.variant.list,window.size,singlevar.result);
    candidate.variant.tabix.vec.tmp <- res.find$candidate.variant.tabix.vec;
    gene.name <- res.find$gene.vec;
    p.value.known.variant.vec <- res.find$p.value.known.variant.vec;
    all.variant.tabix.vec <- res.find$all.variant.tabix.vec;
    raw.data.all <- list();
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file,NULL,candidate.variant.tabix.vec.tmp));

    res.list <- list();
    candidate.variant.tabix.vec <- character(0);
    for(kk in 1:length(raw.data.all))
      {
        raw.data <- raw.data.all[[kk]];
        if(length(extra.pars$QC.par)>0) 
          raw.data <- QC(raw.data,extra.pars$QC.par,cov=1);        
        ix.var <- integer(0);
        if(ANNO=='gene') ix.var <- 1:length(raw.data$ustat[[1]]);
        if(ANNO!='gene')
          {
            if(length(raw.data$ustat[[1]])>0)
              {
                ix.var <- c(ix.var,grep(ANNO,raw.data$anno));
              }
          }
        candidate.variant.tabix.vec[kk] <- get.tabix.range(raw.data$pos[ix.var])
      }
    res <- conditional.rareMETALS.basic(candidate.variant.tabix.vec,gene.name,score.stat.file,cov.file,known.variant.vec,test,maf.cutoff,no.boot,alternative,alpha,ix.gold,out.digits,callrate.cutoff,hwe.cutoff,gene.file,p.value.known.variant.vec);
    return(res);    
  }

