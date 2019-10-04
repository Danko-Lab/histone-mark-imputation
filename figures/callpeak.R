library(dREG);
library(parallel);

source("/local/workdir/zw355/proj/prj15-histone/scripts/hist.param.R");
source("/local/workdir/zw355/proj/prj15-histone/scripts/hist.svm.pred.R");
source("/local/workdir/zw355/proj/prj15-histone/scripts/hist.svm.com.R")


macs_bdgbroadcall<- function(file.pred.bdg, file.org.bdg, par.C, par.c, par.l)
{
    file.pred.tmp <- tempfile(fileext=".bed");
    cmd <- paste("source /programs/bin/util/setup_macs2.sh; export PYTHONPATH=/programs/MACS2-2.1.1/lib64/python2.7/site-packages:$PYTHONPATH ; /programs/MACS2-2.1.1/bin/macs2  bdgbroadcall -i ", file.pred.bdg, "-C", par.C, "-c", par.c, "-l", par.l, " -o ", file.pred.tmp );
print(cmd);
   system(cmd);

    file.org.tmp <- tempfile(fileext=".bed");
    cmd <- paste("source /programs/bin/util/setup_macs2.sh; export PYTHONPATH=/programs/MACS2-2.1.1/lib64/python2.7/site-packages:$PYTHONPATH ; /programs/MACS2-2.1.1/bin/macs2  bdgbroadcall -i ", file.org.bdg, "-C", par.C, "-c", par.c, "-l", par.l, " -o ", file.org.tmp );
print(cmd);
    system(cmd);

    if (file.exists(file.pred.tmp) && file.exists(file.org.tmp))
    {
      tb.pred <- tb.org <- c();
      try( tb.pred <- read.table(file.pred.tmp, skip=1), silent=T );
      try( tb.org <- read.table(file.org.tmp, skip=1), silent=T );
      
      tb0 <- tb1 <- tb2 <- c();
      if(NROW(tb.pred)>0 && NROW(tb.org)>0)
      {
         tb0 <- read.table(pipe(paste("bedtools intersect -a " , file.pred.tmp, "-b", file.org.tmp, "-wa")))
         tb1 <- read.table(pipe(paste("bedtools intersect -a " , file.pred.tmp, "-b", file.org.tmp, "-v")))
         tb2 <- read.table(pipe(paste("bedtools intersect -b " , file.pred.tmp, "-a", file.org.tmp, "-v")))
	  }
	  
	  unlink(file.org.tmp);
	  unlink(file.pred.tmp);
	  ret <- c( NROW(tb.org),  NROW(tb2), NROW(tb0), NROW(tb1), NROW(tb.pred) );

print(ret);
      return(ret);
      
    }
    else
      return(rep(NA,5));
}

macs_bdgpeakcall<- function(file.bedgraph, file.org.peak, par.c, par.l)
{
    file.tmp <- tempfile(fileext=".bed");
    cmd <- paste("source /programs/bin/util/setup_macs2.sh; export PYTHONPATH=/programs/MACS2-2.1.1/lib64/python2.7/site-packages:$PYTHONPATH ; /programs/MACS2-2.1.1/bin/macs2  bdgpeakcall -i ", file.bedgraph, "-c", par.c, "-l", par.l, " -o ", file.tmp );
    system(cmd);

    if (file.exists(file.tmp))
    {
      tb.call <- read.table(file.tmp, skip=1);
      tb.true <- read.table(file.org.peak);
      tb0 <- read.table(pipe(paste("bedtools intersect -a " , file.tmp, "-b", file.org.peak, "-wa")))
      tb1 <- read.table(pipe(paste("bedtools intersect -a " , file.org.peak, "-b", file.tmp, "-wa")))
      
      tb.call <- unique(tb.call[,c(1:3)] )
      tb.true <- unique(tb.true[,c(1:3)] )
      tb0 <- unique(tb0[,c(1:3)] )
      tb1 <- unique(tb1[,c(1:3)] )
      
      return(c( NROW(tb.call) - NROW(tb0),NROW(tb0), NROW(tb1), NROW(tb.true) - NROW(tb1)));
    }
    else
      return(rep(NA,4));
}


broad_peakcall <- function(file.imputed, file.org.peak, cores=20)
{
  library(tools);

  if(file_ext(file.imputed) != "gz" )
  {
     file.input <- tempfile(fileext=".bed");
     system(paste("bigWigToBedGraph", file.imputed,  file.input))
     tb <- read.table(file.input);
     unlink(file.input)
  }   
  else
     tb <- read.table(file.imputed);
   
  tb[,3] <- tb[,3] + 8;
  tb <- tb[-NROW(tb),]
  file.bedgraph <- tempfile(fileext=".bed");
  write.bed(tb, file.bedgraph);   
  
  gc();
   
  par <- expand.grid(par.C=2:10, par.c=3:11)
  par <- par [ par$par.C<par$par.c, ]

  r <- do.call("rbind", mclapply(1:NROW(par), function(i) {
       r0 <- macs_bdgbroadcall (file.bedgraph, file.org.peak, par$par.C[i], par$par.c[i])
       r0 <- c(par.C=par$par.C[i], par.c=par$par.c[i], r0);
       gc();
       return(r0);
   }, mc.cores=cores ) )    
  
  gc();
  unlink(file.bedgraph);
  return(r); 
}       

narrow_peakcall <- function(file.imputed, file.org.peak, cores=20)
{
  library(tools);

  if(file_ext(file.imputed) != "gz" )
  {
     file.input <- tempfile(fileext=".bed");
     system(paste("bigWigToBedGraph", file.imputed,  file.input))
     tb <- read.table(file.input);
     unlink(file.input)
  }   
  else
     tb <- read.table(file.imputed);
   
  tb[,3] <- tb[,3] + 8;
  tb <- tb[-NROW(tb),]
  file.bedgraph <- tempfile(fileext=".bed");
  write.bed(tb, file.bedgraph);   
  
  r <- do.call("rbind", mclapply(3:11, function(i) {
       r0 <- macs_bdgpeakcall (file.bedgraph, file.org.peak, i)
       r0 <- c(par.C=0, par.c=i, r0);
       return(r0);
   }, mc.cores=cores ) )    
  
  unlink(file.bedgraph);
  return(r); 
}       


broad_peakcall<-function(file.imputed.bw, file.org.bw, par.c=5, par.C=2, par.l=500)
{
  library(tools);

  file.pred.bdg <- tempfile(fileext=".bed");
  system(paste("bigWigToBedGraph", file.imputed.bw,  file.pred.bdg))

  file.org.bdg <- tempfile(fileext=".bed");
  system(paste("bigWigToBedGraph", file.org.bw,  file.org.bdg))
   
  r0 <- macs_bdgbroadcall ( file.pred.bdg, file.org.bdg, par.C, par.c, par.l)

  unlink(file.pred.bdg);
  unlink(file.org.bdg);

  return(r0); 
}       

r.GM <- list(
  H3k27me3 <- broad_peakcall( file.gm.H3k27me3.pred.bw , file.gm.H3k27me3.bw),
  H3k27ac  <- broad_peakcall( file.gm.H3k27ac.pred.bw  , file.gm.H3k27ac.bw ), 
  H3k36me3 <- broad_peakcall( file.gm.H3k36me3.pred.bw , file.gm.H3k36me3.bw ), 
  H3k4me1  <- broad_peakcall( file.gm.H3k4me1.pred.bw  , file.gm.H3k4me1.bw ),
  H3k4me2  <- broad_peakcall( file.gm.H3k4me2.pred.bw  , file.gm.H3k4me2.bw ),
  H3k4me3  <- broad_peakcall( file.gm.H3k4me3.pred.bw  , file.gm.H3k4me3.bw ), 
  H3k9ac   <- broad_peakcall( file.gm.H3k9ac.pred.bw   , file.gm.H3k9ac.bw ),
  H3k9me3  <- broad_peakcall( file.gm.H3k9me3.pred.bw  , file.gm.H3k9me3.bw ), 
  H4k20me1 <- broad_peakcall( file.gm.H4k20me1.pred.bw , file.gm.H4k20me1.bw ))

r.K562 <- list(
  H3k27me3 <- broad_peakcall( file.k562.H3k27me3.pred.bw , file.k562.H3k27me3.bw),
  H3k122ac <- broad_peakcall( file.k562.H3k122ac.pred.bw , file.k562.H3k122ac.bw ),
  H3k27ac  <- broad_peakcall( file.k562.H3k27ac.pred.bw  , file.k562.H3k27ac.bw ), 
  H3k36me3 <- broad_peakcall( file.k562.H3k36me3.pred.bw , file.k562.H3k36me3.bw ), 
  H3k4me1  <- broad_peakcall( file.k562.H3k4me1.pred.bw  , file.k562.H3k4me1.bw ),
  H3k4me2  <- broad_peakcall( file.k562.H3k4me2.pred.bw  , file.k562.H3k4me2.bw ),
  H3k4me3  <- broad_peakcall( file.k562.H3k4me3.pred.bw  , file.k562.H3k4me3.bw ), 
  H3k9ac   <- broad_peakcall( file.k562.H3k9ac.pred.bw   , file.k562.H3k9ac.bw ),
  H3k9me3  <- broad_peakcall( file.k562.H3k9me3.pred.bw  , file.k562.H3k9me3.bw ), 
  H4k20me1 <- broad_peakcall( file.k562.H4k20me1.pred.bw , file.k562.H4k20me1.bw ))

ret.K562 <- as.data.frame(do.call("rbind", r.K562));
rownames(ret.K562) <- c(
  "H3k27me3",
  "H3k122ac",
  "H3k27ac" ,
  "H3k36me3",
  "H3k4me1" ,
  "H3k4me2" ,
  "H3k4me3" ,
  "H3k9ac"  ,
  "H3k9me3" ,
  "H4k20me1");
  
ret.GM <- as.data.frame(do.call("rbind", r.GM));
rownames(ret.GM) <- c(
  "H3k27me3",
  "H3k27ac" ,
  "H3k36me3",
  "H3k4me1" ,
  "H3k4me2" , 
  "H3k4me3" , 
  "H3k9ac"  , 
  "H3k9me3" , 
  "H4k20me1" );


colnames(ret.GM) <- c("org", "org_v", "ovelap", "pred_v", "pred")
colnames(ret.K562) <- c("org", "org_v", "ovelap", "pred_v", "pred")

ret.GM$org_ratio   <- ret.GM$org_v/ret.GM$org;
ret.GM$pred_ratio  <- ret.GM$pred_v/ret.GM$pred;
ret.GM$org_over    <- ret.GM$org - ret.GM$org_v;
ret.GM$pred_over   <- ret.GM$pred - ret.GM$pred_v;
ret.GM <- ret.GM[, c("org", "org_ratio", "org_v", "org_over", "pred_over", "pred_v", "pred_ratio", "pred")]

ret.K562$org_ratio  <- ret.K562$org_v/ret.K562$org;
ret.K562$pred_ratio <- ret.K562$pred_v/ret.K562$pred;
ret.K562$org_over   <- ret.K562$org - ret.K562$org_v;
ret.K562$pred_over  <- ret.K562$pred - ret.K562$pred_v;
ret.K562 <- ret.K562[, c("org", "org_ratio", "org_v", "org_over", "pred_over", "pred_v", "pred_ratio", "pred")]


save( r.K562, r.GM, ret.GM, ret.K562, file="callpeak.rdata");