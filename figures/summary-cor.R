library(dREG);
library(parallel);
library(lmtest);

source("../scripts/hist.param.R");
source("../scripts/hist.svm.pred.R");
source("../scripts/hist.svm.com.R")


cor.K562.chr22 <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k122ac.bw , file.k562.H3k122ac.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==2) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k27ac.bw  , file.k562.H3k27ac.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==3) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k27me3.bw , file.k562.H3k27me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==4) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k36me3.bw , file.k562.H3k36me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==5) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me1.bw  , file.k562.H3k4me1.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==6) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me2.bw  , file.k562.H3k4me2.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==7) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me3.bw  , file.k562.H3k4me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==8) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k9ac.bw   , file.k562.H3k9ac.pred.bw   , chr.include="chr22",  sample_rate = 1  );
  if(i==9) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k9me3.bw  , file.k562.H3k9me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==10) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H4k20me1.bw , file.k562.H4k20me1.pred.bw , chr.include="chr22",  sample_rate = 1  );
  return(r); }, mc.cores=10) );

cor.GM.chr22  <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k122ac.bw , file.gm.H3k122ac.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==2) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27ac.bw  , file.gm.H3k27ac.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==3) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27me3.bw , file.gm.H3k27me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==4) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k36me3.bw , file.gm.H3k36me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==5) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me1.bw  , file.gm.H3k4me1.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==6) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me2.bw  , file.gm.H3k4me2.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==7) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me3.bw  , file.gm.H3k4me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==8) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9ac.bw   , file.gm.H3k9ac.pred.bw   , chr.include="chr22",  sample_rate = 1  );
  if(i==9) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9me3.bw  , file.gm.H3k9me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==10) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H4k20me1.bw , file.gm.H4k20me1.pred.bw , chr.include="chr22",  sample_rate = 1 );
  return(r); }, mc.cores=10) );
  

cor.CD4.chr22 <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k122ac.bw , file.cd4.H3k122ac.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==2) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k27ac.bw  , file.cd4.H3k27ac.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==3) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k27me3.bw , file.cd4.H3k27me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==4) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k36me3.bw , file.cd4.H3k36me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==5) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k4me1.bw  , file.cd4.H3k4me1.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==6) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k4me2.bw  , file.cd4.H3k4me2.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==7) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k4me3.bw  , file.cd4.H3k4me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==8) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k9ac.bw   , file.cd4.H3k9ac.pred.bw   , chr.include="chr22",  sample_rate = 1  );
  if(i==9) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H3k9me3.bw  , file.cd4.H3k9me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==10) r <- genomewide_cor(file.bw.CD4[1], file.bw.CD4[2], file.cd4.H4k20me1.bw , file.cd4.H4k20me1.pred.bw , chr.include="chr22",  sample_rate = 1 );
  return(r); }, mc.cores=10) );
  
cor.HELA.chr22 <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k122ac.bw , file.hela.H3k122ac.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==2) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k27ac.bw  , file.hela.H3k27ac.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==3) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k27me3.bw , file.hela.H3k27me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==4) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k36me3.bw , file.hela.H3k36me3.pred.bw , chr.include="chr22",  sample_rate = 1  ); 
  if(i==5) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k4me1.bw  , file.hela.H3k4me1.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==6) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k4me2.bw  , file.hela.H3k4me2.pred.bw  , chr.include="chr22",  sample_rate = 1  );
  if(i==7) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k4me3.bw  , file.hela.H3k4me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==8) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k9ac.bw   , file.hela.H3k9ac.pred.bw   , chr.include="chr22",  sample_rate = 1  );
  if(i==9) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H3k9me3.bw  , file.hela.H3k9me3.pred.bw  , chr.include="chr22",  sample_rate = 1  ); 
  if(i==10) r <- genomewide_cor(file.bw.HELA[1], file.bw.HELA[2], file.hela.H4k20me1.bw , file.hela.H4k20me1.pred.bw , chr.include="chr22",  sample_rate = 1 );     
  return(r); }, mc.cores=10) );
 
save(cor.K562.chr22, cor.GM.chr22, cor.CD4.chr22, cor.HELA.chr22, file="summary-cor.Rdata");   

cor.K562.all  <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k122ac.bw , file.k562.H3k122ac.pred.bw , chr.include=NULL,  sample_rate = 1 ); 
  if(i==2) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k27ac.bw  , file.k562.H3k27ac.pred.bw  , chr.include=NULL,  sample_rate = 1 ); 
  if(i==3) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k27me3.bw , file.k562.H3k27me3.pred.bw , chr.include=NULL,  sample_rate = 1 ); 
  if(i==4) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k36me3.bw , file.k562.H3k36me3.pred.bw , chr.include=NULL,  sample_rate = 1 ); 
  if(i==5) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me1.bw  , file.k562.H3k4me1.pred.bw  , chr.include=NULL,  sample_rate = 1 );
  if(i==6) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me2.bw  , file.k562.H3k4me2.pred.bw  , chr.include=NULL,  sample_rate = 1 );
  if(i==7) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me3.bw  , file.k562.H3k4me3.pred.bw  , chr.include=NULL,  sample_rate = 1 ); 
  if(i==8) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k9ac.bw   , file.k562.H3k9ac.pred.bw   , chr.include=NULL,  sample_rate = 1 );
  if(i==9) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H3k9me3.bw  , file.k562.H3k9me3.pred.bw  , chr.include=NULL,  sample_rate = 1 ); 
  if(i==10) r <- genomewide_cor(file.bw.G1[1], file.bw.G1[2], file.k562.H4k20me1.bw , file.k562.H4k20me1.pred.bw , chr.include=NULL,  sample_rate = 1); 
  return(r); }, mc.cores=10) );

cor.GM.all  <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k122ac.bw , file.gm.H3k122ac.pred.bw , chr.include=NULL,  sample_rate = 1 ); 
  if(i==2) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27ac.bw  , file.gm.H3k27ac.pred.bw  , chr.include=NULL,  sample_rate = 1 ); 
  if(i==3) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27me3.bw , file.gm.H3k27me3.pred.bw , chr.include=NULL,  sample_rate = 1 ); 
  if(i==4) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k36me3.bw , file.gm.H3k36me3.pred.bw , chr.include=NULL,  sample_rate = 1 ); 
  if(i==5) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me1.bw  , file.gm.H3k4me1.pred.bw  , chr.include=NULL,  sample_rate = 1 );
  if(i==6) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me2.bw  , file.gm.H3k4me2.pred.bw  , chr.include=NULL,  sample_rate = 1 );
  if(i==7) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me3.bw  , file.gm.H3k4me3.pred.bw  , chr.include=NULL,  sample_rate = 1 ); 
  if(i==8) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9ac.bw   , file.gm.H3k9ac.pred.bw   , chr.include=NULL,  sample_rate = 1 );
  if(i==9) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9me3.bw  , file.gm.H3k9me3.pred.bw  , chr.include=NULL,  sample_rate = 1 ); 
  if(i==10) r <- genomewide_cor(file.bw.GM[1], file.bw.GM[2], file.gm.H4k20me1.bw , file.gm.H4k20me1.pred.bw , chr.include=NULL,  sample_rate = 1);
  return(r); }, mc.cores=10) );
  

save(cor.K562.chr22, cor.GM.chr22, cor.CD4.chr22, cor.HELA.chr22, cor.GM.all, cor.K562.all, file="summary-cor.Rdata");   
   
genomewide_cor_eachchr <- function(file.bw.plus, file.bw.minus, file.bw.org, file.bw.imputed, cor.method, sample_rate)   
{
   r <- do.call("rbind", mclapply(1:22, function(i){
     r0 <- genomewide_cor( file.bw.plus, file.bw.minus, file.bw.org, file.bw.imputed, chr.include=paste0("chr", i),  sample_rate = 1 );
     return(r0);
     }, mc.cores=11));

   if( cor.method=="pearson") return (r[,1] ) else return (r[,2] );
}
   
cor.GM.each.chr  <- do.call("rbind", lapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k122ac.bw , file.gm.H3k122ac.pred.bw , "pearson", sample_rate = 1 ); 
  if(i==2) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27ac.bw  , file.gm.H3k27ac.pred.bw  , "pearson", sample_rate = 1 ); 
  if(i==3) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27me3.bw , file.gm.H3k27me3.pred.bw , "pearson", sample_rate = 1 ); 
  if(i==4) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k36me3.bw , file.gm.H3k36me3.pred.bw , "pearson", sample_rate = 1 ); 
  if(i==5) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me1.bw  , file.gm.H3k4me1.pred.bw  , "pearson", sample_rate = 1 );
  if(i==6) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me2.bw  , file.gm.H3k4me2.pred.bw  , "pearson", sample_rate = 1 );
  if(i==7) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me3.bw  , file.gm.H3k4me3.pred.bw  , "pearson", sample_rate = 1 ); 
  if(i==8) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9ac.bw   , file.gm.H3k9ac.pred.bw   , "pearson", sample_rate = 1 );
  if(i==9) r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9me3.bw  , file.gm.H3k9me3.pred.bw  , "pearson", sample_rate = 1 ); 
  if(i==10)r <- genomewide_cor_eachchr(file.bw.GM[1], file.bw.GM[2], file.gm.H4k20me1.bw , file.gm.H4k20me1.pred.bw , "pearson", sample_rate = 1);
  return(r); } ) );


save(cor.K562.chr22, cor.GM.chr22, cor.CD4.chr22, cor.HELA.chr22, cor.GM.all, cor.K562.all, cor.GM.each.chr , file="summary-cor.Rdata");   


genomewide_cor_eachchr <- function(file.bw.plus, file.bw.minus, file.bw.org, file.bw.imputed, cor.method, sample_rate)   
{
   r <- do.call("rbind", mclapply(1:22, function(i){
     r0 <- genomewide_cor( file.bw.plus, file.bw.minus, file.bw.org, file.bw.imputed, chr.include=paste0("chr", i),  sample_rate = 1 );
     return(r0);
     }, mc.cores=11));

   if( cor.method=="pearson") return (r[,1] ) else return (r[,2] );
}
   
genomewide_cor_peaks <- function(file.peak, file.bw.org, file.bw.imputed, cor.method="pearson", dist=1000)
{
   if( is.na(file.bw.org) || is.na(file.bw.imputed) )
       return(c( pearson=NA, spearman=NA ));
   
   if( !all(file.exists(c(file.bw.org, file.bw.imputed))))
       return(c( pearson=NA, spearman=NA ));

    exclude_chromosome = "_|chrM|chrY|chrX"; 
    tb <- read.table(file.peak);
    tb <- tb [ grep( exclude_chromosome , tb[,1], invert=TRUE), ]
    tb[,2] <- tb[,2]-100; 
    if (sum(tb[,2]<0)>0) tb[tb[,2]<0, 2] <- 0;
    tb[,3] <- tb[,3]+100; 

    file.tmp <- tempfile(fileext=".bed");
    write.bed(tb, file.tmp);

    tb1 <- read.table(pipe(paste("cat ", file.tmp, " | bedtools merge -i - -d 1000 | sort-bed -")));
    write.bed(tb1, file.tmp);

    tb0 <- read.table(pipe(paste("cat ", file.tmp, " | sort-bed - | bedtools complement -i - -g ../scripts/hg19.chromInfo " )));
    tb_bed <- rbind(tb0[,c(1:3)], tb1[,c(1:3)]);
    tb_bed <- tb_bed[ order(tb_bed[,1], tb_bed[,2]),]
    
    err <- try(load.bigWig(file.bw.imputed), silent=TRUE);
    if(class(err)=="try-error")
       pred <- get_bedgraph_read(tb_bed, file.bw.imputed)
    else
       pred <- get_histone_read(tb_bed, file.bw.imputed);
   
    orgl <- get_histone_read(tb_bed, file.bw.org)
   
    orgl <- orgl/(tb_bed[,3]-tb_bed[,2]);
    pred <- pred/(tb_bed[,3]-tb_bed[,2]);

library(TSdist);
    r.cor1 <- cor(pred, orgl, method ="pearson");
    r.cor2 <- cor(pred, orgl, method ="spearman");
    r.gt <- grangertest(orgl, pred, order=1)
    r.ccordist <- CCorDistance(orgl, pred );
    r.cordist  <- CorDistance(orgl, pred);
library(TSclust);
    r.DWT <- diss.DWT( rbind(orgl, pred) );

    cat("COR=", r.cor1, r.cor2, r.ccordist, r.cordist, r.gt$"Pr(>F)"[2], "\n" );
 
    return(c(pearson=r.cor1,  spearman=r.cor2, r.ccordist=r.ccordist, r.cordist=r.cordist,  r.DWT = r.DWT, pv=r.gt$"Pr(>F)"[2] ));    
}
   
cor.GM.peaks  <- do.call("rbind", lapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor_peaks( file.gm.H3k122ac.peak , file.gm.H3k122ac.bw , file.gm.H3k122ac.pred.bw ); 
  if(i==2) r <- genomewide_cor_peaks( file.gm.H3k27ac.peak  , file.gm.H3k27ac.bw  , file.gm.H3k27ac.pred.bw  ); 
  if(i==3) r <- genomewide_cor_peaks( file.gm.H3k27me3.peak , file.gm.H3k27me3.bw , file.gm.H3k27me3.pred.bw ); 
  if(i==4) r <- genomewide_cor_peaks( file.gm.H3k36me3.peak , file.gm.H3k36me3.bw , file.gm.H3k36me3.pred.bw ); 
  if(i==5) r <- genomewide_cor_peaks( file.gm.H3k4me1.peak  , file.gm.H3k4me1.bw  , file.gm.H3k4me1.pred.bw  );
  if(i==6) r <- genomewide_cor_peaks( file.gm.H3k4me2.peak  , file.gm.H3k4me2.bw  , file.gm.H3k4me2.pred.bw  );
  if(i==7) r <- genomewide_cor_peaks( file.gm.H3k4me3.peak  , file.gm.H3k4me3.bw  , file.gm.H3k4me3.pred.bw  ); 
  if(i==8) r <- genomewide_cor_peaks( file.gm.H3k9ac.peak   , file.gm.H3k9ac.bw   , file.gm.H3k9ac.pred.bw   );
  if(i==9) r <- genomewide_cor_peaks( file.gm.H3k9me3.peak  , file.gm.H3k9me3.bw  , file.gm.H3k9me3.pred.bw  ); 
  if(i==10)r <- genomewide_cor_peaks( file.gm.H4k20me1.peak , file.gm.H4k20me1.bw , file.gm.H4k20me1.pred.bw );
  return(r); } ) );


cor.K562.peaks  <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor_peaks( file.k562.H3k122ac.peak , file.k562.H3k122ac.bw , file.k562.H3k122ac.pred.bw ); 
  if(i==2) r <- genomewide_cor_peaks( file.k562.H3k27ac.peak ,  file.k562.H3k27ac.bw  , file.k562.H3k27ac.pred.bw  ); 
  if(i==3) r <- genomewide_cor_peaks( file.k562.H3k27me3.peak , file.k562.H3k27me3.bw , file.k562.H3k27me3.pred.bw ); 
  if(i==4) r <- genomewide_cor_peaks( file.k562.H3k36me3.peak , file.k562.H3k36me3.bw , file.k562.H3k36me3.pred.bw ); 
  if(i==5) r <- genomewide_cor_peaks( file.k562.H3k4me1.peak ,  file.k562.H3k4me1.bw  , file.k562.H3k4me1.pred.bw  );
  if(i==6) r <- genomewide_cor_peaks( file.k562.H3k4me2.peak ,  file.k562.H3k4me2.bw  , file.k562.H3k4me2.pred.bw  );
  if(i==7) r <- genomewide_cor_peaks( file.k562.H3k4me3.peak ,  file.k562.H3k4me3.bw  , file.k562.H3k4me3.pred.bw  ); 
  if(i==8) r <- genomewide_cor_peaks( file.k562.H3k9ac.peak ,   file.k562.H3k9ac.bw   , file.k562.H3k9ac.pred.bw   );
  if(i==9) r <- genomewide_cor_peaks( file.k562.H3k9me3.peak ,  file.k562.H3k9me3.bw  , file.k562.H3k9me3.pred.bw  ); 
  if(i==10) r <- genomewide_cor_peaks(file.k562.H4k20me1.peak , file.k562.H4k20me1.bw , file.k562.H4k20me1.pred.bw ); 
  return(r); }, mc.cores=2) );


genomewide_cor_window <- function( file.plus, file.minus, file.bw.org, file.bw.imputed, windows=10000)
{
   if( is.na(file.bw.org) || is.na(file.bw.imputed) )
       return(c( pearson=NA, spearman=NA ));
   
   if( !all(file.exists(c(file.bw.org, file.bw.imputed))))
       return(c( pearson=NA, spearman=NA ));

    chrom.info.table <- get.chromosome.info( file.plus, file.minus );

    dnase_bed <- as.data.frame(rbindlist( lapply(1:NROW(chrom.info.table), function(i){ 
       df <- data.frame(chrom.info.table[i,1], seq(1, (chrom.info.table[i,2]-1), windows), seq(1, (chrom.info.table[i,2]-1), windows)+(windows-1)) 
       df <- df[-NROW(df),];
       return(df);
    }) ));
  
    colnames(dnase_bed)<-NULL

    exclude_chromosome = "_|chrM|chrY|chrX"; 
    dnase_bed <- dnase_bed [ grep( exclude_chromosome , dnase_bed[,1], invert=TRUE), ]
    
    err <- try(load.bigWig(file.bw.imputed), silent=TRUE);
    if(class(err)=="try-error")
       pred <- get_bedgraph_read(dnase_bed, file.bw.imputed)
    else
       pred <- get_histone_read(dnase_bed, file.bw.imputed);
   
    orgl <- get_histone_read(dnase_bed, file.bw.org)
   
    orgl <- orgl/(dnase_bed[,3]-dnase_bed[,2]);
    pred <- pred/(dnase_bed[,3]-dnase_bed[,2]);

library(TSdist);
    r.cor1 <- cor(pred, orgl, method ="pearson");
    r.cor2 <- cor(pred, orgl, method ="spearman");
    r.gt <- grangertest(orgl, pred, order=1)
    r.ccordist <- CCorDistance(orgl, pred );
    r.cordist  <- CorDistance(orgl, pred);
library(TSclust);
    r.DWT <- diss.DWT( rbind(orgl, pred) );

    cat("COR=", r.cor1, r.cor2, r.ccordist, r.cordist, r.DWT, r.gt$"Pr(>F)"[2], "\n" );
 
    return(c(pearson=r.cor1,  spearman=r.cor2, r.ccordist=r.ccordist, r.cordist=r.cordist,  r.DWT = r.DWT, pv=r.gt$"Pr(>F)"[2] ));    
}  


cor.K562.wins  <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k122ac.bw , file.k562.H3k122ac.pred.bw ); 
  if(i==2) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k27ac.bw  , file.k562.H3k27ac.pred.bw  ); 
  if(i==3) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k27me3.bw , file.k562.H3k27me3.pred.bw ); 
  if(i==4) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k36me3.bw , file.k562.H3k36me3.pred.bw ); 
  if(i==5) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me1.bw  , file.k562.H3k4me1.pred.bw  );
  if(i==6) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me2.bw  , file.k562.H3k4me2.pred.bw  );
  if(i==7) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k4me3.bw  , file.k562.H3k4me3.pred.bw  ); 
  if(i==8) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k9ac.bw   , file.k562.H3k9ac.pred.bw   );
  if(i==9) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H3k9me3.bw  , file.k562.H3k9me3.pred.bw  ); 
  if(i==10) r <- genomewide_cor_window(file.bw.G1[1], file.bw.G1[2], file.k562.H4k20me1.bw , file.k562.H4k20me1.pred.bw); 
  return(r); }, mc.cores=5) );

cor.GM.wins  <- do.call("rbind", mclapply(1:10, function(i){
  r <- NA;
  if(i==1) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k122ac.bw , file.gm.H3k122ac.pred.bw ); 
  if(i==2) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27ac.bw  , file.gm.H3k27ac.pred.bw  ); 
  if(i==3) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k27me3.bw , file.gm.H3k27me3.pred.bw ); 
  if(i==4) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k36me3.bw , file.gm.H3k36me3.pred.bw ); 
  if(i==5) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me1.bw  , file.gm.H3k4me1.pred.bw  );
  if(i==6) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me2.bw  , file.gm.H3k4me2.pred.bw  );
  if(i==7) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k4me3.bw  , file.gm.H3k4me3.pred.bw  ); 
  if(i==8) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9ac.bw   , file.gm.H3k9ac.pred.bw   );
  if(i==9) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H3k9me3.bw  , file.gm.H3k9me3.pred.bw  ); 
  if(i==10) r <- genomewide_cor_window(file.bw.GM[1], file.bw.GM[2], file.gm.H4k20me1.bw , file.gm.H4k20me1.pred.bw );
  return(r); }, mc.cores=5) );


save(cor.K562.chr22, cor.GM.chr22, cor.CD4.chr22, cor.HELA.chr22, cor.GM.all, cor.K562.all, cor.GM.each.chr, cor.K562.wins, cor.GM.wins, cor.K562.peaks, cor.GM.peaks, file="summary-cor.Rdata");   

