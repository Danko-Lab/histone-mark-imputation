library(bigWig);
library(MASS);
library(parallel)

source("../scripts/hist.param.R");

if(1)
{
	file_H3k27ac_bw  = file.k562.H3k27ac.pred.bw
	file_H3k27me3_bw = file.k562.H3k27me3.pred.bw
	file_H3k36me3_bw = file.k562.H3k36me3.pred.bw
	file_H3k9me3_bw  = file.k562.H3k9me3.pred.bw
	file_H3k4me1_bw  = file.k562.H3k4me1.pred.bw
	file_H3k4me3_bw  = file.k562.H3k4me3.pred.bw

	file_signal_prefix_bed <- "~/temp/K562.chromHMM.pred.adjust"
	file_signal_pdf <- "K562.hist.pred.adjust.dist.pdf"
	file_signal_make <- "spf_k562/k562"
	file_signal_rdata <- "K562.hist.pred.adjust.rdata"
	path_chromHMM <- "spf_k562"
	path_chromHMM2 <- "MS_K562_pred_adjust_seg"
	file.chrom.ret.bed.gz <- "~/temp/MS_K562_pred_adjust_seg.bed.gz"
}


write.bed<-function ( df.bed, file.bed, compress=FALSE )
{
	options("scipen"=100, "digits"=4);
	temp <- tempfile(fileext=".bed");
	write.table( df.bed, file=temp, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");
	if(compress)
 	   system(paste0("sort-bed ", temp,  " | bgzip > ",  file.bed ))
	else
 	   system(paste0("sort-bed ", temp,  " > ",  file.bed ));
 	   
	invisible(unlink(temp));
}


write.signal.bed<-function(df.bed, file.zip, compress=TRUE)
{
     df.bed <- df.bed[which(df.bed[,4]!=0),]
     write.bed(df.bed, file.zip, compress=compress);
}

read.bigwig.fast <- function(file.hist, tb.bed, op="sum", inc=1 )
{
  interval <- unique(c( seq( 1, NROW(tb.bed)+1, by = 1000*100 ), NROW(tb.bed)+1))

  ret <- do.call("c", mclapply(1:(length(interval)-1), function(x)
  {
    bw <- load.bigWig(file.hist);
    batch_indx<- c( interval[x]:(interval[x+1]-1) ) 
    dat <- bed.region.bpQuery.bigWig( bw, tb.bed[batch_indx, ], op=op);
    unload.bigWig(bw);
    return( c(dat)*inc );
  }, mc.cores=24));
  
  return(ret);
}

tb.chrom <- read.table("/fs/cbsudanko/storage/data/hg19/hg19.chromInfo");
tb.chrom <- tb.chrom[grep("_|chrM|chrY|chrX", tb.chrom[,1], invert=TRUE),]

tb.bed <-do.call("rbind", lapply(1:NROW(tb.chrom), function(i){
   return(data.frame(chr=tb.chrom[i,1], start=seq(1, tb.chrom[i,2]-200, 200), stop=seq(1, tb.chrom[i,2]-200, 200)+200 ));
}));

tb.h3k27ac.signal <- tb.h3k27me3.signal <- tb.h3k36me3.signal <- tb.h3k9me3.signal <- tb.h3k4me1.signal <- tb.h3k4me3.signal <- rep( 0, NROW(tb.bed));

if(1)
{
	tb <- read.table("/workdir/cgd24/histoneImputation/k562.H3k27ac.imp.bed.gz");
	idx <- match(paste(tb$V1,tb$V3, sep=":"), paste(tb.bed[,1], tb.bed[,3]-1, sep=":") )
	tb.h3k27ac.signal[ idx ]<- 1

	tb <- read.table("/workdir/cgd24/histoneImputation/k562.H3k27me3.imp.bed.gz");
	idx <- match(paste(tb$V1,tb$V3, sep=":"), paste(tb.bed[,1], tb.bed[,3]-1, sep=":") )
	tb.h3k27me3.signal[ idx ]<- 1

	tb <- read.table("/workdir/cgd24/histoneImputation/k562.H3k36me3.imp.bed.gz");
	idx <- match(paste(tb$V1,tb$V3, sep=":"), paste(tb.bed[,1], tb.bed[,3]-1, sep=":") )
	tb.h3k36me3.signal[ idx ]<- 1

	tb <- read.table("/workdir/cgd24/histoneImputation/k562.H3k9me3.imp.bed.gz");
	idx <- match(paste(tb$V1,tb$V3, sep=":"), paste(tb.bed[,1], tb.bed[,3]-1, sep=":") )
	tb.h3k9me3.signal[ idx ]<- 1

	tb <- read.table("/workdir/cgd24/histoneImputation/k562.H3k4me3.imp.bed.gz");
	idx <- match(paste(tb$V1,tb$V3, sep=":"), paste(tb.bed[,1], tb.bed[,3]-1, sep=":") )
	tb.h3k4me3.signal[ idx ]<- 1

	tb <- read.table("/workdir/cgd24/histoneImputation/k562.H3k4me1.imp.bed.gz");
	idx <- match(paste(tb$V1,tb$V3, sep=":"), paste(tb.bed[,1], tb.bed[,3]-1, sep=":") )
	tb.h3k4me1.signal[ idx ]<- 1
}


tb.signal <- data.frame( H3K27ac=tb.h3k27ac.signal, H3K27me3=tb.h3k27me3.signal, H3K36me3=tb.h3k36me3.signal, H3K4me1=tb.h3k4me1.signal, H3K4me3=tb.h3k4me3.signal, H3K9me3=tb.h3k9me3.signal );

L <- lapply(as.character(tb.chrom[,1]), function(chr)
{
   file.out <- paste(file_signal_make, chr, "binary.txt", sep="_");
   idx <- sort(which(tb.bed[,1]==chr))
cat(chr, "=", min(idx), max(idx),NROW(idx),"\n");
   cat("k562", chr, "\n", file=file.out, sep="\t", append =FALSE)
   write.table(tb.signal[idx,],  file=file.out, row.names=FALSE, col.names=TRUE, append=TRUE, quote=FALSE, sep="\t")
});


system(paste("java -mx6400M -jar src/ChromHMM.jar MakeSegmentation model_18_core_K27ac.txt", path_chromHMM, path_chromHMM2) )

tb <- read.table(paste(path_chromHMM2, "k562_18_segments.bed", sep="/"))
tb$V5 <- substring(tb$V4,2);
for(i in 1:18) cat("state=i", i, sum(as.numeric(tb$V3-tb$V2)[tb$V5==i])/sum(as.numeric(tb$V3-tb$V2)), "\n" )
write.bed(tb[,c(1,2,3,5)], file.chrom.ret.bed.gz, compress=TRUE)


write.signal.bed(data.frame(tb.bed, tb.h3k27ac.signal),  paste(file_signal_prefix_bed, "H3K27ac.bed.gz", sep="."))
write.signal.bed(data.frame(tb.bed, tb.h3k27me3.signal), paste(file_signal_prefix_bed, "H3K27me3.bed.gz", sep="."))
write.signal.bed(data.frame(tb.bed, tb.h3k36me3.signal), paste(file_signal_prefix_bed, "H3K36me3.bed.gz", sep="."))
write.signal.bed(data.frame(tb.bed, tb.h3k9me3.signal),  paste(file_signal_prefix_bed, "H3K9me3.bed.gz", sep="."))
write.signal.bed(data.frame(tb.bed, tb.h3k4me1.signal),  paste(file_signal_prefix_bed, "H3K4me1.bed.gz", sep="."))
write.signal.bed(data.frame(tb.bed, tb.h3k4me3.signal),  paste(file_signal_prefix_bed, "H3K4me3.bed.gz", sep="."))


