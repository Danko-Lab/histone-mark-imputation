library(bigWig);
library(MASS);
library(parallel)


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


write.signal.bed<-function(df.bed, file.zip)
{
     df.bed <- df.bed[which(df.bed[,4]!=0),]
     write.bed(df.bed, file.zip, compress=TRUE);
}

read.bigwig.peak.bin <- function(file.hist.peak, file.tb.bed )
{
  tb <- read.table(pipe(paste("zcat ", file.hist.peak, " | bedtools coverage -a ", file.tb.bed, " -b -")))
  return( (as.numeric(tb$V4>0)));
}


file.k562.H3k27ac.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz";
file.k562.H3k27me3.peak <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k27me3StdAln.bed.gz";
file.k562.H3k36me3.peak <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k36me3StdAln.bed.gz";
file.k562.H3k4me1.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k4me1StdAln.bed.gz";
file.k562.H3k4me3.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k4me3StdAln.bed.gz";
file.k562.H3k9me3.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k9me3StdAln.bed.gz";

tb.chrom <- read.table("/fs/cbsudanko/storage/data/hg19/hg19.chromInfo");
tb.chrom <- tb.chrom[grep("_|chrM|chrY|chrX", tb.chrom[,1], invert=TRUE),]

tb.bed <-do.call("rbind", lapply(1:NROW(tb.chrom), function(i){
   return(data.frame(chr=tb.chrom[i,1], start=seq(1, tb.chrom[i,2]-200, 200), stop=seq(1, tb.chrom[i,2]-200, 200)+200-1 ));
}));

write.bed(tb.bed, file="temp.tb.bed");
tb.bed <- read.table("temp.tb.bed");

tb.h3k27ac.peak <- tb.h3k27me3.peak <- tb.h3k36me3.peak <- tb.h3k9me3.peak <- tb.h3k4me1.peak <- tb.h3k4me3.peak <- c();

tb.h3k27ac.peak  <- read.bigwig.peak.bin( file.k562.H3k27ac.peak,  "temp.tb.bed");
tb.h3k27me3.peak <- read.bigwig.peak.bin( file.k562.H3k27me3.peak, "temp.tb.bed");
tb.h3k36me3.peak <- read.bigwig.peak.bin( file.k562.H3k36me3.peak, "temp.tb.bed");
tb.h3k4me1.peak  <- read.bigwig.peak.bin( file.k562.H3k4me1.peak,  "temp.tb.bed");
tb.h3k4me3.peak  <- read.bigwig.peak.bin( file.k562.H3k4me3.peak,  "temp.tb.bed");
tb.h3k9me3.peak  <- read.bigwig.peak.bin( file.k562.H3k9me3.peak,  "temp.tb.bed");

tb.signal.peak <- data.frame( H3K4me3=tb.h3k4me3.peak, H3K27ac=tb.h3k27ac.peak, H3K4me1=tb.h3k4me1.peak,  H3K36me3=tb.h3k36me3.peak, H3K9me3=tb.h3k9me3.peak, H3K27me3=tb.h3k27me3.peak );
save(tb.signal.peak, file="k562.hist.peak.rdata");

tb.signal.peak <- tb.signal.peak[, c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3", "H3K9me3", "H3K27me3")]


L <- lapply(as.character(tb.chrom[,1]), function(chr)
{
   file.out <- paste0("Rmake0/k562_", chr, "_binary.txt");
   idx <- sort(which(tb.bed[,1]==chr))
   cat("k562", chr, "\n", file=file.out, sep="\t", append =FALSE)
cat(chr, "=", min(idx), max(idx),NROW(idx),"\n");   
   write.table(tb.signal.peak[idx,],  file=file.out, row.names=FALSE, col.names=TRUE, append=TRUE, quote=FALSE, sep="\t")
});

write.signal.bed(data.frame(tb.bed, tb.signal.peak$H3K27ac),  "~/temp/h3k27ac.chromHMM.signal.bed.gz")
write.signal.bed(data.frame(tb.bed, tb.signal.peak$H3K27me3),  "~/temp/h3k27me3.chromHMM.signal.bed.gz")
write.signal.bed(data.frame(tb.bed, tb.signal.peak$H3K36me3),  "~/temp/h3k36me3.chromHMM.signal.bed.gz")
write.signal.bed(data.frame(tb.bed, tb.signal.peak$H3K9me3),  "~/temp/h3k9me3.chromHMM.signal.bed.gz")
write.signal.bed(data.frame(tb.bed, tb.signal.peak$H3K4me1),  "~/temp/h3k4me1.chromHMM.signal.bed.gz")
write.signal.bed(data.frame(tb.bed, tb.signal.peak$H3K4me3),  "~/temp/h3k4me3.chromHMM.signal.bed.gz")

unlink("temp.tb.bed");

