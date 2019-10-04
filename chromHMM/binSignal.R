library(bigWig);
library(MASS);
library(parallel);

if(0)
{
	file_H3k27ac_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig"
	file_H3k27me3_bw = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig"
	file_H3k36me3_bw = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig"
	file_H3k9me3_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig"
	file_H3k4me1_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig"
	file_H3k4me3_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig"

	file_signal_prefix_bed <- "~/temp/K562.chromHMM.org.signal"
	file_signal_pdf <- "K562.hist.org.dist.pdf"
	file_signal_make <- "so_k562/k562"
	file_signal_rdata <- "K562.hist.org.rdata"
	path_chromHMM <- "so_k562"
	path_chromHMM2 <- "MS_K562_org_seg"
	file.chrom.ret.bed.gz <- "~/temp/MS_K562_org_seg.bed.gz"
}

if(0)
{
	file_H3k27ac_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig"
	file_H3k27me3_bw = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27me3StdSig.bigWig"
	file_H3k36me3_bw = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k36me3StdSig.bigWig"
	file_H3k9me3_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k9me3StdSig.bigWig"
	file_H3k4me1_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me1StdSig.bigWig"
	file_H3k4me3_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me3StdSig.bigWig"

	file_signal_prefix_bed <- "~/temp/GM12878.chromHMM.org.signal"
	file_signal_pdf <- "GM12878.hist.org.dist.pdf"
	file_signal_make <- "so_gm12878/gm12878"
	file_signal_rdata <- "GM12878.hist.org.rdata"
	path_chromHMM <- "so_gm12878"
	path_chromHMM2 <- "MS_GM12878_org_seg"
	file.chrom.ret.bed.gz <- "~/temp/MS_GM12878_org_seg.bed.gz"
}

if(1)
{
	file_H3k27ac_bw  = "/workdir/zw355/proj/prj15-histone/UwHistone/GSM2877103_ChIP-seq_K562_H3K27ac_rep1.bw"
	file_H3k27me3_bw = "/workdir/zw355/proj/prj15-histone/UwHistone/wgEncodeUwHistoneK562H3k27me3StdRaw.merge.bigWig"
	file_H3k36me3_bw = "/workdir/zw355/proj/prj15-histone/UwHistone/wgEncodeUwHistoneK562H3k36me3StdRaw.merge.bigWig"
	file_H3k9me3_bw  = "/workdir/zw355/proj/prj15-histone/UwHistone/GSM607494_K562_H3K9me3.bigWIG"
	#file_H3k4me1_bw  = "/workdir/zw355/proj/prj15-histone/UwHistone/GSM788085_wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig"
	file_H3k4me1_bw  = "/workdir/zw355/proj/prj15-histone/UwHistone/GSM3315712_k562_h3k4me1_nownt.bw"
	file_H3k4me3_bw  = "/workdir/zw355/proj/prj15-histone/UwHistone/wgEncodeUwHistoneK562H3k04me3StdZnf4c50c4Raw.merge.bigWig"

	file_signal_prefix_bed <- "~/temp/K562.UWOther.chromHMM.org.signal"
	file_signal_pdf <- "K562.UWother.hist.org.dist.pdf"
	file_signal_make <- "suw_k562/k562"
	file_signal_rdata <- "K562.UWother.hist.org.rdata"
	path_chromHMM <- "suw_k562"
	path_chromHMM2 <- "MS_K562_UWother_seg"
	file.chrom.ret.bed.gz <- "~/temp/MS_K562_UWother_seg.bed.gz"
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

write.signal.bed<-function(df.bed, file.zip)
{
     df.bed <- df.bed[which(df.bed[,4]!=0),]
     write.bed(df.bed, file.zip, compress=TRUE);
}

read.bigwig.fast <- function(file.hist, tb.bed, op="sum" )
{
  interval <- unique(c( seq( 1, NROW(tb.bed)+1, by = 1000*100 ), NROW(tb.bed)+1))

  ret <- do.call("c", mclapply(1:(length(interval)-1), function(x)
  {
    bw <- load.bigWig(file.hist);
    batch_indx<- c( interval[x]:(interval[x+1]-1) ) 
    dat <- bed.region.bpQuery.bigWig( bw, tb.bed[batch_indx, ], op=op);
    unload.bigWig(bw);
    return( c(dat) );
  }, mc.cores=24));
  
  return(ret);
}

tb.chrom <- read.table("/fs/cbsudanko/storage/data/hg19/hg19.chromInfo");
tb.chrom <- tb.chrom[grep("_|chrM|chrY|chrX", tb.chrom[,1], invert=TRUE),]

tb.bed <-do.call("rbind", lapply(1:NROW(tb.chrom), function(i){
   return(data.frame(chr=tb.chrom[i,1], start=seq(1, tb.chrom[i,2]-200, 200), stop=seq(1, tb.chrom[i,2]-200, 200)+200 ));
}));

tb.h3k27ac.sum <- tb.h3k27me3.sum <- tb.h3k36me3.sum <- tb.h3k9me3.sum <- tb.h3k4me1.sum <- tb.h3k4me3.sum <- c();
tb.h3k27ac.max <- tb.h3k27me3.max <- tb.h3k36me3.max <- tb.h3k9me3.max <- tb.h3k4me1.max <- tb.h3k4me3.max <- c();

tb.h3k27ac.sum <- read.bigwig.fast( file_H3k27ac_bw, tb.bed);
tb.h3k27ac.max <- read.bigwig.fast( file_H3k27ac_bw, tb.bed, op="max");

tb.h3k27me3.sum <- read.bigwig.fast( file_H3k27me3_bw, tb.bed);
tb.h3k27me3.max <- read.bigwig.fast( file_H3k27me3_bw, tb.bed, op="max");

tb.h3k36me3.sum <- read.bigwig.fast( file_H3k36me3_bw, tb.bed);
tb.h3k36me3.max <- read.bigwig.fast( file_H3k36me3_bw, tb.bed, op="max");

tb.h3k9me3.sum <- read.bigwig.fast( file_H3k9me3_bw, tb.bed);
tb.h3k9me3.max <- read.bigwig.fast( file_H3k9me3_bw, tb.bed, op="max");

tb.h3k4me1.sum <- read.bigwig.fast( file_H3k4me1_bw, tb.bed);
tb.h3k4me1.max <- read.bigwig.fast( file_H3k4me1_bw, tb.bed, op="max");

tb.h3k4me3.sum <- read.bigwig.fast( file_H3k4me3_bw, tb.bed);
tb.h3k4me3.max <- read.bigwig.fast( file_H3k4me3_bw, tb.bed, op="max");

tb.signal.sum <- data.frame( h3k27ac=tb.h3k27ac.sum, h3k27me3=tb.h3k27me3.sum, h3k36me3=tb.h3k36me3.sum, h3k9me3=tb.h3k9me3.sum, h3k4me1=tb.h3k4me1.sum, h3k4me3=tb.h3k4me3.sum);
tb.signal.max <- data.frame( h3k27ac=tb.h3k27ac.max, h3k27me3=tb.h3k27me3.max, h3k36me3=tb.h3k36me3.max, h3k9me3=tb.h3k9me3.max, h3k4me1=tb.h3k4me1.max, h3k4me3=tb.h3k4me3.max);
save(tb.signal.sum, tb.signal.max, file=file_signal_rdata);
#load(file_signal_rdata);

q1 <- 0.0001
q2 <- 0.05

if(1)
{
tb.h3k27ac.signal <- tb.h3k27me3.signal <- tb.h3k36me3.signal <- tb.h3k9me3.signal <- tb.h3k4me1.signal <- tb.h3k4me3.signal <- rep( 0, NROW(tb.bed));

	F = fitdistr(round(tb.signal.max$h3k27ac), densfun="Poisson")
	Q1 = qpois(1-q1 , F$estimate);
	tb.h3k27ac.signal[ round(tb.signal.max$h3k27ac) >= Q1 & tb.signal.sum$h3k27ac >= quantile(tb.signal.sum$h3k27ac, 1-q2) ] <- 1;

	F = fitdistr(round(tb.signal.max$h3k27me3), densfun="Poisson")
	Q2 = qpois(1-q1, F$estimate);
	tb.h3k27me3.signal[ round(tb.signal.max$h3k27me3) >= Q2 & tb.signal.sum$h3k27me3 >= quantile(tb.signal.max$h3k27me3, 1-q2)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k36me3), densfun="Poisson")
	Q3 = qpois(1- q1, F$estimate);
	tb.h3k36me3.signal[ round(tb.signal.max$h3k36me3) >= Q3 & tb.signal.sum$h3k36me3 >= quantile(tb.signal.sum$h3k36me3, 1-q2)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k9me3), densfun="Poisson")
	Q4 = qpois(1- q1, F$estimate);
	tb.h3k9me3.signal[ round(tb.signal.max$h3k9me3) >= Q4 & tb.signal.sum$h3k9me3 >= quantile(tb.signal.sum$h3k9me3, 1-q2)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k4me1), densfun="Poisson")
	Q5 = qpois(1- q1, F$estimate);
	tb.h3k4me1.signal[ round(tb.signal.max$h3k4me1) >= Q5 & tb.signal.sum$h3k4me1 >= quantile(tb.signal.sum$h3k4me1, 1-q2)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k4me3), densfun="Poisson")
	Q6 = qpois(1-q1, F$estimate);
	tb.h3k4me3.signal[ round(tb.signal.max$h3k4me3) >= Q6 & tb.signal.sum$h3k4me3 >= quantile(tb.signal.sum$h3k4me3, 1-q2)] <- 1;
}
if(0)
{
#state=i 1 0.003282 
#state=i 2 0.01525 
#state=i 3 0.01191 
#state=i 4 0.006077 
#state=i 5 0.02713 
#state=i 6 0.03948 
#state=i 7 0.0788 
#state=i 8 0.05893 
#state=i 9 0.002455 
#state=i 10 0.004538 
#state=i 11 0.01076 
#state=i 12 0.1799 
#state=i 13 0.03522 
#state=i 14 0.06677 
#state=i 15 0.08527 
#state=i 16 0.2259 
#state=i 17 0.07749 
#state=i 18 0.07075 

tb.h3k27ac.signal <- tb.h3k27me3.signal <- tb.h3k36me3.signal <- tb.h3k9me3.signal <- tb.h3k4me1.signal <- tb.h3k4me3.signal <- rep( 0, NROW(tb.bed));


	F = fitdistr(round(tb.signal.sum$h3k27ac), densfun="Poisson")
	Q1 = qpois(1-q1 , F$estimate);
	tb.h3k27ac.signal[ round(tb.signal.sum$h3k27ac) >= Q1  ] <- 1;

	F = fitdistr(round(tb.signal.sum$h3k27me3), densfun="Poisson")
	Q2 = qpois(1-q1, F$estimate);
	tb.h3k27me3.signal[ round(tb.signal.sum$h3k27me3) >= Q2] <- 1;


	F = fitdistr(round(tb.signal.sum$h3k36me3), densfun="Poisson")
	Q3 = qpois(1- q1, F$estimate);
	tb.h3k36me3.signal[ round(tb.signal.sum$h3k36me3) >= Q3] <- 1;


	F = fitdistr(round(tb.signal.sum$h3k9me3), densfun="Poisson")
	Q4 = qpois(1- q1, F$estimate);
	tb.h3k9me3.signal[ round(tb.signal.sum$h3k9me3) >= Q4] <- 1;


	F = fitdistr(round(tb.signal.sum$h3k4me1), densfun="Poisson")
	Q5 = qpois(1- q1, F$estimate);
	tb.h3k4me1.signal[ round(tb.signal.sum$h3k4me1) >= Q5] <- 1;


	F = fitdistr(round(tb.signal.sum$h3k4me3), densfun="Poisson")
	Q6 = qpois(1-q1, F$estimate);
	tb.h3k4me3.signal[ round(tb.signal.sum$h3k4me3) >= Q6] <- 1;
}

tb.signal <- data.frame( H3K27ac=tb.h3k27ac.signal, 
                         H3K27me3=tb.h3k27me3.signal, 
                         H3K36me3=tb.h3k36me3.signal, 
                         H3K4me1=tb.h3k4me1.signal, 
                         H3K4me3=tb.h3k4me3.signal, 
                         H3K9me3=tb.h3k9me3.signal );

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

pdf(file_signal_pdf)
hist(round(tb.h3k27ac.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k27me3.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k36me3.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k9me3.max), breaks=100000, xlim=c(0,100))
hist(round(tb.h3k4me1.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k4me3.max), breaks=10000, xlim=c(0,100))

hist(round(tb.h3k27ac.sum), breaks=10000, xlim=c(0,100*100))
hist(round(tb.h3k27me3.sum), breaks=10000, xlim=c(0,100*100))
hist(round(tb.h3k36me3.sum), breaks=10000, xlim=c(0,100*100))
hist(round(tb.h3k9me3.sum), breaks=100000, xlim=c(0,100*100))
hist(round(tb.h3k4me1.sum), breaks=10000, xlim=c(0,100*100))
hist(round(tb.h3k4me3.sum), breaks=10000, xlim=c(0,100*100))

dev.off();
 



