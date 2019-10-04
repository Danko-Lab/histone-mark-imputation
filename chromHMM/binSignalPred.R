library(bigWig);
library(MASS);
library(parallel)

source("../scripts/hist.param.R");

if(0)
{
	file_H3k27ac_bw  = file.k562.H3k27ac.pred.bw
	file_H3k27me3_bw = file.k562.H3k27me3.pred.bw
	file_H3k36me3_bw = file.k562.H3k36me3.pred.bw
	file_H3k9me3_bw  = file.k562.H3k9me3.pred.bw
	file_H3k4me1_bw  = file.k562.H3k4me1.pred.bw
	file_H3k4me3_bw  = file.k562.H3k4me3.pred.bw

	file_signal_prefix_bed <- "~/temp/K562.chromHMM.pred.filter"
	file_signal_pdf <- "K562.hist.pred.filter.dist.pdf"
	file_signal_make <- "spf_k562/k562"
	file_signal_rdata <- "K562.hist.pred.filter.rdata"
	path_chromHMM <- "spf_k562"
	path_chromHMM2 <- "MS_K562_pred_filter_seg"
	file.chrom.ret.bed.gz <- "~/temp/MS_K562_pred_filter_seg.bed.gz"
}


if(0)
{

	file_H3k27ac_bw  = file.k562.H3k27ac.pred.raw
	file_H3k27me3_bw = file.k562.H3k27me3.pred.raw
	file_H3k36me3_bw = file.k562.H3k36me3.pred.raw
	file_H3k9me3_bw  = file.k562.H3k9me3.pred.raw
	file_H3k4me1_bw  = file.k562.H3k4me1.pred.raw
	file_H3k4me3_bw  = file.k562.H3k4me3.pred.raw

	file_signal_prefix_bed <- "~/temp/K562.chromHMM.pred.signal"
	file_signal_pdf <- "K562.hist.pred.dist.pdf"
	file_signal_make <- "sp_k562/k562"
	file_signal_rdata <- "K562.hist.pred.rdata"
	path_chromHMM <- "sp_k562"
	path_chromHMM2 <- "MS_K562_pred_seg"
	file.chrom.ret.bed.gz <- "~/temp/MS_K562_pred_seg.bed.gz"
}

if(1)
{
	file_H3k27ac_bw  = file.gm.H3k27ac.pred.raw
	file_H3k27me3_bw = file.gm.H3k27me3.pred.raw
	file_H3k36me3_bw = file.gm.H3k36me3.pred.raw
	file_H3k9me3_bw  = file.gm.H3k9me3.pred.raw
	file_H3k4me1_bw  = file.gm.H3k4me1.pred.raw
	file_H3k4me3_bw  = file.gm.H3k4me3.pred.raw

	file_signal_prefix_bed <- "~/temp/GM12878.chromHMM.pred.signal"
	file_signal_pdf <- "GM12878.hist.pred.dist.pdf"
	file_signal_make <- "sp_gm12878/gm12878"
	file_signal_rdata <- "GM12878.hist.pred.rdata"
	path_chromHMM <- "sp_gm12878"
	path_chromHMM2 <- "MS_GM12878_pred_seg"
	file.chrom.ret.bed.gz <- "~/temp/MS_GM12878_pred_seg.bed.gz"
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


tb.h3k27ac.sum <- tb.h3k27me3.sum <- tb.h3k36me3.sum <- tb.h3k9me3.sum <- tb.h3k4me1.sum <- tb.h3k4me3.sum <- c();
tb.h3k27ac.max <- tb.h3k27me3.max <- tb.h3k36me3.max <- tb.h3k9me3.max <- tb.h3k4me1.max <- tb.h3k4me3.max <- c();

tb.h3k27ac.sum <- read.bigwig.fast( file_H3k27ac_bw, tb.bed);
tb.h3k27ac.max <- read.bigwig.fast( file_H3k27ac_bw, tb.bed, op="max", inc=10);

tb.h3k27me3.sum <- read.bigwig.fast( file_H3k27me3_bw, tb.bed);
tb.h3k27me3.max <- read.bigwig.fast( file_H3k27me3_bw, tb.bed, op="max", inc=10);

tb.h3k36me3.sum <- read.bigwig.fast( file_H3k36me3_bw, tb.bed);
tb.h3k36me3.max <- read.bigwig.fast( file_H3k36me3_bw, tb.bed, op="max", inc=10);

tb.h3k9me3.sum <- read.bigwig.fast( file_H3k9me3_bw, tb.bed);
tb.h3k9me3.max <- read.bigwig.fast( file_H3k9me3_bw, tb.bed, op="max", inc=10);

tb.h3k4me1.sum <- read.bigwig.fast( file_H3k4me1_bw, tb.bed);
tb.h3k4me1.max <- read.bigwig.fast( file_H3k4me1_bw, tb.bed, op="max", inc=10);

tb.h3k4me3.sum <- read.bigwig.fast( file_H3k4me3_bw, tb.bed);
tb.h3k4me3.max <- read.bigwig.fast( file_H3k4me3_bw, tb.bed, op="max", inc=10);

tb.signal.sum <- data.frame( h3k27ac=tb.h3k27ac.sum, h3k27me3=tb.h3k27me3.sum, h3k36me3=tb.h3k36me3.sum, h3k9me3=tb.h3k9me3.sum, h3k4me1=tb.h3k4me1.sum, h3k4me3=tb.h3k4me3.sum);
tb.signal.max <- data.frame( h3k27ac=tb.h3k27ac.max, h3k27me3=tb.h3k27me3.max, h3k36me3=tb.h3k36me3.max, h3k9me3=tb.h3k9me3.max, h3k4me1=tb.h3k4me1.max, h3k4me3=tb.h3k4me3.max);
save(tb.signal.sum, tb.signal.max, file=file_signal_rdata);
#load(file=file_signal_rdata);

q1 <- 0.0001
q2 <- 0.05

tb.h3k27ac.signal <- tb.h3k27me3.signal <- tb.h3k36me3.signal <- tb.h3k9me3.signal <- tb.h3k4me1.signal <- tb.h3k4me3.signal <- rep( 0, NROW(tb.bed));
if(1)
{
	F = fitdistr(tb.signal.max$h3k27ac, densfun="normal")
	Q1 = qpois(1-q1 ,  F$estimate[1], F$estimate[2]);
	tb.h3k27ac.signal[ tb.signal.max$h3k27ac >= Q1 & tb.signal.sum$h3k27ac >= quantile(tb.signal.sum$h3k27ac, 1-q2) ] <- 1;

	F = fitdistr(tb.signal.max$h3k27me3, densfun="normal")
	Q2 = qpois(1-q1,  F$estimate[1], F$estimate[2]);
	tb.h3k27me3.signal[ tb.signal.max$h3k27me3 >= Q2 & tb.signal.sum$h3k27me3 >= quantile(tb.signal.max$h3k27me3, 1-q2)] <- 1;


	F = fitdistr(tb.signal.max$h3k36me3, densfun="normal")
	Q3 = qpois(1- q1,  F$estimate[1], F$estimate[2]);
	tb.h3k36me3.signal[tb.signal.max$h3k36me3 >= Q3 & tb.signal.sum$h3k36me3 >= quantile(tb.signal.sum$h3k36me3, 1-q2)] <- 1;


	F = fitdistr(tb.signal.max$h3k9me3, densfun="normal")
	Q4 = qpois(1- q1,  F$estimate[1], F$estimate[2]);
	tb.h3k9me3.signal[ tb.signal.max$h3k9me3 >= Q4 & tb.signal.sum$h3k9me3 >= quantile(tb.signal.sum$h3k9me3, 1-q2)] <- 1;


	F = fitdistr(tb.signal.max$h3k4me1, densfun="normal")
	Q5 = qpois(1- q1,  F$estimate[1], F$estimate[2]);
	tb.h3k4me1.signal[ tb.signal.max$h3k4me1 >= Q5 & tb.signal.sum$h3k4me1 >= quantile(tb.signal.sum$h3k4me1, 1-q2)] <- 1;


	F = fitdistr(tb.signal.max$h3k4me3, densfun="normal")
	Q6 = qpois(1-q1,  F$estimate[1], F$estimate[2]);
	tb.h3k4me3.signal[ tb.signal.max$h3k4me3 >= Q6 & tb.signal.sum$h3k4me3 >= quantile(tb.signal.sum$h3k4me3, 1-q2)] <- 1;
}


if(1)
{
q1 <- 0.0001
q2 <- 0.05


	F = fitdistr(round(tb.signal.max$h3k27ac), densfun="Poisson")
	Q1 = qpois(1-q1 ,  F$estimate);
	tb.h3k27ac.signal[ tb.signal.max$h3k27ac >= Q1 & tb.signal.sum$h3k27ac >= quantile(tb.signal.sum$h3k27ac, 1-q2) ] <- 1;

	F = fitdistr(round(tb.signal.max$h3k27me3), densfun="Poisson")
	Q2 = qpois(1-0.001,  F$estimate);
	tb.h3k27me3.signal[ tb.signal.max$h3k27me3 >= Q2 & tb.signal.sum$h3k27me3 >= quantile(tb.signal.max$h3k27me3, 0.9)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k36me3), densfun="Poisson")
	Q3 = qpois(1- q1,  F$estimate);
	tb.h3k36me3.signal[tb.signal.max$h3k36me3 >= Q3 & tb.signal.sum$h3k36me3 >= quantile(tb.signal.sum$h3k36me3, 1-q2)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k9me3), densfun="Poisson")
	Q4 = qpois(1- 0.001,  F$estimate);
	tb.h3k9me3.signal[ tb.signal.max$h3k9me3 >= Q4 & tb.signal.sum$h3k9me3 >= quantile(tb.signal.sum$h3k9me3, 0.9)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k4me1), densfun="Poisson")
	Q5 = qpois(1- q1,  F$estimate);
	tb.h3k4me1.signal[ tb.signal.max$h3k4me1 >= Q5 & tb.signal.sum$h3k4me1 >= quantile(tb.signal.sum$h3k4me1, 1-q2)] <- 1;


	F = fitdistr(round(tb.signal.max$h3k4me3), densfun="Poisson")
	Q6 = qpois(1-q1,  F$estimate);
	tb.h3k4me3.signal[ tb.signal.max$h3k4me3 >= Q6 & tb.signal.sum$h3k4me3 >= quantile(tb.signal.sum$h3k4me3, 1-q2)] <- 1;
}


if(0)
{
	F = fitdistr( tb.h3k27ac.sum, densfun="normal")
	Q1 = qnorm(1-q1 , F$estimate[1], F$estimate[2]);
	tb.h3k27ac.signal[ tb.h3k27ac.sum >= Q1 ] <- 1;

	F = fitdistr(tb.h3k27me3.sum, densfun="normal")
	Q2 = qnorm(1-q1, F$estimate[1], F$estimate[2]);
	tb.h3k27me3.signal[ tb.h3k27me3.sum >= Q2] <- 1;


	F = fitdistr(tb.h3k36me3.sum, densfun="normal")
	Q3 = qnorm(1- q1, F$estimate[1], F$estimate[2]);
	tb.h3k36me3.signal[ tb.h3k36me3.sum >= Q3] <- 1;


	F = fitdistr(tb.h3k9me3.sum, densfun="normal")
	Q4 = qnorm(1- q1, F$estimate[1], F$estimate[2]);
	tb.h3k9me3.signal[ tb.h3k9me3.sum >= Q4] <- 1;


	F = fitdistr(tb.h3k4me1.sum, densfun="normal")
	Q5 = qnorm(1- q1, F$estimate[1], F$estimate[2]);
	tb.h3k4me1.signal[ tb.h3k4me1.sum >= Q5] <- 1;


	F = fitdistr(tb.h3k4me3.sum, densfun="normal")
	Q6 = qnorm(1-q1, F$estimate[1], F$estimate[2]);
	tb.h3k4me3.signal[ tb.h3k4me3.sum >= Q6] <- 1;
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

pdf(file_signal_pdf);
hist(round(tb.h3k27ac.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k27me3.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k36me3.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k9me3.max), breaks=100000, xlim=c(0,100))
hist(round(tb.h3k4me1.max), breaks=10000, xlim=c(0,100))
hist(round(tb.h3k4me3.max), breaks=10000, xlim=c(0,100))

hist( tb.h3k27ac.sum, breaks=10000, xlim=c(0,100))
hist( tb.h3k27me3.sum, breaks=10000, xlim=c(0,100))
hist( tb.h3k36me3.sum, breaks=10000, xlim=c(0,100))
hist( tb.h3k9me3.sum, breaks=100000, xlim=c(0,100))
hist( tb.h3k4me1.sum, breaks=10000, xlim=c(0,100))
hist( tb.h3k4me3.sum, breaks=10000, xlim=c(0,100))

dev.off();
 
