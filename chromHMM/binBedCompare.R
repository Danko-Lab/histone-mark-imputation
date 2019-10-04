file.k562.H3k27ac.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz";
file.k562.H3k27me3.peak <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k27me3StdAln.bed.gz";
file.k562.H3k36me3.peak <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k36me3StdAln.bed.gz";
file.k562.H3k4me1.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k4me1StdAln.bed.gz";
file.k562.H3k4me3.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k4me3StdAln.bed.gz";
file.k562.H3k9me3.peak  <- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k9me3StdAln.bed.gz";



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
	

tb.chrom <- read.table("/fs/cbsudanko/storage/data/hg19/hg19.chromInfo");
tb.chrom <- tb.chrom[grep("_|chrM|chrY|chrX", tb.chrom[,1], invert=TRUE),]

tb.bed <-do.call("rbind", lapply(1:NROW(tb.chrom), function(i){
   return(data.frame(chr=tb.chrom[i,1], start=seq(1, tb.chrom[i,2]-200, 200), stop=seq(1, tb.chrom[i,2]-200, 200)+200 ));
}));

	
load(file_signal_rdata);

q1 <- 0.0001
q2 <- 0.05


get_jaccard <- function(tb.bed, signal.max, signal.sum, file.peak, q1, q2)
{
    tb.signal <- rep( 0, NROW(tb.bed));
	F = fitdistr(signal.max, densfun="normal")
	Q1 = qpois(1-q1 ,  F$estimate[1], F$estimate[2]);
	tb.signal[ signal.max >= Q1 & signal.sum >= quantile(signal.sum, 1-q2) ] <- 1;

    tmp.bed = tempfile(fileext=".bed");
    df.bed <- data.frame(tb.bed, tb.signal);
    df.bed <- df.bed[which(df.bed[,4]!=0),]
	write.bed(df.bed[,1:3], tmp.bed, compress=FALSE);
    
    tb.peak <- read.table(file.peak)
    tb.intersect <- read.table(pipe(paste("zcat", file.peak, "| awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3}' - | bedtools intersect -a - -b", tmp.bed)), header=FALSE)
    
    #unlink(tmp.bed);
    #return(sum(tb[,1]*tb[,4]))

    unlink(tmp.bed);
    x.merge <- sum(tb.peak$V3-tb.peak$V2) + sum(df.bed$stop-df.bed$start); 
    x.inter <- sum(tb.intersect$V3-tb.intersect$V2)
    return(x.inter/(x.merge-x.inter))

}


df.jac <- c()

for(q1 in c(0.1, 0.01, 0.001, 0.0001, 0.00001) )
{
    jac1 <- get_jaccard(tb.bed, tb.signal.max$h3k27ac, tb.signal.sum$h3k27ac, file.k562.H3k27ac.peak, q1, q2);
    jac2 <- get_jaccard(tb.bed, tb.signal.max$h3k27me3, tb.signal.sum$h3k27me3, file.k562.H3k27me3.peak, q1, q2);
    jac3 <- get_jaccard(tb.bed, tb.signal.max$h3k36me3, tb.signal.sum$h3k36me3, file.k562.H3k36me3.peak, q1, q2);
    jac4 <- get_jaccard(tb.bed, tb.signal.max$h3k9me3, tb.signal.sum$h3k9me3, file.k562.H3k9me3.peak, q1, q2);
    jac5 <- get_jaccard(tb.bed, tb.signal.max$h3k4me3, tb.signal.sum$h3k4me3, file.k562.H3k4me3.peak, q1, q2);
    jac6 <- get_jaccard(tb.bed, tb.signal.max$h3k4me1, tb.signal.sum$h3k4me1, file.k562.H3k4me1.peak, q1, q2);
   
    df.jac <- rbind( df.jac , c(q1, jac1, jac2, jac3, jac4, jac5, jac6) )
}
