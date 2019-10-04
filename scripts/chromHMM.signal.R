source("hist.param.R");
source("hist.svm.com.R")

library(bigWig)
library(MASS);
library(parallel)
library(tools)

file.hg19.chroInfo <- "/fs/cbsudanko/storage/data/hg19/hg19.chromInfo";

convert_bed_to_bw <- function(file.bed.gz, file.chr.info)
{
    file.bed.tmp <- tempfile(fileext=".bed");
    file.bw.tmp <- tempfile(fileext=".bw");
    #system(paste("zcat  ", file.bed.gz, ">" , file.bed.tmp));
    system(paste("zcat  ", file.bed.gz, "| awk '{ print tolower( substr( $0, 1, 1 ) ) substr( $0, 2 ); }' - | sort-bed - > " , file.bed.tmp));
    system(paste("bedGraphToBigWig ", file.bed.tmp, file.chr.info, file.bw.tmp));
    return(file.bw.tmp)
}

write.signal.bed<-function(df.bed, file.zip)
{
    df.bed <- df.bed[which(df.bed[,4]!=0),]
    write.bed(df.bed, file.zip, compress=TRUE);
}

read.bigwig.fast <- function(file.hist, tb.bed, op="sum" )
{
  interval <- unique(c( seq( 1, NROW(tb.bed)+1, by = 1000*1000 ), NROW(tb.bed)+1))

  ret <- do.call("c", mclapply(1:(length(interval)-1), function(x)
  {
    bw <- load.bigWig(file.hist);
    batch_indx<- c( interval[x]:(interval[x+1]-1) ) 
    dat <- bed.region.bpQuery.bigWig( bw, tb.bed[batch_indx, ], op=op);
    unload.bigWig(bw);
    return( c(dat) );
  }, mc.cores=4));
  
  return(ret);
}


get_chromHMM <- function(file.chrom.prefix, file.H3k27ac.bw,  file.H3k27me3.bw, file.H3k36me3.bw, file.H3k4me1.bw, file.H3k4me3.bw, file.H3k9me3.bw, chr=NULL)
{
    org.file.H3k27ac.bw = file.H3k27ac.bw;
    org.file.H3k27me3.bw = file.H3k27me3.bw; 
    org.file.H3k36me3.bw = file.H3k36me3.bw; 
    org.file.H3k4me1.bw = file.H3k4me1.bw; 
    org.file.H3k4me3.bw = file.H3k4me3.bw; 
    org.file.H3k9me3.bw = file.H3k9me3.bw;

    if(file_ext(org.file.H3k27ac.bw)=="gz")  file.H3k27ac.bw  <- convert_bed_to_bw(file.H3k27ac.bw, file.hg19.chroInfo);
    if(file_ext(org.file.H3k27me3.bw)=="gz") file.H3k27me3.bw <- convert_bed_to_bw(file.H3k27me3.bw, file.hg19.chroInfo);
    if(file_ext(org.file.H3k36me3.bw)=="gz") file.H3k36me3.bw <- convert_bed_to_bw(file.H3k36me3.bw, file.hg19.chroInfo);
    if(file_ext(org.file.H3k4me1.bw)=="gz")  file.H3k4me1.bw  <- convert_bed_to_bw(file.H3k4me1.bw, file.hg19.chroInfo);
    if(file_ext(org.file.H3k4me3.bw)=="gz")  file.H3k4me3.bw  <- convert_bed_to_bw(file.H3k4me3.bw, file.hg19.chroInfo);
    if(file_ext(org.file.H3k9me3.bw)=="gz")  file.H3k9me3.bw  <- convert_bed_to_bw(file.H3k9me3.bw, file.hg19.chroInfo);
    
    tb.chrom <- read.table(file.hg19.chroInfo);
    tb.chrom <- tb.chrom[grep("_|chrM|chrY|chrX", tb.chrom[,1], invert=TRUE),]

    tb.all.bed <-do.call("rbind", lapply(1:NROW(tb.chrom), function(i){
       return(data.frame(chr=tb.chrom[i,1], start=seq(1, tb.chrom[i,2]-200, 200), stop=seq(1, tb.chrom[i,2]-200, 200)+200 ));
    }));

    tb.bed <- bedTools.intersect(tb.all.bed, bed.proseq, options="-wa")
    tb.bed.idx <- match(paste(tb.bed[,1], tb.bed[,2], sep=":") , paste(tb.all.bed[,1], tb.all.bed[,2], sep=":") );

    tb.h3k27ac.sum <- read.bigwig.fast( file.H3k27ac.bw, tb.bed);
    tb.h3k27ac.max <- read.bigwig.fast( file.H3k27ac.bw, tb.bed, op="max");
    tb.h3k27me3.sum <- read.bigwig.fast(file.H3k27me3.bw, tb.bed);
    tb.h3k27me3.max <- read.bigwig.fast(file.H3k27me3.bw, tb.bed, op="max");
    tb.h3k36me3.sum <- read.bigwig.fast(file.H3k36me3.bw, tb.bed);
    tb.h3k36me3.max <- read.bigwig.fast(file.H3k36me3.bw, tb.bed, op="max");
    tb.h3k9me3.sum <- read.bigwig.fast( file.H3k9me3.bw, tb.bed);
    tb.h3k9me3.max <- read.bigwig.fast( file.H3k9me3.bw, tb.bed, op="max");
    tb.h3k4me1.sum <- read.bigwig.fast( file.H3k4me1.bw, tb.bed);
    tb.h3k4me1.max <- read.bigwig.fast( file.H3k4me1.bw, tb.bed, op="max");
    tb.h3k4me3.sum <- read.bigwig.fast( file.H3k4me3.bw, tb.bed);
    tb.h3k4me3.max <- read.bigwig.fast( file.H3k4me3.bw, tb.bed, op="max");

if(1)
{
    pdf(paste0(file.chrom.prefix, ".hist.pdf"))
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
}

    tb.signal.sum <- data.frame(h3k27ac=tb.h3k27ac.sum, h3k27me3=tb.h3k27me3.sum, h3k36me3=tb.h3k36me3.sum, h3k9me3=tb.h3k9me3.sum, h3k4me1=tb.h3k4me1.sum, h3k4me3=tb.h3k4me3.sum)
    tb.signal.max <- data.frame(h3k27ac=tb.h3k27ac.max, h3k27me3=tb.h3k27me3.max, h3k36me3=tb.h3k36me3.max, h3k9me3=tb.h3k9me3.max, h3k4me1=tb.h3k4me1.max, h3k4me3=tb.h3k4me3.max)

    q1 <- 0.001
    q2 <- 0.05
    tb.h3k27ac.signal <- tb.h3k27me3.signal <- tb.h3k36me3.signal <- tb.h3k9me3.signal <- tb.h3k4me1.signal <- tb.h3k4me3.signal <- rep( 2, NROW(tb.bed));
    tb.h3k27ac.all.signal <- tb.h3k27me3.all.signal <- tb.h3k36me3.all.signal <- tb.h3k9me3.all.signal <- tb.h3k4me1.all.signal <- tb.h3k4me3.all.signal <- rep( 2, NROW(tb.all.bed));

if(0)
{
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
else
{
    F = fitdistr(round(tb.signal.sum$h3k27ac), densfun="Poisson")
    Q1 = qpois(1-q1 , F$estimate);
    tb.h3k27ac.signal[ round(tb.signal.sum$h3k27ac) >= Q1  ] <- 1;
    tb.h3k27ac.all.signal[tb.bed.idx] <- tb.h3k27ac.signal;

    F = fitdistr(round(tb.signal.sum$h3k27me3), densfun="Poisson")
    Q2 = qpois(1-q1, F$estimate);
    tb.h3k27me3.signal[ round(tb.signal.sum$h3k27me3) >= Q2] <- 1;
    tb.h3k27me3.all.signal[tb.bed.idx] <- tb.h3k27me3.signal;

    F = fitdistr(round(tb.signal.sum$h3k36me3), densfun="Poisson")
    Q3 = qpois(1- q1, F$estimate);
    tb.h3k36me3.signal[ round(tb.signal.sum$h3k36me3) >= Q3] <- 1;
    tb.h3k36me3.all.signal[tb.bed.idx] <- tb.h3k36me3.signal;

    F = fitdistr(round(tb.signal.sum$h3k9me3), densfun="Poisson")
    Q4 = qpois(1- q1, F$estimate);
    tb.h3k9me3.signal[ round(tb.signal.sum$h3k9me3) >= Q4] <- 1;
    tb.h3k9me3.all.signal[tb.bed.idx] <- tb.h3k9me3.signal;

    F = fitdistr(round(tb.signal.sum$h3k4me1), densfun="Poisson")
    Q5 = qpois(1- q1, F$estimate);
    tb.h3k4me1.signal[ round(tb.signal.sum$h3k4me1) >= Q5] <- 1;
    tb.h3k4me1.all.signal[tb.bed.idx] <- tb.h3k4me1.signal;

    F = fitdistr(round(tb.signal.sum$h3k4me3), densfun="Poisson")
    Q6 = qpois(1-q1, F$estimate);
    tb.h3k4me3.signal[ round(tb.signal.sum$h3k4me3) >= Q6] <- 1;
    tb.h3k4me3.all.signal[tb.bed.idx] <- tb.h3k4me3.signal;
}

if(0)
{
    write.signal.bed(data.frame(tb.all.bed, tb.h3k27ac.all.signal),  paste0(file.chrom.prefix, ".H3K27ac.bed.gz"))
    write.signal.bed(data.frame(tb.all.bed, tb.h3k27me3.all.signal), paste0(file.chrom.prefix, ".h3k27me3.bed.gz"))
    write.signal.bed(data.frame(tb.all.bed, tb.h3k36me3.all.signal), paste0(file.chrom.prefix, ".h3k36me3.bed.gz"))
    write.signal.bed(data.frame(tb.all.bed, tb.h3k9me3.all.signal),  paste0(file.chrom.prefix, ".h3k9me3.bed.gz"))
    write.signal.bed(data.frame(tb.all.bed, tb.h3k4me1.all.signal),  paste0(file.chrom.prefix, ".h3k4me1.bed.gz"))
    write.signal.bed(data.frame(tb.all.bed, tb.h3k4me3.all.signal),  paste0(file.chrom.prefix, ".h3k4me3.bed.gz"))
}

    tb.signal <- data.frame( H3K27ac=tb.h3k27ac.all.signal, 
                             H3K27me3=tb.h3k27me3.all.signal, 
                             H3K36me3=tb.h3k36me3.all.signal, 
                             H3K4me1=tb.h3k4me1.all.signal, 
                             H3K4me3=tb.h3k4me3.all.signal, 
                             H3K9me3=tb.h3k9me3.all.signal );

 
    file.chrom.input      <- paste0(file.chrom.prefix, ".input/")
    dir.create(file.chrom.input);
    file.chrom.output     <- paste0(file.chrom.prefix, ".output/")
    dir.create(file.chrom.output);
    file.chrom.ret.bed.gz <- paste0(file.chrom.prefix, ".chromHMM.bed.gz")  

    L <- lapply(as.character(tb.chrom[,1]), function(chr)
    {
       file.out <- paste(file.chrom.input, "/chrom_", chr, "_binary.txt", sep="");
       idx <- sort(which(tb.all.bed[,1]==chr))
    cat(chr, "=", min(idx), max(idx), NROW(idx), file.out, "\n");
       cat("chrom", chr, "\n", file=file.out, sep="\t", append =FALSE)
       write.table(tb.signal[idx,],  file=file.out, row.names=FALSE, col.names=TRUE, append=TRUE, quote=FALSE, sep="\t")
    });

    system(paste("java -mx6400M -jar ../chromHMM/src/ChromHMM.jar MakeSegmentation ../chromHMM/model_18_core_K27ac.txt", file.chrom.input, file.chrom.output) )

    tb <- read.table(paste(file.chrom.output, "chrom_18_segments.bed", sep="/"))
    tb$V5 <- substring(tb$V4,2);
    for(i in 1:18) cat("state=i", i, sum(as.numeric(tb$V3-tb$V2)[tb$V5==i])/sum(as.numeric(tb$V3-tb$V2)), "\n" )
    write.bed(tb[,c(1,2,3,5)], file=file.chrom.ret.bed.gz, compress=TRUE)

    if(file_ext(org.file.H3k27ac.bw)=="gz")  unlink(file.H3k27ac.bw)
    if(file_ext(org.file.H3k27me3.bw)=="gz") unlink(file.H3k27me3.bw)
    if(file_ext(org.file.H3k36me3.bw)=="gz") unlink(file.H3k36me3.bw)
    if(file_ext(org.file.H3k4me1.bw)=="gz")  unlink(file.H3k4me1.bw)
    if(file_ext(org.file.H3k4me3.bw)=="gz")  unlink(file.H3k4me3.bw)
    if(file_ext(org.file.H3k9me3.bw)=="gz")  unlink(file.H3k9me3.bw)
    
    unlink(file.chrom.input, TRUE)
    unlink(file.chrom.output, TRUE)
}

if(1)
{
   prj.folder <- "/workdir/zw355/proj/prj15-histone/GBM20/"
   folder.chrmm <- "./GBM20.chrmm.001.2nd/"
   dir.create(folder.chrmm);
   tb <- read.table(paste0(prj.folder, "GBM20.file.tab"), stringsAsFactors=F);
   files.bigwig.plus <- tb[,1]
   bed.proseq <- read.table(paste0(prj.folder, "/dreg.call/UMU22.merge.bed.gz"), stringsAsFactors=F);

   for(i in 1:NROW(tb))
   {
      file.H3k27ac.bw  <- paste(prj.folder, "/GBM20.H3k27ac.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
      file.H3k27me3.bw <- paste(prj.folder, "/GBM20.H3k27me3.", basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
      file.H3k36me3.bw <- paste(prj.folder, "/GBM20.H3k36me3.", basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
      file.H3k4me1.bw  <- paste(prj.folder, "/GBM20.H3k4me1.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
      file.H3k4me3.bw  <- paste(prj.folder, "/GBM20.H3k4me3.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
      file.H3k9me3.bw  <- paste(prj.folder, "/GBM20.H3k9me3.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");

      prefix = paste(folder.chrmm, basename(files.bigwig.plus[i]), sep="");
      get_chromHMM(prefix, file.H3k27ac.bw,  file.H3k27me3.bw, file.H3k36me3.bw, file.H3k4me1.bw, file.H3k4me3.bw, file.H3k9me3.bw, chr=NULL)
  }

}


