library(bigWig)
library(dREG)
library(data.table);
library(snowfall);
library(minerva);
library(parallel);

#source("../script/hist.svm.main.R");
#source("../script/hist.svm.com.R");
source("../script/hist.param.R");

calculate_cor <- function (y.org, y.pred , range=2000 )
{
   if(range!=10)
   {
       n.per.range <- round( range/10 );

       y.range.org <-  unlist(lapply(1:floor(NROW(y.org)/n.per.range),  function(i) {sum(y.org [(i-1)*n.per.range+c(1:n.per.range)])}));
       y.range.pred <- unlist(lapply(1:floor(NROW(y.pred)/n.per.range), function(i) {sum(y.pred[(i-1)*n.per.range+c(1:n.per.range)])}));
   }
   else
   {
        y.range.org <- y.org;
        y.range.pred <- y.pred;
   }

   return(list(r.pearson = cor( y.range.org, y.range.pred, method = "pearson" ),
      r.spearman = cor( y.range.org, y.range.pred, method = "spearman" ),
      r.mad = mad( y.range.org -  y.range.pred ) )  );
}


get_histone_read<-function(bedA, file.histone, ref.bed=NULL, method="mean", normalized=TRUE)
{
   bedA <- bedA[,c(1:3)];
   na.idx <- which( is.na(bedA[,1]) |  is.na(bedA[,2]) |  is.na(bedA[,3])  )
   if(length(na.idx)>0)
         stop("NA in bedA\n");

   bw.hist <- load.bigWig(file.histone);

   y <- c();
   for(i in 1:ceiling(NROW(bedA)/1000000))
   {
       start <- 1+(i-1)*1000000;
       end <- if(i*1000000>NROW(bedA)) NROW(bedA) else i*1000000;
      y <- try( c(y, unlist(bed.region.bpQuery.bigWig( bw=bw.hist, bed= bedA[start:end,]) )) );
      if(class(y)=="try-error")
      {
          browser();
         return(NULL);
         }
   }

   if(!is.null(ref.bed))
   {
      file.map <- tempfile(fileext=".bed")
      write.table( cbind(bedA, y), file=file.map, quote=F, row.names=F, col.names=F, sep="\t");
      file.ref <- tempfile(fileext=".bed")
      write.table( ref.bed, file=file.map, quote=F, row.names=F, col.names=F, sep="\t");
      tb <- read.table(pipe( paste("bedmap --echo --mean ", file.ref,  file.map )));
 browser();
      y <- tb[,4];

      if(normalized)
         y <- y / (bw.hist$basesCovered * bw.hist$mean ) * 1e9
   }

   unload.bigWig(bw.hist);

   return(y);
}

compare_signal <- function( out.prefix, file.bed, file.bw1, file.bw2=NULL, ref.gene.bed=NULL, scatter=TRUE, scatter.log=NULL, range=2000, est.range=c(10, 20, 50, 100, seq(200, 10000, 200)) )
{
    cat("source 1 file=", file.bed, "\n");
    cat("source 2 file=", file.bw1, "\n");
    cat("source 3 file=", file.bw2, "\n");

    tb.s1 <- read.table(file.bed);
    if(is.null(file.bw2))
    {
        y.org <- tb.s1[,4];
        y.bw1 <- get_histone_read( tb.s1[, c(1:3)], file.bw1, ref.gene.bed);
    }
    else
    {
        y.bw1 <- get_histone_read( tb.s1[, c(1:3)], file.bw1, ref.gene.bed);
        y.org <- get_histone_read( tb.s1[, c(1:3)], file.bw2, ref.gene.bed);
    }

    if(scatter)
    {
        r <- rbindlist( mclapply(est.range, function(x) { calculate_cor( y.org, y.bw1, x) }, mc.cores=12  ) )

        pdf(paste(out.prefix, "trend", "pdf", sep="."), width=6, height=6);
        plot(est.range, unlist(r[,1]), type="l", col="red", cex=0.6, xlim=range(est.range), ylim=c(-0.2,1) )
        lines(est.range, unlist(r[,2]), cex=0.6, col="blue")
        lines(est.range, unlist(r[,3]/max(r[,3])), cex=0.6, col="black", lty="22");
        
        legend("topleft", c("pearson", "spearman", "mad"), col=c("red", "blue", "black"), lty=c("solid", "solid", "22") );
        dev.off();
    }

    #r.mic1 = mine( y.org,  y.bw1, n.cores=12, alpha=0.4 )$MIC;
    r.mic1 = 0;

    n.per.range <- round( range/10 );
    y.range.pred <-  unlist(lapply(1:floor(NROW(y.org)/n.per.range),  function(i) {mean(y.org [(i-1)*n.per.range+c(1:n.per.range)])}));
    y.range.org <- unlist(lapply(1:floor(NROW(y.bw1)/n.per.range), function(i) {mean(y.bw1[ (i-1)*n.per.range+c(1:n.per.range) ])}));

    write.table(cbind(y.range.org,     y.range.pred), file=paste(out.prefix, ".org.pred.2k.tab", sep=""), quote=F, row.names=F, col.names=F, sep="\t");

    r.cor1 = cor( y.range.org, y.range.pred, method = "pearson" );
    r.cor2 = cor( y.range.org, y.range.pred, method = "spearman" );
    r.mad = mad( y.range.org -  y.range.pred );
    r.gmic = mine( y.range.org,  y.range.pred, n.cores=12, )$MIC;

    cat("COR=", round(r.cor1,3),  round(r.cor2,3), round(r.mad,3), "MIC[g|1]=", round(r.gmic,3), round(r.mic1,3), "\n" );

    if(scatter)
    {
        source("/home/zw355/src/Rplot/denScatter.R");
        png(paste(out.prefix, "all", "png", sep="."), width=900, height=900);
        densScatterplot( y.range.org, y.range.pred, main=paste(out.prefix, "(r=", round(r.cor1,2),"/", round(r.cor2,2), ")"), xlab="histone", ylab="predict");
        dev.off();

        png(paste(out.prefix, "all", "log", "png", sep="."), width=900, height=900);
        densScatterplot( y.range.org, y.range.pred, uselog=TRUE, main=paste(out.prefix, "(r=", round(r.cor1,2),"/", round(r.cor2,2), ")"), xlab="log(histone)", ylab="log(predict)");
        dev.off();
    }
}


## Read in refseq genes.
refGene <- read.table("/fs/cbsudanko/storage/projects/mcf7tamres/annotations/refGene.bed.gz")
refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene <- refGene[(refGene$V3-refGene$V2)>5000,]
bodies <- refGene
bodies <- bodies [bodies$V1!="chrY",]
bodies$V2[bodies$V6 == "+"] <-bodies$V2[bodies$V6 == "+"]+1000
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-1000

if(1)
{
    compare_signal("H3k27ac.S1.V3.G1.chr22",   "../k562/H3k27ac.S1.V3.G1_chr22.bed.gz",    file.k562.H3k27ac.bw,   range=2000, scatter=FALSE );
    compare_signal("H3k122ac.S1.V3.G1.chr22",  "../k562/H3k122ac.S1.V3.G1_chr22.bed.gz",   file.k562.H3k122ac.bw,  range=2000, scatter=FALSE );
    compare_signal("H3K27ac.S1.exl21.G1.chr22","../k562/H3K27ac.S1.exl21.G1_chr22.bed.gz", file.k562.H3k27ac.bw,   range=2000, scatter=FALSE  );
    compare_signal("H3K27ac.S1.G1.chr22" ,     "../k562/H3K27ac.S1.G1_chr22.bed.gz",       file.k562.H3k27ac.bw,   range=2000, scatter=FALSE  );
    compare_signal("H3k4me1.S1.V2.G1.chr22",   "../k562/H3k4me1.S1.V2.G1_chr22.bed.gz",    file.k562.H3k4me1.bw,   range=2000, scatter=FALSE  );
    compare_signal("H3k4me2.S1.V2.G1.chr22",   "../k562/H3k4me2.S1.V2.G1_chr22.bed.gz",    file.k562.H3k4me2.bw,   range=2000, scatter=FALSE  );
    compare_signal("H3k9ac.S1.V2.G1.chr22",    "../k562/H3k9ac.S1.V2.G1_chr22.bed.gz",     file.k562.H3k9ac.bw,    range=2000, scatter=FALSE  );
    compare_signal("H4k20me1.S1.V3.G1.chr22",  "../k562/H4k20me1.S1.V3.G1_chr22.bed.gz",   file.k562.H4k20me1.bw,  range=2000, scatter=FALSE  );
    compare_signal("H3k9me3.S1.V2.G1.chr22",   "../k562/H3k9me3.S1.V2.G1_chr22.bed.gz",    file.k562.H3k9me3.bw,   range=2000, scatter=FALSE  );
    compare_signal("H3k4me3.S1.V3.G1.chr22",   "../k562/H3k4me3.S1.V3.G1_chr22.bed.gz",    file.k562.H3k4me3.bw,   range=2000, scatter=FALSE  );
    compare_signal("H3k27me3.S1.V3.G1.chr22",  "../k562/H3k27me3.S1.V3.G1_chr22.bed.gz",   file.k562.H3k27me3.bw,  range=2000, scatter=FALSE  );
    compare_signal("H3k36me3.S1.V3.G1.chr22",  "../k562/H3k36me3.S1.V3.G1_chr22.bed.gz",   file.k562.H3k36me3.bw,  range=2000, scatter=FALSE  );

    compare_signal("H3k27ac.S1.V3.GM12878.chr22",  "../gm12878/H3k27ac.S1.V3.GM12878_chr22.bed.gz",  file.gm.H3k27ac.bigw,  range=2000, scatter=FALSE );
    compare_signal("H3k27me3.S1.V3.GM12878.chr22", "../gm12878/H3k27me3.S1.V3.GM12878_chr22.bed.gz", file.gm.H3k27me3.bigw, range=2000, scatter=FALSE );
    compare_signal("H3k9ac.S1.V2.GM12878.chr22",   "../gm12878/H3k9ac.S1.V2.GM12878_chr22.bed.gz",   file.gm.H3k9ac.bigw,   range=2000, scatter=FALSE );
    compare_signal("H3k9me3.S1.V2.GM12878.chr22",  "../gm12878/H3k9me3.S1.V2.GM12878_chr22.bed.gz",  file.gm.H3k9me3.bigw,  range=2000, scatter=FALSE );
    compare_signal("H3k4me1.S1.V2.GM12878.chr22",  "../gm12878/H3k4me1.S1.V2.GM12878_chr22.bed.gz",  file.gm.H3k4me1.bigw,  range=2000, scatter=FALSE );
    compare_signal("H3k4me2.S1.V2.GM12878.chr22",  "../gm12878/H3k4me2.S1.V2.GM12878_chr22.bed.gz",  file.gm.H3k4me2.bigw,  range=2000, scatter=FALSE );
    compare_signal("H3k4me3.S1.V3.GM12878.chr22",  "../gm12878/H3k4me3.S1.V3.GM12878_chr22.bed.gz",  file.gm.H3k4me3.bigw,  range=2000, scatter=FALSE );
    compare_signal("H4k20me1.S1.V3.GM12878.chr22", "../gm12878/H4k20me1.S1.V3.GM12878_chr22.bed.gz", file.gm.H4k20me1.bigw, range=2000, scatter=FALSE );
    compare_signal("H3k36me3.S1.V3.GM12878.chr22", "../gm12878/H3k36me3.S1.V3.GM12878_chr22.bed.gz", file.gm.H3k36me3.bigw, range=2000, scatter=FALSE );

    #compare_signal("H3k27ac.S1.V3.G1.gene.chr22.2k",   "../k562/H3k27ac.S1.V3.G1_chr22.bed.gz",    histpath(file.H3k27ac.bw),   ref.gene.bed=bodies[bodies$V1=="chr22",], range=2000 );

}


compare_signal("GSM2877103_ChIP-seq_K562_H3K27ac_chr22",
    file.bed="../k562/H3k27ac.S1.V3.G1_chr22.bed.gz",
    "../rep_expm/GSM2877103_ChIP-seq_K562_H3K27ac_rep1.bw",
    "../rep_expm/GSM2877104_ChIP-seq_K562_H3K27ac_rep2.bw",
    range=2000);

compare_signal("GSM2877111_ChIP-seq_K562_Non-hub_KO_H3K27ac_chr22",
    file.bed="../k562/H3k27ac.S1.V3.G1_chr22.bed.gz",
    "../rep_expm/GSM2877111_ChIP-seq_K562_Non-hub_KO_H3K27ac_rep1.bw",
    "../rep_expm/GSM2877112_ChIP-seq_K562_Non-hub_KO_H3K27ac_rep2.bw",
    range=2000);

compare_signal("GSM2877119_ChIP-seq_K562_Hub_KO_H3K27ac_chr22",
    file.bed="../k562/H3k27ac.S1.V3.G1_chr22.bed.gz",
    "../rep_expm/GSM2877119_ChIP-seq_K562_Hub_KO_H3K27ac_rep1.bw",
    "../rep_expm/GSM2877120_ChIP-seq_K562_Hub_KO_H3K27ac_rep2.bw",
    range=2000);
