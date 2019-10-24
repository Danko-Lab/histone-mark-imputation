#/programs/R-3.4.2/bin/R
source("/workdir/zw355/proj/prj15-histone/analyze-dnase-imputation/draw-base.R");

g_noplot <<- FALSE;

setwd("/workdir/zw355/proj/prj15-histone/figures");

file.temp.black  <- write.temp.bed(read.table(file.hg19.black.bed)[,c(1:3)], compress=FALSE )
tb.unmap.bed <- read.table(file.hg19.unmap.bed)[,c(1:3)];
tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
rm(tb.unmap.bed);

zoom_rate=1;

setwd("/workdir/zw355/proj/prj15-histone/analyze-dnase-imputation")

file.removed.reg <<- write.temp.bed(read.table(pipe(paste("cat ", file.temp.black, file.temp.unmap , " | sort-bed - | bedtools merge -i - "))));

if(1)
{
  draw_cor_signal( "dNase.K562.1k",    file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", range=1000, xlim=c(0,50), ylim=c(0,50), out.bed.file= "dnase.K562.1k-data.bed");
  draw_cor_signal( "dNase.K562.100",    file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", range=100, xlim=c(0,10), ylim=c(0,10), out.bed.file= "dnase.K562.100-data.bed");
  draw_cor_signal( "dNase.K562.peaks.1k", file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", range=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.dnase.peakcalling, out.bed.file= "dnase.K562.peaks.1k-data.bed");
  draw_cor_signal( "dNase.K562.peaks.100", file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", range=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.dnase.peakcalling, out.bed.file= "dnase.K562.peaks.100-data.bed" );  
  draw_cor_signal( "raw.dNase.K562.peaks.CTCF.1k",    file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", range=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.K562.CTCF.bed, out.bed.file= "dnase.K562.peaks.CTCF.1k-data.bed" );
  draw_cor_signal( "raw.dNase.K562.peaks.CTCF.100",    file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", range=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.K562.CTCF.bed, out.bed.file= "dnase.K562.peaks.CTCF.100-data.bed"  );
  draw_cor_signal( "dNase.K562.peaks.K27ac.1k", file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", range=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.k562.H3k27ac.peak, out.bed.file= "dnase.K562.peaks.K27ac.1k-data.bed");
  draw_cor_signal( "dNase.K562.peaks.K27ac.100", file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", range=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.k562.H3k27ac.peak, out.bed.file= "dnase.K562.peaks.K27ac.100-data.bed" );
}

## Conservative peaks for TFs.
if(0)
{
  peakpath <- "~/cbsudanko/projects/TFBindingSVR/dTOX/k562/chipseq_peak/"
  narrowPeaks <- dir(peakpath, pattern="conserv.bed.gz")
  
    
  png("CompleteData.png", width=2500, height=2500)
  par(mfrow=c(8, 8))

  for(np in narrowPeaks) {
   
    file.bed.path = paste(peakpath, np, sep="")
    #draw_cor_signal( paste(np, "dNase.K562.peaks.1k", sep=""), file.dnase.bw, file.dNase.pred.bw, NULL, chr="chr22", range=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.bed.path );
    draw_cor_signal( paste(np, "dNase.K562.peaks.100", sep=""), file.dnase.bw, file.dNase.pred.bw, NULL, chr="chr22", range=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.bed.path );
  }
  
  dev.off()
}

## SYDH narrowPeaks
if(0)
{
  peakpath <- "~/cbsudanko/data/hg19/k562/sydh_tfs/"
  narrowPeaks <- dir(peakpath, pattern="narrowPeak")
  for(np in narrowPeaks) {
    file.bed.path = paste(peakpath, np, sep="")
    draw_cor_signal( paste(np, "dNase.K562.peaks.100", sep=""), file.dnase.bw, file.dNase.pred.bw, NULL, chr="chr22", range=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.bed.path );
  }
}

## HAIB broadPeaks
if(0)
{
  peakpath <- "~/cbsudanko/data/hg19/k562/haib_tfs/"
  narrowPeaks <- dir(peakpath, pattern="broadPeak")
  for(np in narrowPeaks) {
    file.bed.path = paste(peakpath, np, sep="")
    draw_cor_signal( paste(np, "dNase.K562.peaks.100", sep=""), file.dnase.bw, file.dNase.pred.bw, NULL, chr="chr22", range=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.bed.path );
  }
}
  
