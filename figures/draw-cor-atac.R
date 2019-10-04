source("draw-cor-base.R");
source("draw-cor-param.R");




if(0)
{
  file.removed.reg <<- write.temp.bed(read.table(pipe(paste("cat ", file.temp.black, file.temp.unmap , " | sort-bed - | bedtools merge -i - "))));
  draw_cor_signal( "ATAC.K562.10k",   file.k562.ATAC.bw,     file.k562.ATAC.pred.bw,  NULL, chr="chr22", win.size=10000, xlim=c(0,6000), ylim=c(0,6000) );

  draw_cor_signal( "raw.ATAC.K562.10k",   file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=10000, xlim=c(0, 6000), ylim=c(0,6000) );
  draw_cor_signal( "raw.ATAC.K562.1k",    file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=1000, xlim=c(0,600), ylim=c(0,600) );
  draw_cor_signal( "raw.ATAC.K562.100",   file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,60), ylim=c(0,60) );

  draw_cor_signal( "raw.ATAC.K562.peaks.10k",   file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=10000, xlim=c(0, 6000), ylim=c(0,6000), file.peaks=file.k562.ATAC.peak );
  draw_cor_signal( "raw.ATAC.K562.peaks.1k",    file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=1000, xlim=c(0,600), ylim=c(0,600), file.peaks=file.k562.ATAC.peak );
  draw_cor_signal( "raw.ATAC.K562.peaks.100",   file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,60), ylim=c(0,60), file.peaks=file.k562.ATAC.peak );

  draw_cor_signal( "raw.ATAC.K562.peaks.CTCF.10k",   file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=10000, xlim=c(0, 6000), ylim=c(0,6000), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "raw.ATAC.K562.peaks.CTCF.1k",    file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=1000, xlim=c(0,600), ylim=c(0,600), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "raw.ATAC.K562.peaks.CTCF.100",   file.k562.ATAC.bw,     file.k562.ATAC.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,60), ylim=c(0,60), file.peaks=file.K562.CTCF.bed );

}

