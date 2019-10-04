source("draw-cor-base.R");
source("draw-cor-param.R");


if(0)
{
  file.removed.reg <<- write.temp.bed(read.table(pipe(paste("cat ", file.temp.black, file.temp.unmap , " | sort-bed - | bedtools merge -i - "))));

  draw_cor_signal( "dNase.K562.10k",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=10000, xlim=c(0,50), ylim=c(0,50) );
  draw_cor_signal( "dNase.K562.1k",    file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=1000, xlim=c(0,50), ylim=c(0,50) );
  draw_cor_signal( "dNase.K562.100",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10) );

  #file.tmp.bed <- tempfile(fileext=".bed");  
  #system(paste("zcat ", file.K562.CTCF.bed, " | cat", file.dnase.peakcalling, " - | sort-bed - | bedtools merge -i - | sort-bed - | bgzip > ", file.tmp.bed))

  draw_cor_signal( "dNase.K562.peaks.10k",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=10000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.dnase.peakcalling );
  draw_cor_signal( "dNase.K562.peaks.1k",    file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.dnase.peakcalling );
  draw_cor_signal( "dNase.K562.peaks.100",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.dnase.peakcalling );

  draw_cor_signal( "dNase.K562.peaks.CTCF.10k",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=10000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "dNase.K562.peaks.CTCF.1k",    file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "dNase.K562.peaks.CTCF.100",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.K562.CTCF.bed );

  draw_cor_signal( "raw.dNase.K562.peaks.10k",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=10000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.dnase.peakcalling );
  draw_cor_signal( "raw.dNase.K562.peaks.1k",    file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.dnase.peakcalling );
  draw_cor_signal( "raw.dNase.K562.peaks.100",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.dnase.peakcalling );

  draw_cor_signal( "raw.dNase.K562.peaks.CTCF.10k",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=10000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "raw.dNase.K562.peaks.CTCF.1k",    file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=1000, xlim=c(0,50), ylim=c(0,50), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "raw.dNase.K562.peaks.CTCF.100",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.K562.CTCF.bed );

  draw_cor_signal( "dNase.K562.100",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10) );
  draw_cor_signal( "dNase.K562.CTCF.100",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "dNase.K562.dNase.peaks.100",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.dnase.peakcalling );
  draw_cor_signal( "dNase.K562.H3k27ac.peaks.100",   file.dnase.bw,     file.dNase.pred.bw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.k562.H3k27ac.peak );


  draw_cor_signal( "raw.dNase.K562.100",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10) );
  draw_cor_signal( "raw.dNase.K562.CTCF.100",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.K562.CTCF.bed );
  draw_cor_signal( "raw.dNase.K562.dNase.peaks.100",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.dnase.peakcalling );
  draw_cor_signal( "raw.dNase.K562.H3k27ac.peaks.100",   file.dnase.bw,     file.dNase.pred.raw,  NULL, chr="chr22", win.size=100, xlim=c(0,10), ylim=c(0,10), file.peaks=file.k562.H3k27ac.peak );
}
