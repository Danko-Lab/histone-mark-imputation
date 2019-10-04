source("draw-base.R");

g_noplot <<- TRUE;

file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )
tb.unmap.bed <- read.table(file.unmap.bed)[,c(1:3)];
tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
rm(tb.unmap.bed);

zoom_rate=1;

if(1)
{
  # G1/chr22/10K
  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )
  r0 <- draw_cor_signal("G1.chr22.H3k27ac.S1.V3.1k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.bw,   file.k562.H3k27ac.ctrl.bw, chr="chr22", range=10000 );
  r1 <- draw_cor_signal("G1.chr22.H3k27ac.S1.V3.1k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.bw,   file.k562.H3k27ac.ctrl.bw, chr="chr22", range=1000 );
  r2 <- draw_cor_signal("G1.chr22.H3k27ac.S1.V3.1k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.bw,   file.k562.H3k27ac.ctrl.bw, chr="chr22", range=100 );
  r3 <- draw_cor_signal("G1.chr22.H3k27ac.S1.V3.1k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.bw,   file.k562.H3k27ac.ctrl.bw, chr="chr22", range=10 );

  df.G1.chr22 <- rbind(r0, r1, r2, r3);
}

