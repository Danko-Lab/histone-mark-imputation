source("draw-cor-base.R");

g_noplot <<- FALSE;

file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )
tb.unmap.bed <- read.table(file.unmap.bed)[,c(1:3)];
tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
rm(tb.unmap.bed);

if(1)
{
  file.IceChip.H3k4me3.bw  <- "/local/workdir/zw355/proj/prj15-histone/pred-iceChip/GSK562_AB-12209_H3K4me3_ICeChIP.hg19.bw"
  file.IceChip.H3k4me3.peak <- "/local/workdir/zw355/proj/prj15-histone/pred-iceChip/GSK562_AB-12209_H3K4me3_ICeChIP.hg19.macs2.peak"
  file.IceChip.H3k4me3.pred.raw<- "/local/workdir/zw355/proj/prj15-histone/pred-iceChip/raw.IceChip.H3k4me3.S1.bigWig"

  add_spikes_to_removed_regs( file.IceChip.H3k4me3.peak, file.IceChip.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.IceChip.H3k4me3.S1.10k",  file.IceChip.H3k4me3.bw,  file.IceChip.H3k4me3.pred.raw,    
                  NULL,   chr="chr22", win.size=10000, xlim=c(0,2000), ylim=c(0,2000), file.peaks=NULL);

  file.IceChip.H3k4me3.bw  <- "/local/workdir/zw355/proj/prj15-histone/pred-iceChip/GSK562_AB-12209_H3K4me3_ICeChIP.hg19.bw"
  file.IceChip.H3k4me3.peak <- "/local/workdir/zw355/proj/prj15-histone/pred-iceChip/GSK562_AB-12209_H3K4me3_ICeChIP.hg19.macs2.peak"
  file.IceChip.H3k4me3.pred.raw<- "/local/workdir/zw355/proj/prj15-histone/pred-iceChip/raw.IceChip.H3k4me3.S1.bigWig"

  add_spikes_to_removed_regs( file.IceChip.H3k4me3.peak, file.IceChip.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.IceChip.H3k4me3.S1.10k",  file.IceChip.H3k4me3.bw,  file.IceChip.H3k4me3.pred.raw,    
                  NULL,   chr="chr22", win.size=10000, xlim=c(0,2000), ylim=c(0,2000), file.peaks=NULL);

}
