source("draw-cor-base.R");

g_noplot <<- FALSE;

# no black list in mm9
#file.blacklist <- "/fs/cbsudanko/storage/data/mm10/mm10.blacklist.bed.gz"
#file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )

file.temp.black <- NULL
file.temp.unmap <- NULL

if(1)
{
  file.prophaseLZ.H3k4me3.pred.raw<- "../prophase1/LZ.mm10.H3k4me3.S1.V3.bw"
  file.prophaseP.H3k4me3.pred.raw<- "../prophase1/P.mm10.H3k4me3.S1.V3.bw"
  file.prophaseD.H3k4me3.pred.raw<- "../prophase1/D.mm10.H3k4me3.S1.V3.bw"
  
  file.compare.H3k4me3.bw <- "../prophase1/GSM3734414_ZY.H3K4me3.merge.bw"
  file.compare.H3k4me3.peak <- "../prophase1/GSM3734414_ZY.H3K4me3.merge.mm10.peak.bed"
  
  add_spikes_to_removed_regs( file.compare.H3k4me3.peak, file.compare.H3k4me3.bw )
  r2 <- draw_cor_signal("raw.prophaseLZ.H3k4me3.S1.1k",  file.compare.H3k4me3.bw,  file.prophaseLZ.H3k4me3.pred.raw,    
                  NULL,   chr="chr1", win.size=1000, xlim=c(0,20000), ylim=c(0,20000), file.peaks=file.compare.H3k4me3.peak);

  r2 <- draw_cor_signal("raw.prophaseD.H3k4me3.S1.1k",  file.compare.H3k4me3.bw,  file.prophaseD.H3k4me3.pred.raw,    
                  NULL,   chr="chr1", win.size=1000, xlim=c(0,20000), ylim=c(0,20000), file.peaks=file.compare.H3k4me3.peak);

  r2 <- draw_cor_signal("raw.prophaseP.H3k4me3.S1.1k",  file.compare.H3k4me3.bw,  file.prophaseP.H3k4me3.pred.raw,    
                  NULL,   chr="chr1", win.size=1000, xlim=c(0,20000), ylim=c(0,20000), file.peaks=file.compare.H3k4me3.peak);
   
}

