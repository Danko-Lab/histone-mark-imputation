source("draw-cor-base.R");

g_noplot <<- FALSE;
file.blacklist <- "/fs/cbsudanko/storage/data/mm10/mm10.blacklist.bed.gz"
file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )

if(0)
{
  file.mESC.H3k27me3.bw <- "../pred-mm9/mm9.GSM307619_ES.H3K27me3.bw"
  file.mESC.H3k27me3.peak <- "../pred-mm9/mm9.GSM307619_ES.H3K27me3.peak"
  file.mESC.H3k27me3.pred.raw<- "../pred-mm9/raw.MM9.mEsc.H3k27me3.S1.pred.bw"

  file.mESC.H3k4me3.bw <- "../pred-mm9/mm9.GSM307618_ES.H3K4me3.bw"
  file.mESC.H3k4me3.peak <- "../pred-mm9/mm9.GSM307618_ES.H3K4me3.peak"
  file.mESC.H3k4me3.pred.raw<- "../pred-mm9/raw.MM9.mEsc.H3k4me3.S1.pred.bw"
  
  add_spikes_to_removed_regs( file.mESC.H3k27me3.peak, file.mESC.H3k27me3.bw )
  r1 <- draw_cor_signal("raw.mESC.H3k27me3.S1.1k",  file.mESC.H3k27me3.bw,  file.mESC.H3k27me3.pred.raw,    
                  NULL,   chr="chr1", win.size=1000, xlim=c(0,2000), ylim=c(0,2000), file.peaks=NULL);
  
  add_spikes_to_removed_regs( file.mESC.H3k27me3.peak, file.mESC.H3k27me3.bw )
  r1 <- draw_cor_signal("raw.mESC.H3k27me3.S1.10k",  file.mESC.H3k27me3.bw,  file.mESC.H3k27me3.pred.raw,    
                  NULL,   chr="chr1", win.size=10000, xlim=c(0,15000), ylim=c(0,15000), file.peaks=NULL);

  add_spikes_to_removed_regs( file.mESC.H3k4me3.peak, file.mESC.H3k4me3.bw )
  r2 <- draw_cor_signal("raw.mESC.H3k4me3.S1.1k",  file.mESC.H3k4me3.bw,  file.mESC.H3k4me3.pred.raw,    
                  NULL,   chr="chr1", win.size=1000, xlim=c(0,2000), ylim=c(0,2000), file.peaks=NULL);
    
  add_spikes_to_removed_regs( file.mESC.H3k4me3.peak, file.mESC.H3k4me3.bw )
  r2 <- draw_cor_signal("raw.mESC.H3k4me3.S1.10k",  file.mESC.H3k4me3.bw,  file.mESC.H3k4me3.pred.raw,    
                  NULL,   chr="chr1", win.size=10000, xlim=c(0,15000), ylim=c(0,15000), file.peaks=NULL);
}
