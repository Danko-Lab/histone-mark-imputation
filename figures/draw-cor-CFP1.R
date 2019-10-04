source("draw-cor-base.R");

g_noplot <<- FALSE;
file.blacklist <- "/fs/cbsudanko/storage/data/mm10/mm10.blacklist.bed.gz"
file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )
file.temp.unmap <- NULL;

file.CFP1.H3k4me3.peak  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.merge.all.narrowPeak.bed.gz"

file.CFP1.wt.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.wt.merge.bw"
file.CFP1.ko.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.ko.1.bw"
file.CFP1.C169A.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.C169A.merge.bw"
file.CFP1.wt_resuce.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.wt_rescue.merge.bw"

file.CFP1.wt.H3k4me3.pred.raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-cfp1/raw.CFP1.wt.H3k4me3.S1.peak.bw"
file.CFP1.ko.H3k4me3.pred.raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-cfp1/raw.CFP1.ko.H3k4me3.S1.peak.bw"
file.CFP1.C169A.H3k4me3.pred.raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-cfp1/raw.CFP1.C169A.H3k4me3.S1.peak.bw"
file.CFP1.wt_resuce.H3k4me3.pred.raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-cfp1/raw.CFP1.wt_resuce.H3k4me3.S1.peak.bw"

if(0)
{
  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.wt.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.wt.H3k4me3.S1.10k",  file.CFP1.wt.H3k4me3.bw,  file.CFP1.wt.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=10000, xlim=c(0,1200000), ylim=c(0,1200000), file.peaks=file.CFP1.H3k4me3.peak);

  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.ko.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.ko.H3k4me3.S1.10k",  file.CFP1.ko.H3k4me3.bw,  file.CFP1.ko.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=10000, xlim=c(0,1200000), ylim=c(0,1200000), file.peaks=file.CFP1.H3k4me3.peak);

  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.C169A.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.C169A.H3k4me3.S1.10k",  file.CFP1.C169A.H3k4me3.bw,  file.CFP1.C169A.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=10000, xlim=c(0,1200000), ylim=c(0,1200000), file.peaks=file.CFP1.H3k4me3.peak);


  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.wt_resuce.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.wt_resuce.H3k4me3.S1.10k",  file.CFP1.wt_resuce.H3k4me3.bw,  file.CFP1.wt_resuce.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=10000, xlim=c(0,1200000), ylim=c(0,1200000), file.peaks=file.CFP1.H3k4me3.peak);
}

if(1)
{
  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.wt.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.wt.H3k4me3.S1.1k",  file.CFP1.wt.H3k4me3.bw,  file.CFP1.wt.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=1000, xlim=c(0,120000), ylim=c(0,120000), file.peaks=file.CFP1.H3k4me3.peak);

  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.ko.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.ko.H3k4me3.S1.1k",  file.CFP1.ko.H3k4me3.bw,  file.CFP1.ko.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=1000, xlim=c(0,120000), ylim=c(0,120000), file.peaks=file.CFP1.H3k4me3.peak);

  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.C169A.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.C169A.H3k4me3.S1.1k",  file.CFP1.C169A.H3k4me3.bw,  file.CFP1.C169A.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=1000, xlim=c(0,120000), ylim=c(0,120000), file.peaks=file.CFP1.H3k4me3.peak);


  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.wt_resuce.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.wt_resuce.H3k4me3.S1.1k",  file.CFP1.wt_resuce.H3k4me3.bw,  file.CFP1.wt_resuce.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=1000, xlim=c(0,120000), ylim=c(0,120000), file.peaks=file.CFP1.H3k4me3.peak);

}


if(1)
{
  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.wt.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.wt.H3k4me3.S1.0.1k",  file.CFP1.wt.H3k4me3.bw,  file.CFP1.wt.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=100, xlim=c(0,12000), ylim=c(0,12000), file.peaks=file.CFP1.H3k4me3.peak);

  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.ko.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.ko.H3k4me3.S1.0.1k",  file.CFP1.ko.H3k4me3.bw,  file.CFP1.ko.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=100, xlim=c(0,12000), ylim=c(0,12000), file.peaks=file.CFP1.H3k4me3.peak);

  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.C169A.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.C169A.H3k4me3.S1.0.1k",  file.CFP1.C169A.H3k4me3.bw,  file.CFP1.C169A.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=100, xlim=c(0,12000), ylim=c(0,12000), file.peaks=file.CFP1.H3k4me3.peak);


  add_spikes_to_removed_regs( file.CFP1.H3k4me3.peak, file.CFP1.wt_resuce.H3k4me3.bw )
  r0 <- draw_cor_signal("raw.CFP1.wt_resuce.H3k4me3.S1.0.1k",  file.CFP1.wt_resuce.H3k4me3.bw,  file.CFP1.wt_resuce.H3k4me3.pred.raw,    
                  NULL,   chr=NULL, win.size=100, xlim=c(0,12000), ylim=c(0,12000), file.peaks=file.CFP1.H3k4me3.peak);

}