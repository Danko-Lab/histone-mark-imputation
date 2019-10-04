source("draw-cor-base.R");
source("draw-cor-param.R");


if(0)
{
  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )

  r1 <- draw_cor_signal("raw.G1.H3k27ac.S3.V1.peaks.10k",  file.k562.H3k27ac.bw,  "../pred-k562/raw.H3k27ac.S3.V1.G1.bw",    file.k562.H3k27ac.ctrl.bw,   chr="chr22", win.size=10000, file.peaks=file.k562.H3k27ac.peak );
  r2 <- draw_cor_signal("raw.G1.H3k27ac.S3.V1.peaks.1k",   file.k562.H3k27ac.bw,  "../pred-k562/raw.H3k27ac.S3.V1.G1.bw",    file.k562.H3k27ac.ctrl.bw,   chr="chr22", win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k27ac.peak );
  r3 <- draw_cor_signal("raw.G1.H3k27ac.S3.V1.peaks.500",  file.k562.H3k27ac.bw,  "../pred-k562/raw.H3k27ac.S3.V1.G1.bw",    file.k562.H3k27ac.ctrl.bw,   chr="chr22", win.size=500, xlim=c(0,7500), ylim=c(0,7500),file.peaks=file.k562.H3k27ac.peak );
  r4 <- draw_cor_signal("raw.G1.H3k27ac.S3.V1.peaks.100",  file.k562.H3k27ac.bw,  "../pred-k562/raw.H3k27ac.S3.V1.G1.bw",    file.k562.H3k27ac.ctrl.bw,   chr="chr22", win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.k562.H3k27ac.peak );

  r1 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks.10k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=10000, file.peaks=file.k562.H3k27ac.peak );
  r2 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks.1k",   file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k27ac.peak );
  r3 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks.500",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=500, xlim=c(0,7500), ylim=c(0,7500),file.peaks=file.k562.H3k27ac.peak );
  r4 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks.100",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.k562.H3k27ac.peak );
}

if(0)
{
  add_spikes_to_removed_regs( file.gm.H3k27ac.peak, file.gm.H3k27ac.bw )
  r1 <- draw_cor_signal("GM.all.H3k27ac.S1.V3.peaks.10k",   file.gm.H3k27ac.bw,   file.gm.H3k27ac.pred.bw,     file.gm.H3k27ac.ctrl.bw,  chr=NULL, win.size=10000, file.peaks=file.gm.H3k27ac.peak  );
  r1 <- draw_cor_signal("GM.all.H3k27ac.S1.V3.peaks.100",   file.gm.H3k27ac.bw,   file.gm.H3k27ac.pred.bw,     file.gm.H3k27ac.ctrl.bw,  chr=NULL, win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.gm.H3k27ac.peak);

  r1 <- draw_cor_signal("raw.GM.H3k27ac.S1.V3.peaks.10k",   file.gm.H3k27ac.bw,   file.gm.H3k27ac.pred.raw,    file.gm.H3k27ac.ctrl.bw,  chr="chr22", win.size=10000, file.peaks=file.gm.H3k27ac.peak );
  r2 <- draw_cor_signal("raw.GM.H3k27ac.S1.V3.peaks.100",   file.gm.H3k27ac.bw,   file.gm.H3k27ac.pred.raw,    file.gm.H3k27ac.ctrl.bw,  chr="chr22", win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.gm.H3k27ac.peak );
  r3 <- draw_cor_signal("raw.GM.H3k27ac.S1.V3.peaks.500",   file.gm.H3k27ac.bw,   file.gm.H3k27ac.pred.raw,    file.gm.H3k27ac.ctrl.bw,  chr="chr22", win.size=500, xlim=c(0,7500), ylim=c(0,7500),file.peaks=file.gm.H3k27ac.peak );
  r4 <- draw_cor_signal("raw.GM.H3k27ac.S1.V3.peaks.1k",    file.gm.H3k27ac.bw,   file.gm.H3k27ac.pred.raw,    file.gm.H3k27ac.ctrl.bw,  chr="chr22", win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.gm.H3k27ac.peak );
}



if(0)
{
  r1 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext100.10k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=10000, file.peaks=file.k562.H3k27ac.peak, peak.mode=1 );
  r2 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext100.1k",   file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k27ac.peak, peak.mode=1 );
  r3 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext100.500",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=500, xlim=c(0,7500), ylim=c(0,7500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1 );
  r4 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext100.100",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1 );

  r1 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext500.10k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=10000, file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=500 );
  r2 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext500.1k",   file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=500 );
  r3 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext500.500",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=500, xlim=c(0,7500), ylim=c(0,7500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=500 );
  r4 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext500.100",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=500 );

  r1 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext1k.10k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=10000, file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=1000 );
  r2 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext1k.1k",   file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=1000 );
  r3 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext1k.500",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=500, xlim=c(0,7500), ylim=c(0,7500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=1000 );
  r4 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext1k.100",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=1000 );

  r1 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext10k.10k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=10000, file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=10000 );
  r2 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext10k.1k",   file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=10000 );
  r3 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext10k.500",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=500, xlim=c(0,7500), ylim=c(0,7500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=10000 );
  r4 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext10k.100",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=100, xlim=c(0,1500), ylim=c(0,1500),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=10000 );

}


if(0)
{
  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )
  r2 <- draw_cor_signal("raw.G1.H3k27ac.S1.V3.peaks_ext10k.1k",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.pred.raw,    file.k562.H3k27ac.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k27ac.peak, peak.mode=1, peak.ext=10000 );

  add_spikes_to_removed_regs( file.k562.H3k4me3.peak, file.k562.H3k4me3.bw )
  r2 <- draw_cor_signal("raw.G1.H3k4me3.S1.V3.peaks_ext10k.1k",  file.k562.H3k4me3.bw,  file.k562.H3k4me3.pred.raw,    file.k562.H3k4me3.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k4me3.peak, peak.mode=1, peak.ext=10000 );

  add_spikes_to_removed_regs( file.k562.H3k4me1.peak, file.k562.H3k4me1.bw )
  r2 <- draw_cor_signal("raw.G1.H3k4me1.S1.V2.peaks_ext10k.1k",  file.k562.H3k4me1.bw,  file.k562.H3k4me1.pred.raw,    file.k562.H3k4me1.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k4me1.peak, peak.mode=1, peak.ext=10000 );

  add_spikes_to_removed_regs( file.k562.H3k4me2.peak, file.k562.H3k4me2.bw )
  r2 <- draw_cor_signal("raw.G1.H3k4me2.S1.V2.peaks_ext10k.1k",  file.k562.H3k4me2.bw,  file.k562.H3k4me2.pred.raw,    file.k562.H3k4me2.ctrl.bw,   chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k4me2.peak, peak.mode=1, peak.ext=10000 );

  add_spikes_to_removed_regs( file.k562.H3k9ac.peak, file.k562.H3k9ac.bw )
  r2 <- draw_cor_signal("raw.G1.H3k9ac.S1.V2.peaks_ext10k.1k",   file.k562.H3k9ac.bw,   file.k562.H3k9ac.pred.raw,     file.k562.H3k9ac.ctrl.bw,    chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k9ac.peak, peak.mode=1, peak.ext=10000 );

  add_spikes_to_removed_regs( file.k562.H3k122ac.peak, file.k562.H3k122ac.bw )
  r2 <- draw_cor_signal("raw.G1.H3k122ac.S1.V3.peaks_ext10k.1k", file.k562.H3k122ac.bw, file.k562.H3k122ac.pred.raw,   file.k562.H3k122ac.ctrl.bw,  chr=NULL, win.size=1000, xlim=c(0,15000), ylim=c(0,15000),file.peaks=file.k562.H3k122ac.peak, peak.mode=1, peak.ext=10000 );

  #distributed
  add_spikes_to_removed_regs( file.k562.H3k36me3.peak, file.k562.H3k36me3.bw )
  r2 <- draw_cor_signal("raw.G1.H3k36me3.S1.V3.10k",   file.k562.H3k36me3.bw,  file.k562.H3k36me3.pred.raw,    file.k562.H3k36me3.ctrl.bw,   chr=NULL, win.size=10000, xlim=c(0,15000), ylim=c(0,15000)  );

  add_spikes_to_removed_regs( file.k562.H3k27me3.peak, file.k562.H3k27me3.bw )
  r2 <- draw_cor_signal("raw.G1.H3k27me3.S1.V3.10k",   file.k562.H3k27me3.bw,  file.k562.H3k27me3.pred.raw,    file.k562.H3k27me3.ctrl.bw,   chr=NULL, win.size=10000, xlim=c(0,15000), ylim=c(0,15000)  );

  add_spikes_to_removed_regs( file.k562.H4k20me1.peak, file.k562.H4k20me1.bw )
  r2 <- draw_cor_signal("raw.G1.H4k20me1.S1.V3.10k",   file.k562.H4k20me1.bw,  file.k562.H4k20me1.pred.raw,    file.k562.H4k20me1.ctrl.bw,   chr=NULL, win.size=10000, xlim=c(0,15000), ylim=c(0,15000)  );

  add_spikes_to_removed_regs( file.k562.H3k9ac.peak,   file.k562.H3k9ac.bw )
  r2 <- draw_cor_signal("raw.G1.H3k9ac.S1.V2.10k",     file.k562.H3k9ac.bw,    file.k562.H3k9ac.pred.raw,      file.k562.H4k20me1.ctrl.bw,   chr=NULL, win.size=10000, xlim=c(0,15000), ylim=c(0,15000)  );
}


if(0)
{
  df.merge <- cbind(G1=df.G1.chr22[,3],  GM=df.GM.chr22[,3], HCT=df.HCT.chr22[,3], HELA=df.HELA.chr22[,3], CD4=df.CD4.chr22[,3], horse=df.horse.chr22[,3]) 
  save(df.merge, df.G1.chr22,  GM=df.GM.chr22, HCT=df.HCT.chr22, HELA=df.HELA.chr22, CD4=df.CD4.chr22, horse=df.horse.chr22, file="df.merge.10k.rdata");
  rownames(df.merge) <- c("H3k122ac", "H3k27ac", "H3k27me3", "H3k36me3", "H3k4me1", "H3k4me2", "H3k4me3", "H3k9ac", "H3k9me3", "H4k20me1")
  #df.merge[is.na(df.merge)] = 0
  df.merge <- t(df.merge^2);
  pdf("R2-cells.pdf")
  nba_heatmap <- heatmap(df.merge, Rowv=NA, Colv=NA, na.rm=TRUE, col =  rev(heat.colors(256)), scale="none", margins=c(5,10))
  dev.off();
}

if(0)
{
  df.merge <- cbind(G1=df.G1.chr22.100k[,3],  GM=df.GM.chr22.100k[,3], HCT=df.HCT.chr22.100k[,3], HELA=df.HELA.chr22.100k[,3], CD4=df.CD4.chr22.100k[,3], horse=df.horse.chr22.100k[,3]) 
  save(df.merge, df.G1.chr22.100k,  GM=df.GM.chr22.100k, HCT=df.HCT.chr22.100k, HELA=df.HELA.chr22.100k, CD4=df.CD4.chr22.100k, horse=df.horse.chr22.100k, file="df.merge.100k.rdata");
  rownames(df.merge) <- c("H3k122ac", "H3k27ac", "H3k27me3", "H3k36me3", "H3k4me1", "H3k4me2", "H3k4me3", "H3k9ac", "H3k9me3", "H4k20me1")
  df.merge <- df.merge[c(3,4,10),]
  df.merge <- t(df.merge^2);
  pdf("R2-100K-cells.pdf")
  nba_heatmap <- heatmap(df.merge, Rowv=NA, Colv=NA, na.rm=TRUE, col =  rev(heat.colors(256)), scale="none", margins=c(5,10))
  dev.off();
}
