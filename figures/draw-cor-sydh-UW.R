source("draw-cor-base.R");

g_noplot <<- FALSE

file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )
tb.unmap.bed <- read.table(file.unmap.bed)[,c(1:3)];
tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
rm(tb.unmap.bed);

zoom_rate=1;

file.UW.GM.H3k27me3   <- "../UwHistone/wgEncodeUwHistoneGm12878H3k27me3StdRawRep2.bigWig"
file.UW.GM.H3k36me3   <- "../UwHistone/wgEncodeUwHistoneGm12878H3k36me3StdRawRep2.bigWig"
file.UW.GM.H3k4me3    <- "../UwHistone/wgEncodeUwHistoneGm12878H3k4me3StdRawRep2.bigWig"
file.UW.HELA.H3k27me3 <- "../UwHistone/wgEncodeUwHistoneHelas3H3k27me3StdRawRep2.bigWig"
file.UW.HELA.H3k36me3 <- "../UwHistone/wgEncodeUwHistoneHelas3H3k36me3StdRawRep2.bigWig"
file.UW.HELA.H3k4me3  <- "../UwHistone/wgEncodeUwHistoneHelas3H3k4me3StdRawRep2.bigWig"
file.UW.HCT.H3k4me3   <- "../UwHistone/wgEncodeUwHistoneHct116H3k4me3StdRawRep2.bigWig"
file.UW.K562.H3k4me3  <- "../UwHistone/wgEncodeUwHistoneK562H3k04me3StdZnf2c10c5RawRep2.bigWig"
file.UW.K562.H3k27me3 <- "../UwHistone/wgEncodeUwHistoneK562H3k27me3StdRawRep2.bigWig"
file.UW.K562.H3k36me3 <- "../UwHistone/wgEncodeUwHistoneK562H3k36me3StdRawRep2.bigWig"

file.SYDH.HCT.H3k4me1     <- "../SydhHistone/wgEncodeSydhHistoneHct116H3k04me1UcdSig.bigWig"
file.SYDH.HCT.H3k27ac     <- "../SydhHistone/wgEncodeSydhHistoneHct116H3k27acUcdSig.bigWig"
file.SYDH.K562.H3k27me3   <- "../SydhHistone/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig"
file.SYDH.K562.H3k4me1    <- "../SydhHistone/wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig"
file.SYDH.K562.H3k4me3    <- "../SydhHistone/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig"
file.SYDH.K562.H3k9ac     <- "../SydhHistone/wgEncodeSydhHistoneK562H3k9acbUcdSig.bigWig"


if(1)
{
  # G1/chr22/10K
  add_spikes_to_removed_regs( file.k562.H3k27me3.peak, file.k562.H3k27me3.bw )
  r2 <- draw_cor_signal("UW.G1.chr22.H3k27me3.S1.V3.10k",  file.k562.H3k27me3.bw,  file.UW.K562.H3k27me3,   file.k562.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k36me3.peak, file.k562.H3k36me3.bw )
  r3 <- draw_cor_signal("UW.G1.chr22.H3k36me3.S1.V3.10k",  file.k562.H3k36me3.bw,  file.UW.K562.H3k36me3,   file.k562.H3k36me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me3.peak, file.k562.H3k4me3.bw )
  r6 <- draw_cor_signal("UW.G1.chr22.H3k4me3.S1.V3.10k",   file.k562.H3k4me3.bw,   file.UW.K562.H3k4me3,    file.k562.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9ac.peak, file.k562.H3k9ac.bw )

  df.UW.G1.chr22.10k <- rbind(r2, r3, r6);
}

if(1)
{
  # GM/chr22/10K
  add_spikes_to_removed_regs( file.gm.H3k27me3.peak, file.gm.H3k27me3.bw )
  r2 <- draw_cor_signal("UW.GM.chr22.H3k27me3.S1.V3.10k",  file.gm.H3k27me3.bw,    file.UW.GM.H3k27me3 ,  file.gm.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k36me3.peak, file.gm.H3k36me3.bw )
  r3 <- draw_cor_signal("UW.GM.chr22.H3k36me3.S1.V3.10k",  file.gm.H3k36me3.bw,    file.UW.GM.H3k27me3 ,  file.gm.H3k36me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k4me3.peak, file.gm.H3k4me3.bw )
  r6 <- draw_cor_signal("UW.GM.chr22.H3k4me3.S1.V3.10k",   file.gm.H3k4me3.bw,     file.UW.GM.H3k4me3 ,   file.gm.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );

  df.UW.GM.chr22.10k <- rbind( r2, r3, r6);
}


if(1)
{
  # HELA/chr22/10K
  add_spikes_to_removed_regs( file.hela.H3k27me3.peak, file.hela.H3k27me3.bw )
  r2 <- draw_cor_signal("UW.HELA.chr22.H3k27me3.S1.V3.10k",  file.hela.H3k27me3.bw, file.UW.HELA.H3k27me3  , file.hela.H3k27me3.ctrl.bw, xlim=c(0, 20000), ylim=c(0, 1000), chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.hela.H3k36me3.peak, file.hela.H3k36me3.bw )
  r3 <- draw_cor_signal("UW.HELA.chr22.H3k36me3.S1.V3.10k",  file.hela.H3k36me3.bw, file.UW.HELA.H3k36me3 , file.hela.H3k36me3.ctrl.bw, xlim=c(0, 20000), ylim=c(0, 1000), chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.hela.H3k4me3.peak, file.hela.H3k4me3.bw )
  r6 <- draw_cor_signal("UW.HELA.chr22.H3k4me3.S1.V3.10k",   file.hela.H3k4me3.bw,  file.UW.HELA.H3k4me3 ,  file.hela.H3k4me3.ctrl.bw,  xlim=c(0, 20000), ylim=c(0, 1000),chr="chr22", range=10000*zoom_rate );

  df.UW.HELA.chr22.10k <- rbind( r2, r3, r6);

}

if(1)
{
# HCT/chr22/10K
  add_spikes_to_removed_regs( file.hct.H3k4me3.peak, file.hct.H3k4me3.bw )
  r6 <- draw_cor_signal("UW.HCT.chr22.H3k4me3.S1.V3.10k",   file.hct.H3k4me3.bw,   file.UW.HCT.H3k4me3,  file.hct.H3k4me3.ctrl.bw,   chr="chr22", range=10000*zoom_rate );
  df.UW.HCT.chr22.10k <- rbind( NA, NA, r6);
}


if(1)
{

  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )
  r1 <- draw_cor_signal("SYDH.HCT.chr22.H3k27ac.S1.V3.10k",  file.k562.H3k27ac.bw,  file.SYDH.HCT.H3k27ac,   file.k562.H3k27ac.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me1.peak, file.k562.H3k4me1.bw )
  r4 <- draw_cor_signal("SYDH.HCT.chr22.H3k4me1.S1.V3.10k",   file.k562.H3k4me1.bw,   file.SYDH.HCT.H3k4me1,    file.k562.H3k4me1.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  df.SYDH.HCT.chr22.10k <- rbind(r1, NA, r4, NA, NA);

  # G1/chr22/10K
  add_spikes_to_removed_regs( file.k562.H3k27me3.peak, file.k562.H3k27me3.bw )
  r2 <- draw_cor_signal("SYDH.G1.chr22.H3k27me3.S1.V3.10k",  file.k562.H3k27me3.bw,  file.SYDH.K562.H3k27me3,   file.k562.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me1.peak, file.k562.H3k4me1.bw )
  r4 <- draw_cor_signal("SYDH.G1.chr22.H3k4me1.S1.V3.10k",   file.k562.H3k4me1.bw,   file.SYDH.K562.H3k4me1,    file.k562.H3k4me1.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me3.peak, file.k562.H3k4me3.bw )
  r6 <- draw_cor_signal("SYDH.G1.chr22.H3k4me3.S1.V3.10k",   file.k562.H3k4me3.bw,   file.SYDH.K562.H3k4me3,    file.k562.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9ac.peak, file.k562.H3k9ac.bw )
  r7 <- draw_cor_signal("SYDH.G1.chr22.H3k9ac.S1.V2.10k",    file.k562.H3k9ac.bw,    file.SYDH.K562.H3k9ac,     file.k562.H3k9ac.ctrl.bw,   chr="chr22", range=10000*zoom_rate );

  df.SYDH.G1.chr22.10k <- rbind(NA, r2, r4, r6, r7);
}

if(0)
{
save(df.SYDH.G1.chr22.10k, df.UW.G1.chr22.10k, df.UW.GM.chr22.10k, df.UW.HELA.chr22.10k, df.UW.HCT.chr22.10k , file="df.UW.SYDH.compare.rdata")

#> df.SYDH.G1.chr22.10k                                                                                                                                                                                                                                                        
#                              plot rem.count pearson spearman   mad                                                                                                                                                                                                           
#1 SYDH.G1.chr22.H3k27me3.S1.V3.10k        10  0.9192   0.8574 10957                                                                                                                                                                                                           
#2  SYDH.G1.chr22.H3k4me1.S1.V3.10k        13  0.6283   0.8113 19309                                                                                                                                                                                                           
#3  SYDH.G1.chr22.H3k4me3.S1.V3.10k        18  0.9030   0.5238  4483                                                                                                                                                                                                           
#4   SYDH.G1.chr22.H3k9ac.S1.V2.10k        20  0.8207   0.7219 16157                                                                                                                                                                                                           
#> df.UW.G1.chr22.10k                                                                                                                                                                                                                                                          
#                            plot rem.count pearson spearman  mad                                                                                                                                                                                                              
#1 UW.G1.chr22.H3k27me3.S1.V3.10k        10  0.7548   0.6865 2938                                                                                                                                                                                                              
#2 UW.G1.chr22.H3k36me3.S1.V3.10k        15  0.8503   0.7773 3358                                                                                                                                                                                                              
#3  UW.G1.chr22.H3k4me3.S1.V3.10k        18  0.9355   0.5915 1400                                                                                                                                                                                                              
#> df.UW.GM.chr22.10k                                                                                                                                                                                                                                                          
#                            plot rem.count pearson spearman  mad                                                                                                                                                                                                              
#1 UW.GM.chr22.H3k27me3.S1.V3.10k         1  0.8296   0.7951 3655                                                                                                                                                                                                              
#2 UW.GM.chr22.H3k36me3.S1.V3.10k         0 -0.4373  -0.5812 7061                                                                                                                                                                                                              
#3  UW.GM.chr22.H3k4me3.S1.V3.10k         0  0.9198   0.6023 1353                                                                                                                                                                                                              
#> df.UW.HELA.chr22.10k                                                                                                                                                                                                                                                        
#                              plot rem.count pearson spearman    mad                                                                                                                                                                                                          
#1 UW.HELA.chr22.H3k27me3.S1.V3.10k         6  0.7204   0.6017 1943.0                                                                                                                                                                                                          
#2 UW.HELA.chr22.H3k36me3.S1.V3.10k         8  0.8909   0.7752 2334.6                                                                                                                                                                                                          
#3  UW.HELA.chr22.H3k4me3.S1.V3.10k        20  0.9355   0.5765  808.2                                                                                                                                                                                                          
#> df.UW.HCT.chr22.10k                                                                                                                                                                                                                                                         
#                            plot rem.count pearson spearman  mad                                                                                                                                                                                                              
#1                           <NA>        NA      NA       NA   NA                                                                                                                                                                                                              
#2                           <NA>        NA      NA       NA   NA                                                                                                                                                                                                              
#3 UW.HCT.chr22.H3k4me3.S1.V3.10k         0  0.9563   0.5677 1596           
#

  df.merge <- cbind(G1=df.UW.G1.chr22.10k[,3],  GM=df.UW.GM.chr22.10k[,3], HCT=df.UW.HCT.chr22.10k[,3], HELA=df.UW.HELA.chr22.10k[,3] ) 
  rownames(df.merge) <- c("H3k27me3", "H3k36me3", "H3k4me3" )
  mat_data <- df.merge;
  mat_data <- round(mat_data,2);

  library(gplots)
  library(RColorBrewer)

  # creates a own color palette from red to green
  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.33,length=100),   # for red
  seq(0.34,0.66,length=100),            # for yellow
  seq(0.67,1,length=100))              # for green

# creates a 5 x 5 inch image
  pdf("draw-pearson-UW-10kb.pdf")

 heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Pearson Correlation(Broad vs UW)", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",           # turn off column clustering        
  cexRow=1.5, 
  cexCol=1.5)            

  dev.off()


  df.merge <- cbind(G1=df.SYDH.G1.chr22.10k[,3], HCT=df.SYDH.HCT.chr22.10k[,3]) 
  rownames(df.merge) <- c("H3k27ac", "H3k27me3", "H3k4me1", "H3k4me3", "H3k9ac" )
  mat_data <- df.merge;
  mat_data <- round(mat_data,2);

  library(gplots)
  library(RColorBrewer)

  # creates a own color palette from red to green
  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.33,length=100),   # for red
  seq(0.34,0.66,length=100),            # for yellow
  seq(0.67,1,length=100))              # for green

# creates a 5 x 5 inch image
  pdf("draw-pearson-SYDH-10kb.pdf")

 heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Pearson Correlation(Broad vs SYDH)", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",    # only draw a row dendrogram
  Colv="NA",            # turn off column clustering        
  Rowv="NA",
  cexRow=1.5, 
  cexCol=1.5)            

  dev.off()


}