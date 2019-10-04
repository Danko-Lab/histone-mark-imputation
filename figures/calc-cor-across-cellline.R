source("draw-cor-base.R");
source("draw-cor-param.R");

NCORES <- 10;
win.levels <- c(10, 1)

find_same_mark_G1 <-function(mark)
{
  ms <- c("H3k122ac", "H3k27ac", "H3k27me3", "H3k36me3", "H3k4me1", "H3k4me2", "H3k4me3", "H3k9ac", "H3k9me3", "H4k20me1")
  
  fs <- c(file.k562.H3k122ac.bw,
      file.k562.H3k27ac.bw ,
      file.k562.H3k27me3.bw,
      file.k562.H3k36me3.bw,
      file.k562.H3k4me1.bw ,
      file.k562.H3k4me2.bw ,
      file.k562.H3k4me3.bw ,
      file.k562.H3k9ac.bw  ,
      file.k562.H3k9me3.bw ,
      file.k562.H4k20me1.bw)

   return(fs[which(ms==mark)]);
}

df.G1GM.chr22 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.GM.info), function(i) {
      rx <- draw_cor( paste( df.GM.info[i,"Cell"], df.GM.info[i,"Mark"], df.GM.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.GM.info[i,"Mark"] ), 
                         df.GM.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.GM.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.GM.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.GM.info[i,"file.unmap.bed"],  
                         file.black.bed = df.GM.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = FALSE, 
                         overwrite=FALSE,
                         save.rdata=FALSE);}, mc.cores=NCORES  )
    df.G1GM.chr22 <- rbind(df.G1GM.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}


df.G1CD4.chr22 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.CD4.info), function(i) {
      rx <- draw_cor( paste( df.CD4.info[i,"Cell"], df.CD4.info[i,"Mark"], df.CD4.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.CD4.info[i,"Mark"] ), 
                         df.CD4.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.CD4.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.CD4.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.CD4.info[i,"file.unmap.bed"],  
                         file.black.bed = df.CD4.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = FALSE, 
                         overwrite=FALSE);}, mc.cores=NCORES  )
    df.G1CD4.chr22 <- rbind(df.G1CD4.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}


df.G1HELA.chr22 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.HELA.info), function(i) {
      rx <- draw_cor( paste( df.HELA.info[i,"Cell"], df.HELA.info[i,"Mark"], df.HELA.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.HELA.info[i,"Mark"] ), 
                         df.HELA.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.HELA.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.HELA.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.HELA.info[i,"file.unmap.bed"],  
                         file.black.bed = df.HELA.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = FALSE, 
                         overwrite=FALSE);}, mc.cores=NCORES  )
    df.G1HELA.chr22 <- rbind(df.G1HELA.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
   
}

df.G1HCT.chr22 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.HCT.info), function(i){
      rx <- draw_cor( paste( df.HCT.info[i,"Cell"], df.HCT.info[i,"Mark"], df.HCT.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.HCT.info[i,"Mark"] ), 
                         df.HCT.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.HCT.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.HCT.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.HCT.info[i,"file.unmap.bed"],  
                         file.black.bed = df.HCT.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = FALSE, 
                         overwrite=FALSE);}, mc.cores=NCORES  )

    df.G1HCT.chr22 <- rbind(df.G1HCT.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.G1Alex.K562.chr22 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.Alex.K562.info), function(i){
      rx <- draw_cor( paste( df.Alex.K562.info[i,"Cell"], df.Alex.K562.info[i,"Mark"], df.Alex.K562.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.Alex.K562.info[i,"Mark"] ), 
                         df.Alex.K562.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.Alex.K562.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.Alex.K562.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.Alex.K562.info[i,"file.unmap.bed"],  
                         file.black.bed = df.Alex.K562.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = FALSE, 
                         overwrite=FALSE);}, mc.cores=NCORES  )
    df.G1Alex.K562.chr22 <- rbind( df.G1Alex.K562.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.G1Alex.Mnase.chr22 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.Alex.Mnase.info), function(i){
      rx <- draw_cor( paste( df.Alex.Mnase.info[i,"Cell"], df.Alex.Mnase.info[i,"Mark"], df.Alex.Mnase.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.Alex.Mnase.info[i,"Mark"] ), 
                         df.Alex.Mnase.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.Alex.Mnase.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.Alex.Mnase.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.Alex.Mnase.info[i,"file.unmap.bed"],  
                         file.black.bed = df.Alex.Mnase.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = FALSE, 
                         overwrite=FALSE);}, mc.cores=NCORES  )

    df.G1Alex.Mnase.chr22 <- rbind( df.G1Alex.Mnase.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.G1Horse.chr22 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.Don.Horse.info), function(i){
      rx <- draw_cor( paste( df.Don.Horse.info[i,"Cell"], df.Don.Horse.info[i,"Mark"], df.Don.Horse.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.Don.Horse.info[i,"Mark"] ), 
                         df.Don.Horse.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.Don.Horse.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.Don.Horse.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.Don.Horse.info[i,"file.unmap.bed"],  
                         file.black.bed = df.Don.Horse.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = FALSE, 
                         overwrite=FALSE);}, mc.cores=NCORES  )
    df.G1Horse.chr22 <- rbind( df.G1Horse.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.G1mESC.chr1 <- c();
for(zoom_rate in win.levels)
{
   L=mclapply( 1:NROW(df.MM9.mESC.info), function(i){
      rx <- draw_cor( paste( df.MM9.mESC.info[i,"Cell"], df.MM9.mESC.info[i,"Mark"], df.MM9.mESC.info[i,"Ver"], sep="." ),   
                         find_same_mark_G1( df.MM9.mESC.info[i,"Mark"] ), 
                         df.MM9.mESC.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.MM9.mESC.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.MM9.mESC.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.MM9.mESC.info[i,"file.unmap.bed"],  
                         file.black.bed = df.MM9.mESC.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr1", 
                         gplot = FALSE, 
                         overwrite=FALSE);}, mc.cores=NCORES  )
    df.G1mESC.chr1 <- rbind( df.G1mESC.chr1, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}


save(df.G1mESC.chr1, 
     df.G1GM.chr22, 
     df.G1CD4.chr22, 
     df.G1HELA.chr22, 
     df.G1HCT.chr22, 
     df.G1Alex.K562.chr22, 
     df.G1Alex.Mnase.chr22, 
     df.G1Horse.chr22, file="draw-cor-across-cellline.rdata")



