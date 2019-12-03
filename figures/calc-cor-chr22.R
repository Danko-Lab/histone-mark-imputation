source("draw-cor-base.R");
source("draw-cor-param.R");

NCORES <- 10;

df.G1.chr22 <- c();
for(zoom_rate in c(10, 1, 0.1, 0.01)) 
{
   L=mclapply( 1:NROW(df.G1.info), function(i) {
      rx <- draw_cor( paste( df.G1.info[i,"Cell"], df.G1.info[i,"Mark"], df.G1.info[i,"Ver"], sep="." ),   
                         df.G1.info[i,"file.pred.bw"], 
                         df.G1.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.G1.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.G1.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.G1.info[i,"file.unmap.bed"],  
                         file.black.bed = df.G1.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22",
                         chromStart=23830526,
                         chromEnd=48129895, 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES )
    df.G1.chr22 <- rbind(df.G1.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.GM.chr22 <- c();
for(zoom_rate in c(10, 1, 0.1, 0.01))
{
   L=mclapply( 1:NROW(df.GM.info), function(i) {
      rx <- draw_cor( paste( df.GM.info[i,"Cell"], df.GM.info[i,"Mark"], df.GM.info[i,"Ver"], sep="." ),   
                         df.GM.info[i,"file.pred.bw"], 
                         df.GM.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.GM.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.GM.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.GM.info[i,"file.unmap.bed"],  
                         file.black.bed = df.GM.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = TRUE, 
                         overwrite=TRUE,
                         save.rdata=TRUE);}, mc.cores=NCORES  )
    df.GM.chr22 <- rbind(df.GM.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}


df.CD4.chr22 <- c();
for(zoom_rate in c( 10, 1, 0.1, 0.01))
{
   L=mclapply( 1:NROW(df.CD4.info), function(i) {
      rx <- draw_cor( paste( df.CD4.info[i,"Cell"], df.CD4.info[i,"Mark"], df.CD4.info[i,"Ver"], sep="." ),   
                         df.CD4.info[i,"file.pred.bw"], 
                         df.CD4.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.CD4.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.CD4.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.CD4.info[i,"file.unmap.bed"],  
                         file.black.bed = df.CD4.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES  )
    df.CD4.chr22 <- rbind(df.CD4.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}


df.HELA.chr22 <- c();
for(zoom_rate in c(0.01, 0.1, 1, 10 ))
{
   L=mclapply( 1:NROW(df.HELA.info), function(i) {
      rx <- draw_cor( paste( df.HELA.info[i,"Cell"], df.HELA.info[i,"Mark"], df.HELA.info[i,"Ver"], sep="." ),   
                         df.HELA.info[i,"file.pred.bw"], 
                         df.HELA.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.HELA.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.HELA.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.HELA.info[i,"file.unmap.bed"],  
                         file.black.bed = df.HELA.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES  )
    df.HELA.chr22 <- rbind(df.HELA.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
   
}

df.HCT.chr22 <- c();
for(zoom_rate in c(0.01, 0.1, 1, 10 ))
{
   L=mclapply( 1:NROW(df.HCT.info), function(i){
      rx <- draw_cor( paste( df.HCT.info[i,"Cell"], df.HCT.info[i,"Mark"], df.HCT.info[i,"Ver"], sep="." ),   
                         df.HCT.info[i,"file.pred.bw"], 
                         df.HCT.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.HCT.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.HCT.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.HCT.info[i,"file.unmap.bed"],  
                         file.black.bed = df.HCT.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES  )

    df.HCT.chr22 <- rbind(df.HCT.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.Alex.K562.chr22 <- c();
for(zoom_rate in c(0.01, 0.1, 1, 10))
{
   L=mclapply( 1:NROW(df.Alex.K562.info), function(i){
      rx <- draw_cor( paste( df.Alex.K562.info[i,"Cell"], df.Alex.K562.info[i,"Mark"], df.Alex.K562.info[i,"Ver"], sep="." ),   
                         df.Alex.K562.info[i,"file.pred.bw"], 
                         df.Alex.K562.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.Alex.K562.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.Alex.K562.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.Alex.K562.info[i,"file.unmap.bed"],  
                         file.black.bed = df.Alex.K562.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22",
                         chromStart=23830526,
                         chromEnd=48129895, 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES  )
    df.Alex.K562.chr22 <- rbind( df.Alex.K562.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.Alex.Mnase.chr22 <- c();
for(zoom_rate in c(0.01, 0.1, 1, 10))
{
   L=mclapply( 1:NROW(df.Alex.Mnase.info), function(i){
      rx <- draw_cor( paste( df.Alex.Mnase.info[i,"Cell"], df.Alex.Mnase.info[i,"Mark"], df.Alex.Mnase.info[i,"Ver"], sep="." ),   
                         df.Alex.Mnase.info[i,"file.pred.bw"], 
                         df.Alex.Mnase.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.Alex.Mnase.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.Alex.Mnase.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.Alex.Mnase.info[i,"file.unmap.bed"],  
                         file.black.bed = df.Alex.Mnase.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22",
                         chromStart=23830526,
                         chromEnd=48129895, 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES  )

    df.Alex.Mnase.chr22 <- rbind( df.Alex.Mnase.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.Horse.chr22 <- c();
for(zoom_rate in c(0.01, 0.1, 1, 10 ))
{
   L=mclapply( 1:NROW(df.Don.Horse.info), function(i){
      rx <- draw_cor( paste( df.Don.Horse.info[i,"Cell"], df.Don.Horse.info[i,"Mark"], df.Don.Horse.info[i,"Ver"], sep="." ),   
                         df.Don.Horse.info[i,"file.pred.bw"], 
                         df.Don.Horse.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.Don.Horse.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.Don.Horse.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.Don.Horse.info[i,"file.unmap.bed"],  
                         file.black.bed = df.Don.Horse.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr22", 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES  )
    df.Horse.chr22 <- rbind( df.Horse.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}

df.mESC.chr1 <- c();
for(zoom_rate in c(0.01, 0.1, 1, 10 ))
{
   L=mclapply( 1:NROW(df.MM9.mESC.info), function(i){
      rx <- draw_cor( paste( df.MM9.mESC.info[i,"Cell"], df.MM9.mESC.info[i,"Mark"], df.MM9.mESC.info[i,"Ver"], sep="." ),   
                         df.MM9.mESC.info[i,"file.pred.bw"], 
                         df.MM9.mESC.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.MM9.mESC.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.MM9.mESC.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.MM9.mESC.info[i,"file.unmap.bed"],  
                         file.black.bed = df.MM9.mESC.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr1", 
                         gplot = TRUE, 
                         overwrite=TRUE);}, mc.cores=NCORES  )
    df.mESC.chr1 <- rbind( df.mESC.chr1, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
}


save(df.mESC.chr1, df.G1.chr22, df.GM.chr22, df.CD4.chr22, df.HELA.chr22, df.HCT.chr22, df.Alex.K562.chr22, df.Alex.Mnase.chr22, df.Horse.chr22, file="draw-cor-chr22.rdata")



