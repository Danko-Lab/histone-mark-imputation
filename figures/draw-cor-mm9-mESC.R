source("draw-cor-base.R");
source("draw-cor-param.R");

file.mESC.H3k4me3.pred.raw<- "../pred-mm9-mESC/raw.mESC.H3k4me3.trainbyself.S1.bw"
file.mESC.H3k27me3.pred.raw<- "../pred-mm9-mESC/raw.mESC.H3k27me3.trainbyself.S1.bw"

for(zoom_rate in c(1, 10 ))
{
      r1 <- draw_cor( paste( df.MM9.mESC.info[1,"Cell"], df.MM9.mESC.info[1,"Mark"], "trainbyself",df.MM9.mESC.info[1,"Ver"], sep="." ),   
                         #df.MM9.mESC.info[1,"file.pred.bw"], 
                         file.mESC.H3k27me3.pred.raw, 
                         df.MM9.mESC.info[1,"file.exp.bw"], 
                         file.ctrl.bw = df.MM9.mESC.info[1,"file.ctrl.bw"],  
                         file.exp.peak = df.MM9.mESC.info[1,"file.exp.peak"],  
                         file.unmap.bed = df.MM9.mESC.info[1,"file.unmap.bed"],  
                         file.black.bed = df.MM9.mESC.info[1,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr1", 
                         gplot = TRUE, 
                         overwrite=TRUE)

}
for(zoom_rate in c(1, 10 ))
{
    r2 <- draw_cor( paste( df.MM9.mESC.info[2,"Cell"], df.MM9.mESC.info[2,"Mark"], "trainbyself", df.MM9.mESC.info[2,"Ver"], sep="." ),   
                         #df.MM9.mESC.info[2,"file.pred.bw"],   
                         file.mESC.H3k4me3.pred.raw, 
                         df.MM9.mESC.info[2,"file.exp.bw"], 
                         file.ctrl.bw = df.MM9.mESC.info[2,"file.ctrl.bw"],  
                         file.exp.peak = df.MM9.mESC.info[2,"file.exp.peak"],  
                         file.unmap.bed = df.MM9.mESC.info[2,"file.unmap.bed"],  
                         file.black.bed = df.MM9.mESC.info[2,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*zoom_rate,
                         chr = "chr1", 
                         gplot = TRUE, 
                         overwrite=TRUE)
}

