source("draw-cor-base.R");
source("draw-cor-param.R");

# G1/all chrs/10K
df.G1.all.10k <- c();
for(i in 1:NROW(df.G1.info))
{
   rx <- draw_cor( paste( df.G1.info[i,"Cell"], df.G1.info[i,"Mark"], df.G1.info[i,"Ver"], sep="." ),   
                         df.G1.info[i,"file.pred.bw"], 
                         df.G1.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.G1.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.G1.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.G1.info[i,"file.unmap.bed"],  
                         file.black.bed = df.G1.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*10,
                         chr = NULL,
                         overwrite=TRUE, 
                         gplot = TRUE );
   df.G1.all.10k <- rbind(df.G1.all.10k, rx);
}

# G1/all chrs/10K
df.GM.all.10k <- c();
for(i in 1:NROW(df.GM.info))
{
   rx <- draw_cor( paste( df.GM.info[i,"Cell"], df.GM.info[i,"Mark"], df.GM.info[i,"Ver"], sep="." ),   
                         df.GM.info[i,"file.pred.bw"], 
                         df.GM.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.GM.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.GM.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.GM.info[i,"file.unmap.bed"],  
                         file.black.bed = df.GM.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         win.size = 1000*10,
                         chr = NULL, 
                         overwrite=TRUE,
                         gplot = TRUE );
   df.GM.all.10k <- rbind(df.GM.all.10k, rx);
}


# G1/all peaks/win=1K/ext_size=0k
df.G1.peaks.w1k.ext0k <- c();
for(i in 1:NROW(df.G1.info))
{
   rx <- draw_cor( paste( df.G1.info[i,"Cell"], df.G1.info[i,"Mark"], df.G1.info[i,"Ver"], sep="." ),   
                         df.G1.info[i,"file.pred.bw"], 
                         df.G1.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.G1.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.G1.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.G1.info[i,"file.unmap.bed"],  
                         file.black.bed = df.G1.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         file.draw.peak = df.G1.info[i,"file.exp.peak"],  
                         win.size = 1000,
                         peak.ext=0,
                         chr = NULL, 
                         overwrite=TRUE,
                         gplot = TRUE );
   df.G1.peaks.w1k.ext0k <- rbind(df.G1.peaks.w1k.ext0k, rx);
}


# G1/all peaks/win=1K/ext_size=10k
df.G1.peaks.w1k.ext10k <- c();
for(i in 1:NROW(df.G1.info))
{
   rx <- draw_cor( paste( df.G1.info[i,"Cell"], df.G1.info[i,"Mark"], df.G1.info[i,"Ver"], sep="." ),   
                         df.G1.info[i,"file.pred.bw"], 
                         df.G1.info[i,"file.exp.bw"], 
                         file.ctrl.bw = df.G1.info[i,"file.ctrl.bw"],  
                         file.exp.peak = df.G1.info[i,"file.exp.peak"],  
                         file.unmap.bed = df.G1.info[i,"file.unmap.bed"],  
                         file.black.bed = df.G1.info[i,"file.black.bed"],  
                         file.cell.remove = NULL,
                         file.draw.peak = df.G1.info[i,"file.exp.peak"],  
                         win.size = 1000,
                         peak.ext=10000,
                         chr = NULL, 
                         overwrite=TRUE,
                         gplot = TRUE );
   df.G1.peaks.w1k.ext10k <- rbind(df.G1.peaks.w1k.ext10k, rx);
}


save(df.G1.all.10k, df.GM.all.10k, df.G1.peaks.w1k.ext0k, df.G1.peaks.w1k.ext10k, file="draw-cor-hg19-all.rdata" )
