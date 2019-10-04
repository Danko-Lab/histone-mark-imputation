source("draw-cor-base.R");

combine_res<-function()
{
    tb.mark <- read.csv("histone-res1.csv", header=F); 
    tb.public <- read.table("public-res1.tsv", sep="\t", header=T); 

    library(tools)

    tb.public$FILE <- basename(as.vector(tb.public$URL));
    for(i in 1:NROW(tb.public))
    {
        if (file_ext(tb.public$FILE[i])=="gz")
        {
           str.tmp <- file_path_sans_ext(tb.public$FILE[i])
           tb.public$FILE[i] <- paste(file_path_sans_ext(str.tmp), "bigWig", sep=".")
        }   
    }

    tb.public <- tb.public[, c(2,3,7)]
    tb.public <- tb.public[tb.public$FILE!="",]
    tb.public$FILE <- paste0("/local/workdir/zw355/proj/prj15-histone/VerifyHistone/",tb.public$FILE);
    tb.public$TID <- 0;

    tb.mark <- tb.mark[,-2];
    colnames(tb.mark) <- c("TID", "Cell", "Mark", "FILE");
    tb.mark$FILE <- trimws(tb.mark$FILE);
    return(rbind(tb.mark, tb.public))
}


calculate_matrix<-function(tb.hist, cell, marker, file.hist.ctrl, file.hist.peak, chr="chr22", win.size=10000)
{
    tb.hist0 <- tb.hist[tb.hist$Mark==marker & tb.hist$Cell==cell,,drop=F ]
    if(NROW(tb.hist0)<2) return(NULL);
    
    x <- do.call("rbind", lapply(1:(NROW(tb.hist0$FILE)-1), function(i){return(data.frame(Var1=tb.hist0$FILE[i], Var2=tb.hist0$FILE[(i+1):NROW(tb.hist0$FILE)])) }));

    tb.cor <-  mclapply(1:NROW(x), function(i) {
         rx <- draw_cor("Noname",  
                         file.pred.bw = as.character(x$Var2[i]), 
                         file.exp.bw = as.character(x$Var1[i]), 
                         file.ctrl.bw = file.hist.ctrl, 
                         file.exp.peak = file.hist.peak,  
                         file.unmap.bed = file.hg19.unmap.bed,
                         file.black.bed = file.hg19.black.bed, 
                         chr=chr, 
                         win.size=win.size,
                         gplot=FALSE);
         return(rx);
      }, mc.cores=12, mc.preschedule=FALSE);   
    
    r <- data.frame(Cell=cell, Mark=marker, do.call("rbind",tb.cor)) ;
    rownames(r) <- NULL;
    
    gc(reset=TRUE);
    return(r);
}

get_cor_xkb<-function(win.rate=1)
{
  tb.hist <- combine_res();
  #Do we need to remove training dataset?
  #tb.hist <- tb.hist[tb.hist$TID==0,]

  # G1/chr22
  r1 <- calculate_matrix(tb.hist, "K562", "H3k27ac",  file.k562.H3k27ac.ctrl.bw,  file.k562.H3k27ac.peak,  chr="chr22", win.size=1000*win.rate);
  r2 <- calculate_matrix(tb.hist, "K562", "H3k27me3", file.k562.H3k27me3.ctrl.bw, file.k562.H3k27me3.peak, chr="chr22", win.size=1000*win.rate);
  r3 <- calculate_matrix(tb.hist, "K562", "H3k36me3", file.k562.H3k36me3.ctrl.bw, file.k562.H3k36me3.peak, chr="chr22", win.size=1000*win.rate);
  r4 <- calculate_matrix(tb.hist, "K562", "H3k4me1",  file.k562.H3k4me1.ctrl.bw,  file.k562.H3k4me1.peak,  chr="chr22", win.size=1000*win.rate);
  r5 <- calculate_matrix(tb.hist, "K562", "H3k4me2",  file.k562.H3k4me2.ctrl.bw,  file.k562.H3k4me2.peak,  chr="chr22", win.size=1000*win.rate);
  r6 <- calculate_matrix(tb.hist, "K562", "H3k4me3",  file.k562.H3k4me3.ctrl.bw,  file.k562.H3k4me3.peak,  chr="chr22", win.size=1000*win.rate);
  r7 <- calculate_matrix(tb.hist, "K562", "H3k9ac",   file.k562.H3k9ac.ctrl.bw,   file.k562.H3k9ac.peak,   chr="chr22", win.size=1000*win.rate);
  r8 <- calculate_matrix(tb.hist, "K562", "H3k9me3",  file.k562.H3k9me3.ctrl.bw,  file.k562.H3k9me3.peak,  chr="chr22", win.size=1000*win.rate);
  r9 <- calculate_matrix(tb.hist, "K562", "H4k20me1", file.k562.H4k20me1.ctrl.bw, file.k562.H4k20me1.peak, chr="chr22", win.size=1000*win.rate);

  df.K562.allcomb <- rbind(r1, r2, r3, r4, r5, r6, r7,r8, r9);


  # GM12878/chr22
  r1 <- calculate_matrix(tb.hist, "GM12878", "H3k27ac",  file.gm.H3k27ac.ctrl.bw,  file.gm.H3k27ac.peak,  chr="chr22", win.size=1000*win.rate);
  r2 <- calculate_matrix(tb.hist, "GM12878", "H3k27me3", file.gm.H3k27me3.ctrl.bw, file.gm.H3k27me3.peak, chr="chr22", win.size=1000*win.rate);
  r3 <- calculate_matrix(tb.hist, "GM12878", "H3k36me3", file.gm.H3k36me3.ctrl.bw, file.gm.H3k36me3.peak, chr="chr22", win.size=1000*win.rate);
  r4 <- calculate_matrix(tb.hist, "GM12878", "H3k4me1",  file.gm.H3k4me1.ctrl.bw,  file.gm.H3k4me1.peak,  chr="chr22", win.size=1000*win.rate);
  r5 <- calculate_matrix(tb.hist, "GM12878", "H3k4me2",  file.gm.H3k4me2.ctrl.bw,  file.gm.H3k4me2.peak,  chr="chr22", win.size=1000*win.rate);
  r6 <- calculate_matrix(tb.hist, "GM12878", "H3k4me3",  file.gm.H3k4me3.ctrl.bw,  file.gm.H3k4me3.peak,  chr="chr22", win.size=1000*win.rate);
  r7 <- calculate_matrix(tb.hist, "GM12878", "H3k9ac",   file.gm.H3k9ac.ctrl.bw,   file.gm.H3k9ac.peak,   chr="chr22", win.size=1000*win.rate);
  r8 <- calculate_matrix(tb.hist, "GM12878", "H3k9me3",  file.gm.H3k9me3.ctrl.bw,  file.gm.H3k9me3.peak,  chr="chr22", win.size=1000*win.rate);
  r9 <- calculate_matrix(tb.hist, "GM12878", "H4k20me1", file.gm.H4k20me1.ctrl.bw, file.gm.H4k20me1.peak, chr="chr22", win.size=1000*win.rate);

  df.GM12878.allcomb <- rbind(r1, r2, r3, r4, r5, r6, r7,r8, r9);
  
  return( list(K562=df.K562.allcomb, GM=df.GM12878.allcomb, win.size=1000*win.rate))
}

if(1)
{
    df.cor.10k <- df.cor.1k <- df.cor.100 <- df.cor.10  <-c();
    df.cor.10k <- get_cor_xkb(10);gc();
    df.cor.1k <- get_cor_xkb(1);gc();
    save(df.cor.10k, df.cor.1k, df.cor.100, df.cor.10, file="calc-comb-GM-K562.rdata");
    df.cor.100 <- get_cor_xkb(0.1);gc();
    df.cor.10 <- get_cor_xkb(0.01);gc();
    save(df.cor.10k, df.cor.1k, df.cor.100, df.cor.10, file="calc-comb-GM-K562.rdata");
}



