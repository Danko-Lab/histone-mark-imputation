##
##
##
if(0)
{
  tb <- read.table("public-res1.tsv", sep="\t", header=TRUE)
  for(i in 1:NROW(tb))
     system(paste0("wget -P /workdir/zw355/proj/prj15-histone/VerifyHistone ", tb[i, "URL"]));
}


source("draw-base.R");

g_noplot <<- TRUE

file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )
tb.unmap.bed <- read.table(file.unmap.bed)[,c(1:3)];
tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
rm(tb.unmap.bed);

zoom_rate=1;

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

tb.hist <- combine_res();
tb.hist <- tb.hist[tb.hist$TID==0,]

calculate_cor<-function(tb.hist, cell, marker, file.hist.broad, file.hist.ctrl, chr="chr22", range=10000)
{
   g_noplot <<- TRUE;

    tb.hist0 <- tb.hist[tb.hist$Mark==marker & tb.hist$Cell==cell,,drop=F ]
    if(NROW(tb.hist0)==0) return(NULL);
    
    tb.cor <-  mclapply(1:NROW(tb.hist0), function(i) {
         rx <- draw_cor_signal("Noname",  file.hist.broad, as.character(tb.hist0$FILE[i]), file.hist.ctrl, chr=chr, range=range );
         return(rx);
      }, mc.cores=10);   
    
    return(cbind(tb.hist0,  do.call("rbind",tb.cor)) );
}  


if(0)
{
  # G1/chr22/10K
  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )
  r1 <- calculate_cor(tb.hist, "K562", "H3k27ac",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k27me3.peak, file.k562.H3k27me3.bw )
  r2 <- calculate_cor(tb.hist, "K562", "H3k27me3", file.k562.H3k27me3.bw, file.k562.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k36me3.peak, file.k562.H3k36me3.bw )
  r3 <- calculate_cor(tb.hist, "K562", "H3k36me3", file.k562.H3k36me3.bw, file.k562.H3k36me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me1.peak, file.k562.H3k4me1.bw )
  r4 <- calculate_cor(tb.hist, "K562", "H3k4me1",  file.k562.H3k4me1.bw,  file.k562.H3k4me1.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me2.peak, file.k562.H3k4me2.bw )
  r5 <- calculate_cor(tb.hist, "K562", "H3k4me2",  file.k562.H3k4me2.bw,  file.k562.H3k4me2.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me3.peak, file.k562.H3k4me3.bw )
  r6 <- calculate_cor(tb.hist, "K562", "H3k4me3",  file.k562.H3k4me3.bw,  file.k562.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9ac.peak, file.k562.H3k9ac.bw )
  r7 <- calculate_cor(tb.hist, "K562", "H3k9ac",   file.k562.H3k9ac.bw,   file.k562.H3k9ac.ctrl.bw,   chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9me3.peak, file.k562.H3k9me3.bw )
  r8 <- calculate_cor(tb.hist, "K562", "H3k9me3",  file.k562.H3k9me3.bw,  file.k562.H3k9me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H4k20me1.peak, file.k562.H4k20me1.bw )
  r9 <- calculate_cor(tb.hist, "K562", "H4k20me1", file.k562.H4k20me1.bw, file.k562.H4k20me1.ctrl.bw, chr="chr22", range=10000*zoom_rate );

  df.K562.chr22 <- rbind(r1, r2, r3, r4, r5, r6, r7,r8, r9);
  
}


if(0)
{
  # GM/chr22/10K
  add_spikes_to_removed_regs( file.gm.H3k27ac.peak, file.gm.H3k27ac.bw )
  r1 <- calculate_cor(tb.hist, "GM12878", "H3k27ac",  file.gm.H3k27ac.bw,  file.gm.H3k27ac.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k27me3.peak, file.gm.H3k27me3.bw )
  r2 <- calculate_cor(tb.hist, "GM12878", "H3k27me3", file.gm.H3k27me3.bw, file.gm.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k36me3.peak, file.gm.H3k36me3.bw )
  r3 <- calculate_cor(tb.hist, "GM12878", "H3k36me3", file.gm.H3k36me3.bw, file.gm.H3k36me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k4me1.peak, file.gm.H3k4me1.bw )
  r4 <- calculate_cor(tb.hist, "GM12878", "H3k4me1",  file.gm.H3k4me1.bw,  file.gm.H3k4me1.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k4me2.peak, file.gm.H3k4me2.bw )
  r5 <- calculate_cor(tb.hist, "GM12878", "H3k4me2",  file.gm.H3k4me2.bw,  file.gm.H3k4me2.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k4me3.peak, file.gm.H3k4me3.bw )
  r6 <- calculate_cor(tb.hist, "GM12878", "H3k4me3",  file.gm.H3k4me3.bw,  file.gm.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k9ac.peak, file.gm.H3k9ac.bw )
  r7 <- calculate_cor(tb.hist, "GM12878", "H3k9ac",   file.gm.H3k9ac.bw,   file.gm.H3k9ac.ctrl.bw,   chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k9me3.peak, file.gm.H3k9me3.bw )
  r8 <- calculate_cor(tb.hist, "GM12878", "H3k9me3",  file.gm.H3k9me3.bw,  file.gm.H3k9me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H4k20me1.peak, file.gm.H4k20me1.bw )
  r9 <- calculate_cor(tb.hist, "GM12878", "H4k20me1", file.gm.H4k20me1.bw, file.gm.H4k20me1.ctrl.bw, chr="chr22", range=10000*zoom_rate );

  df.GM12878.chr22 <- rbind(r1, r2, r3, r4, r5, r6, r7,r8, r9);
}

if(0)
{
	pdf("boxplot-histone-chr22.pdf");

	df.K562.chr22$Mark <- as.character(df.K562.chr22$Mark)
	boxplot(pearson~Mark,data=df.K562.chr22, main="Chr. 22 of K562",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 
	df.GM12878.chr22$Mark <- as.character(df.GM12878.chr22$Mark)
	boxplot(pearson~Mark,data=df.GM12878.chr22, main="Chr. 22 of GM12878",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 
	df.chr22 <- rbind( df.K562.chr22, df.GM12878.chr22);
	boxplot(pearson~Mark,data=df.chr22, main="Chr. 22 of GM12878 + K562",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 

	dev.off();
}


if(0)
{
  # G1/chr22/10K
  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )
  r1 <- calculate_cor(tb.hist, "K562", "H3k27ac",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k27me3.peak, file.k562.H3k27me3.bw )
  r2 <- calculate_cor(tb.hist, "K562", "H3k27me3", file.k562.H3k27me3.bw, file.k562.H3k27me3.ctrl.bw, chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k36me3.peak, file.k562.H3k36me3.bw )
  r3 <- calculate_cor(tb.hist, "K562", "H3k36me3", file.k562.H3k36me3.bw, file.k562.H3k36me3.ctrl.bw, chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me1.peak, file.k562.H3k4me1.bw )
  r4 <- calculate_cor(tb.hist, "K562", "H3k4me1",  file.k562.H3k4me1.bw,  file.k562.H3k4me1.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me2.peak, file.k562.H3k4me2.bw )
  r5 <- calculate_cor(tb.hist, "K562", "H3k4me2",  file.k562.H3k4me2.bw,  file.k562.H3k4me2.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me3.peak, file.k562.H3k4me3.bw )
  r6 <- calculate_cor(tb.hist, "K562", "H3k4me3",  file.k562.H3k4me3.bw,  file.k562.H3k4me3.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9ac.peak, file.k562.H3k9ac.bw )
  r7 <- calculate_cor(tb.hist, "K562", "H3k9ac",   file.k562.H3k9ac.bw,   file.k562.H3k9ac.ctrl.bw,   chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9me3.peak, file.k562.H3k9me3.bw )
  r8 <- calculate_cor(tb.hist, "K562", "H3k9me3",  file.k562.H3k9me3.bw,  file.k562.H3k9me3.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H4k20me1.peak, file.k562.H4k20me1.bw )
  r9 <- calculate_cor(tb.hist, "K562", "H4k20me1", file.k562.H4k20me1.bw, file.k562.H4k20me1.ctrl.bw, chr=NULL, range=10000*zoom_rate );

  df.K562 <- rbind(r1, r2, r3, r4, r5, r6, r7,r8, r9);
  
}


if(0)
{
  # GM/chr22/10K
  add_spikes_to_removed_regs( file.gm.H3k27ac.peak, file.gm.H3k27ac.bw )
  r1 <- calculate_cor(tb.hist, "GM12878", "H3k27ac",  file.gm.H3k27ac.bw,  file.gm.H3k27ac.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k27me3.peak, file.gm.H3k27me3.bw )
  r2 <- calculate_cor(tb.hist, "GM12878", "H3k27me3", file.gm.H3k27me3.bw, file.gm.H3k27me3.ctrl.bw, chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k36me3.peak, file.gm.H3k36me3.bw )
  r3 <- calculate_cor(tb.hist, "GM12878", "H3k36me3", file.gm.H3k36me3.bw, file.gm.H3k36me3.ctrl.bw, chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k4me1.peak, file.gm.H3k4me1.bw )
  r4 <- calculate_cor(tb.hist, "GM12878", "H3k4me1",  file.gm.H3k4me1.bw,  file.gm.H3k4me1.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k4me2.peak, file.gm.H3k4me2.bw )
  r5 <- calculate_cor(tb.hist, "GM12878", "H3k4me2",  file.gm.H3k4me2.bw,  file.gm.H3k4me2.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k4me3.peak, file.gm.H3k4me3.bw )
  r6 <- calculate_cor(tb.hist, "GM12878", "H3k4me3",  file.gm.H3k4me3.bw,  file.gm.H3k4me3.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k9ac.peak, file.gm.H3k9ac.bw )
  r7 <- calculate_cor(tb.hist, "GM12878", "H3k9ac",   file.gm.H3k9ac.bw,   file.gm.H3k9ac.ctrl.bw,   chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H3k9me3.peak, file.gm.H3k9me3.bw )
  r8 <- calculate_cor(tb.hist, "GM12878", "H3k9me3",  file.gm.H3k9me3.bw,  file.gm.H3k9me3.ctrl.bw,  chr=NULL, range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.gm.H4k20me1.peak, file.gm.H4k20me1.bw )
  r9 <- calculate_cor(tb.hist, "GM12878", "H4k20me1", file.gm.H4k20me1.bw, file.gm.H4k20me1.ctrl.bw, chr=NULL, range=10000*zoom_rate );

  df.GM12878 <- rbind(r1, r2, r3, r4, r5, r6, r7,r8, r9);
}

if(0)
{
   pdf("boxplot-histone.pdf");

   df.K562$Mark <- as.character(df.K562$Mark)
   boxplot(pearson~Mark,data=df.K562, main="K562",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 
   df.GM12878$Mark <- as.character(df.GM12878$Mark)
   boxplot(pearson~Mark,data=df.GM12878, main="GM12878",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 
   df.all <- rbind( df.K562, df.GM12878);
   boxplot(pearson~Mark,data=df.all, main="GM12878 + K562",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 

   dev.off();
}


if(0)
{
  # GM/chr22/10K
  add_spikes_to_removed_regs( file.gm.H3k27ac.peak, file.gm.H3k27ac.bw )
  r10 <- calculate_cor(tb.hist, "K562", "H3k27ac",  file.gm.H3k27ac.bw,  file.gm.H3k27ac.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )
  r11 <- calculate_cor(tb.hist, "GM12878", "H3k27ac",  file.k562.H3k27ac.bw,  file.k562.H3k27ac.ctrl.bw,  chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H3k27me3.peak, file.gm.H3k27me3.bw )
  r20 <- calculate_cor(tb.hist, "K562", "H3k27me3", file.gm.H3k27me3.bw, file.gm.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k27me3.peak, file.k562.H3k27me3.bw )
  r21 <- calculate_cor(tb.hist, "GM12878", "H3k27me3", file.k562.H3k27me3.bw, file.k562.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H3k36me3.peak, file.gm.H3k36me3.bw )
  r30 <- calculate_cor(tb.hist, "K562", "H3k36me3", file.gm.H3k36me3.bw, file.gm.H3k36me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k36me3.peak, file.k562.H3k36me3.bw )
  r31 <- calculate_cor(tb.hist, "GM12878", "H3k36me3", file.k562.H3k36me3.bw, file.k562.H3k36me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H3k4me1.peak, file.gm.H3k4me1.bw )
  r40 <- calculate_cor(tb.hist, "K562", "H3k4me1",  file.gm.H3k4me1.bw,  file.gm.H3k4me1.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me1.peak, file.k562.H3k4me1.bw )
  r41 <- calculate_cor(tb.hist, "GM12878", "H3k4me1",  file.k562.H3k4me1.bw,  file.k562.H3k4me1.ctrl.bw,  chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H3k4me2.peak, file.gm.H3k4me2.bw )
  r50 <- calculate_cor(tb.hist, "K562", "H3k4me2",  file.gm.H3k4me2.bw,  file.gm.H3k4me2.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me2.peak, file.k562.H3k4me2.bw )
  r51 <- calculate_cor(tb.hist, "GM12878", "H3k4me2",  file.k562.H3k4me2.bw,  file.k562.H3k4me2.ctrl.bw,  chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H3k4me3.peak, file.gm.H3k4me3.bw )
  r60 <- calculate_cor(tb.hist, "K562", "H3k4me3",  file.gm.H3k4me3.bw,  file.gm.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me3.peak, file.k562.H3k4me3.bw )
  r61 <- calculate_cor(tb.hist, "GM12878", "H3k4me3",  file.k562.H3k4me3.bw,  file.k562.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H3k9ac.peak, file.gm.H3k9ac.bw )
  r70 <- calculate_cor(tb.hist, "K562", "H3k9ac",   file.gm.H3k9ac.bw,   file.gm.H3k9ac.ctrl.bw,   chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9ac.peak, file.k562.H3k9ac.bw )
  r71 <- calculate_cor(tb.hist, "GM12878", "H3k9ac",   file.k562.H3k9ac.bw,   file.k562.H3k9ac.ctrl.bw,   chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H3k9me3.peak, file.gm.H3k9me3.bw )
  r80 <- calculate_cor(tb.hist, "K562", "H3k9me3",  file.gm.H3k9me3.bw,  file.gm.H3k9me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9me3.peak, file.k562.H3k9me3.bw )
  r81 <- calculate_cor(tb.hist, "GM12878", "H3k9me3",  file.k562.H3k9me3.bw,  file.k562.H3k9me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );

  add_spikes_to_removed_regs( file.gm.H4k20me1.peak, file.gm.H4k20me1.bw )
  r90 <- calculate_cor(tb.hist, "K562", "H4k20me1", file.gm.H4k20me1.bw, file.gm.H4k20me1.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H4k20me1.peak, file.k562.H4k20me1.bw )
  r91 <- calculate_cor(tb.hist, "GM12878", "H4k20me1", file.k562.H4k20me1.bw, file.k562.H4k20me1.ctrl.bw, chr="chr22", range=10000*zoom_rate );

  df.mixed.chr22 <- rbind(r10, r11, r20, r21, r30, r31, r40, r41, r50, r51, r60, r61, r70, r71, r80, r81, r90, r91);
  df.mixed.chr22$Mark <- as.character(df.mixed.chr22$Mark)

   pdf("boxplot-histone-mixed.pdf");
   boxplot(pearson~Mark,data=df.mixed.chr22, main="GM12878 + K562",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 
   dev.off();

}  
  
#save(df.all, df.mixed.chr22, file="draw-cor-markers.rdata")


calculate_K562_GM_matrix<-function(tb.hist, marker, file.hist.ctrl, chr="chr22", range=10000)
{
   g_noplot <<- TRUE;

    tb.K562 <- tb.hist[tb.hist$Mark==marker & tb.hist$Cell=="K562",,drop=F ]
    tb.GM   <- tb.hist[tb.hist$Mark==marker & tb.hist$Cell=="GM12878",,drop=F ]
    tb.hist0 <- rbind(tb.K562, tb.GM);
    if(NROW(tb.hist0)<2) return(NULL);
    
    x <- do.call("rbind", lapply(1:(NROW(tb.hist0$FILE)-1), function(i){return(data.frame(Var1=tb.hist0$FILE[i], Var2=tb.hist0$FILE[(i+1):NROW(tb.hist0$FILE)])) }));

    tb.cor <-  mclapply(1:NROW(x), function(i) {
         rx <- draw_cor_signal("Noname",  as.character(x$Var1[i]), as.character(x$Var2[i]), file.hist.ctrl, chr=chr, range=range );
         return(rx);
      }, mc.cores=10);   
    
    r <- data.frame(Cell="K562/GM", Mark=marker, do.call("rbind",tb.cor)) ;
    rownames(r) <- NULL;
    return(r);
}

if(0)
{
  tb.hist <- combine_res();

  # G1/chr22/10K
  add_spikes_to_removed_regs( file.k562.H3k27ac.peak, file.k562.H3k27ac.bw )
  r1 <- calculate_K562_GM_matrix(tb.hist, "H3k27ac",  file.k562.H3k27ac.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k27me3.peak,file.k562.H3k27me3.bw )
  r2 <- calculate_K562_GM_matrix(tb.hist, "H3k27me3", file.k562.H3k27me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k36me3.peak,file.k562.H3k36me3.bw )
  r3 <- calculate_K562_GM_matrix(tb.hist, "H3k36me3", file.k562.H3k36me3.ctrl.bw, chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me1.peak, file.k562.H3k4me1.bw )
  r4 <- calculate_K562_GM_matrix(tb.hist, "H3k4me1",  file.k562.H3k4me1.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me2.peak, file.k562.H3k4me2.bw )
  r5 <- calculate_K562_GM_matrix(tb.hist, "H3k4me2",  file.k562.H3k4me2.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k4me3.peak, file.k562.H3k4me3.bw )
  r6 <- calculate_K562_GM_matrix(tb.hist, "H3k4me3",  file.k562.H3k4me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9ac.peak,  file.k562.H3k9ac.bw )
  r7 <- calculate_K562_GM_matrix(tb.hist, "H3k9ac",   file.k562.H3k9ac.ctrl.bw,   chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H3k9me3.peak, file.k562.H3k9me3.bw )
  r8 <- calculate_K562_GM_matrix(tb.hist, "H3k9me3",  file.k562.H3k9me3.ctrl.bw,  chr="chr22", range=10000*zoom_rate );
  add_spikes_to_removed_regs( file.k562.H4k20me1.peak,file.k562.H4k20me1.bw )
  r9 <- calculate_K562_GM_matrix(tb.hist, "H4k20me1", file.k562.H4k20me1.ctrl.bw, chr="chr22", range=10000*zoom_rate );

  df.K562.GM.allcomb <- rbind(r1, r2, r3, r4, r5, r6, r7,r8, r9);

  pdf("boxplot-histone-allcomb-mixed.pdf");
  boxplot(pearson~Mark,data=df.K562.GM.allcomb, main="(K562+GM) vs (K562+GM)",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5) 
  dev.off();
  
}