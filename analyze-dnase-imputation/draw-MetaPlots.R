## 
## Draws meta plots near reigons where txn does or does NOT correctly predict DNase-1-seq signal.
##

source("/workdir/zw355/proj/prj15-histone/scripts/hist.svm.com.R"); # For write.tmp.bed(df.bed, compress=FALSE)

source("metaPlot.R")
data_pth = "/fs/cbsudanko/storage/data/hg19/k562/"

## Graphing parameters
stp = step_size = 25
hW = halfWindow_size = 2000
th = 0.25 

## Graphical paremeters
lw = 2 ## LWD in plots.

## Get BED coordinates for regions.
## Columns represent: 
##    chr, start, stop, tbs0.idx, org, pred, ratio
bf_ctcf <- read.table("dnase.K562.peaks.CTCF.100-data.bed.gz")
colnames(bf_ctcf) <- c("chr", "start", "stop", "tbs0.idx", "org", "pred", "ratio")

bf_peaks <- read.table("dnase.K562.peaks.100-data.bed.gz")
colnames(bf_peaks) <- c("chr", "start", "stop", "tbs0.idx", "org", "pred", "ratio")

bf_k27ac <- read.table("dnase.K562.peaks.K27ac.100-data.bed.gz")
colnames(bf_k27ac) <- c("chr", "start", "stop", "tbs0.idx", "org", "pred", "ratio")

## Santiy check the threshold
if(0) {
plot(bf_ctcf$org, bf_ctcf_pred)
abline(h=th, col="red")
abline(v=th, col="red")
}
## Looks resonable!

####################################################
##
## Pick the points of interest.

## Upper right
indx <- bf_ctcf$pred > th & bf_ctcf$org > th
file.tmp.bed <- write.temp.bed(bf_ctcf[indx,], compress=FALSE)
ur <- read.table(pipe(paste("bedtools merge -d 1 -i ",  file.tmp.bed)))#, " | sort-bed -")))

## Lower right
indx <- bf_ctcf$pred < th & bf_ctcf$org > th
file.tmp.bed <- write.temp.bed(bf_ctcf[indx,], compress=FALSE)
lr <- read.table(pipe(paste("bedtools merge -d 1 -i ",  file.tmp.bed)))#, " | sort-bed -")))

## K27ac
indx <- bf_k27ac$pred > th & bf_k27ac$org > th
file.tmp.bed <- write.temp.bed(bf_k27ac[indx,], compress=FALSE)
urk <- read.table(pipe(paste("bedtools merge -d 1 -i ",  file.tmp.bed)))#, " | sort-bed -")))

## Now for peaks...

## Upper right
indx <- bf_peaks$pred > th & bf_peaks$org > th
file.tmp.bed <- write.temp.bed(bf_peaks[indx,], compress=FALSE)
urp <- read.table(pipe(paste("bedtools merge -d 1 -i ",  file.tmp.bed)))#, " | sort-bed -")))

## Lower right
indx <- bf_peaks$pred < th & bf_peaks$org > th
file.tmp.bed <- write.temp.bed(bf_peaks[indx,], compress=FALSE)
lrp <- read.table(pipe(paste("bedtools merge -d 1 -i ",  file.tmp.bed)))#, " | sort-bed -")))


####################################################
##
## Generate meta plots for GRO-cap and PRO-seq.

analyzeData <- function( hd1, main_title, hW= halfWindow_size, ylims= NULL ) {

  ## Generate meta plots for run-on assays.
  pth = paste(data_pth, "groseq_tss/", sep="")
  GC       = metaPlot_str(hd1, "groseq_tss_wTAP_plus.bigWig", "groseq_tss_wTAP_minus.bigWig", main="GROcap type1", path=pth, stp=step_size, halfWindow=halfWindow_size)
  GC_noTap = metaPlot_str(hd1, "groseq_tss_noTAP_plus.bigWig", "groseq_tss_noTAP_minus.bigWig", main="GROcap noTAP type1", path=pth, stp=step_size, halfWindow=halfWindow_size)
 
  pth = paste(data_pth, "proseq/", sep="")
  PS  = metaPlot_str(hd1, "K562_unt.sort.bed.gz_plus.bw", "K562_unt.sort.bed.gz_minus.bw", main="PRO-seq type1", path=pth, stp=step_size, halfWindow=halfWindow_size)
 
  ## Generate meta plots for histone marks.
  pth = paste(data_pth, "histones/", sep="")
  h3k27ac = metaPlot(hd1, "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig", main="H3K27ac", path=pth, stp=step_size, halfWindow=halfWindow_size)
  h3k27me3= metaPlot(hd1, "wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig", main="H3K27me3", path=pth, stp=step_size, halfWindow=halfWindow_size)
  h3k4me3 = metaPlot(hd1, "wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig", main="H3K4me3", path=pth, stp=step_size, halfWindow=halfWindow_size)
  h3k4me1 = metaPlot(hd1, "wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig", main="H3K4me1", path=pth, stp=step_size, halfWindow=halfWindow_size)
  ctcf    = metaPlot(hd1, "wgEncodeBroadHistoneK562CtcfStdSig.bigWig", main="CTCF", path=pth, stp=step_size, halfWindow=halfWindow_size)

  ## Generate meta plots for DNase-I
  pth = paste(data_pth, "dnase/", sep="")
  dnase = metaPlot(hd1, "wgEncodeOpenChromDnaseK562SigV2.bigWig", main="DNase", path=pth, stp=step_size, halfWindow=halfWindow_size) 

  ## Now get the ylims if need be...
  if(is.null(ylims)) {
   ylims <- list( c(0,get_ylims_str(list(GC, GC_noTap, PS))), 
                  c(0,get_ylims(list(h3k27ac,h3k27me3,h3k4me3,h3k4me1))),
                  c(0,get_ylims(list(ctcf))),
                  c(0,get_ylims(list(dnase))) )
  }

  ## Draw meta plots for each of these groups...
#  pdf(paste(mainTitle, "_MetaPlots.pdf", sep=""))

  plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims[[1]], xlab="Distance to TSS", ylab="Signal [Normalized counts]", main = "Pol II")
  add_data_str(GC, "#ff6751", "#ff6751", lwd=lw)
  add_data_str(GC_noTap, "#006751", "#006751", lwd=lw)
  add_data_str(PS, "#000000", "#000000", lwd=lw)

  plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims[[2]], xlab="Distance to TSS", ylab="Signal [Normalized counts]", main = "Histone mods")
  add_data(h3k27ac, "#00aa00", lwd=lw)
  add_data(h3k27me3, "#00119e", lwd=lw)
  add_data(h3k4me3, "#7a0000", lwd=lw)
  add_data(h3k4me1, "#7a9999", lwd=lw)

  plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims[[3]], xlab="Distance to TSS", ylab="Signal [Normalized counts]", main = "CTCF")
  add_data(ctcf, "#000000", lwd=lw)

  plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims[[4]], xlab="Distance to TSS", ylab="Signal [Normalized counts]", main = "DNase-I")
  add_data(dnase, "#000000", lwd=lw)

#  dev.off()

  return(ylims)
}

pdf("CTCF_Sites.pdf")

par(mfrow = c(3,4))
yl <- analyzeData(ur, "UpperRight-CTCF")
yl <- analyzeData(lr, "LowerRight-CTCF", ylims= yl)
yl <- analyzeData(urk, "UpperRight-K27ac", ylims= yl)

dev.off()

pdf("DNase_Peaks.pdf")

par(mfrow = c(2,4))
yl <- analyzeData(urp, "UpperRight-CTCF")
yl <- analyzeData(lrp, "LowerRight-CTCF", ylims= yl)

dev.off()


