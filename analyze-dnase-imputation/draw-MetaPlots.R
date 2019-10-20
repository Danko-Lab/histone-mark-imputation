## 
## Draws meta plots near reigons where txn does or does NOT correctly predict DNase-1-seq signal.
##

source("/workdir/zw355/proj/prj15-histone/scripts/hist.svm.com.R"); # For write.tmp.bed(df.bed, compress=FALSE)

source("metaPlot.R")
data_pth = "/fs/cbsudanko/storage/data/hg19/k562/"

## Graphing parameters
step_size = 25
halfWindow_size = 2000
th = 0.25

## Get BED coordinates for regions.
## Columns represent: 
##    chr, start, stop, tbs0.idx, org, pred, ratio
bf_ctcf <- read.table("dnase.K562.peaks.CTCF.100-data.bed.gz")
colnames(bf_ctcf) <- c("chr", "start", "stop", "tbs0.idx", "org", "pred", "ratio")

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


####################################################
##
## Generate meta plots for GRO-cap and PRO-seq.

analyzeData <- function( hd1, main_title ) {

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

  ## Draw meta plots for each of these groups...
  pdf(paste(mainTitle, "_MetaPlots.pdf", sep=""))
  add_data_str(GC, "#cb6751", "#cb6751")
  add_data_str(GC_noTap, "#cb6751", "#cb6751")
  add_data_str(PS, "#cb6751", "#cb6751")

  plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims, xlab="Distance to TSS", ylab="Signal [Normalized counts]")
  add_data(dnase, "#cb6751", "#cb6751")
  
  plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims, xlab="Distance to TSS", ylab="Signal [Normalized counts]")
  add_data(h3k27ac, "#cb6751", "#cb6751")
  add_data(h3k27me3, "#9e6ebd", "#9e6ebd")
  add_data(h3k4me3, "#7aa457", "#7aa457")
  add_data(h3k4me1, "#7aa457", "#7aa457")
 
  plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims, xlab="Distance to TSS", ylab="Signal [Normalized counts]")
  add_data(ctcf, "#cb6751", "#cb6751")
 
  dev.off()

}

analyzeData(ur, "UpperRight-CTCF")
analyzeData(lr, "LowerRight-CTCF")

