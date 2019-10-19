## 
## Draws meta plots near reigons where txn does or does NOT correctly predict DNase-1-seq signal.
##

source("metaPlot.R")
data_pth = "/fs/cbsudanko/storage/data/hg19/k562/"

## Graphing parameters
step_size = 10
halfWindow_size = 300

## Get BED coordinates for regions.
bf <- read.table("dnase.K562.peaks.CTCF.100-data.bed.gz")

## Generate meta plots for GRO-cap and PRO-seq.
pth = paste(data_pth, "groseq_tss/", sep="")
gc       = metaPlot_str(hd1, "groseq_tss_wTAP_plus.bigWig", "groseq_tss_wTAP_minus.bigWig", main="GROcap type1", path=pth, stp=step_size, halfWindow=halfWindow_size)
gc_noTap = metaPlot_str(hd1, "groseq_tss_noTAP_plus.bigWig", "groseq_tss_noTAP_minus.bigWig", main="GROcap noTAP type1", path=pth, stp=step_size, halfWindow=halfWindow_size)

pth = paste(data_pth, "proseq/", sep="")
ps  = metaPlot_str(hd1, "K562_unt.sort.bed.gz_plus.bw", "K562_unt.sort.bed.gz_minus.bw", main="PRO-seq type1", path=pth, stp=step_size, halfWindow=halfWindow_size)


## Generate meta plots for histone marks.
pth = paste(data_pth, "histones/", sep="")
h3k27ac = metaPlot(hd1, "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig", main="H3K27ac", path=pth, stp=step_size, halfWindow=halfWindow_size)
h3k27me3= metaPlot(hd1, "wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig", main="H3K27me3", path=pth, stp=step_size, halfWindow=halfWindow_size)
h3k4me3 = metaPlot(hd1, "wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig", main="H3K4me3", path=pth, stp=step_size, halfWindow=halfWindow_size)
h3k4me1 = metaPlot(hd1, "wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig", main="H3K4me1", path=pth, stp=step_size, halfWindow=halfWindow_size)
ctcf    = metaPlot(hd1, "wgEncodeBroadHistoneK562CtcfStdSig.bigWig", main="CTCF", path=pth, stp=step_size, halfWindow=halfWindow_size)

## Draw meta plots for each of these groups...
pdf("S2_Normalized_TSS.pdf")

plot(-100000, -1000000, xlim= c(-hW, hW), ylim=ylims, xlab="Distance to TSS", ylab="Signal [S2 normalized counts]")
add_data(LZ, "#cb6751", "#cb6751")
add_data(P, "#9e6ebd", "#9e6ebd")
add_data(D, "#7aa457", "#7aa457")

dev.off()



