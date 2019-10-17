## 
## Draws meta plots near reigons where txn does or does NOT correctly predict DNase-1-seq signal.
##

load("metaPlot.R")

data_pth = ""

## Get BED coordinates for regions.
bf <- read.table("dnase.K562.peaks.CTCF.100-data.bed.gz")

## Generate meta plots for GRO-cap and PRO-seq.
pth = paste(data_pth, "groseq_tss/", sep="")
metaPlot_str(hd1, "groseq_tss_wTAP_plus.bigWig", "groseq_tss_wTAP_minus.bigWig", main="GROcap type1", path=pth, stp=10, halfWindow=300)
metaPlot_str(hd1, "groseq_tss_noTAP_plus.bigWig", "groseq_tss_noTAP_minus.bigWig", main="GROcap noTAP type1", path=pth, stp=10, halfWindow=300)

pth = paste(data_pth, "proseq/", sep="")
metaPlot_str(hd1, "K562_unt.sort.bed.gz_plus.bw", "K562_unt.sort.bed.gz_minus.bw", main="PRO-seq type1", path=pth, stp=10, halfWindow=300)


## Generate meta plots for histone marks.
pth = paste(data_pth, "histones/", sep="")
metaPlot(hd1, "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig", main="H3K27ac type1", path=pth, stp=10, halfWindow=300)



