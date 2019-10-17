require(bigWig)

metaPlot <- function(bed, HP_FILE, path="./", stp= 100, halfWindow= 10000, ...) {
	bed <- center.bed(bed, halfWindow, halfWindow)
	HP <- load.bigWig(paste(path, HP_FILE, sep=""))
	H_meta <- metaprofile.bigWig(bed, HP, step=stp)

        N = length(H_meta$middle)
        x = 1:N*stp 
	
	plot.metaprofile(H_meta, X0= 0, ...)
	abline(v=halfWindow, lty="dotted")
}

metaPlot_str <- function(bed, HP, HM, path="./", stp= 100, halfWindow= 10000, ...) {
        bed <- center.bed(bed, halfWindow, halfWindow)

        HP <- load.bigWig(paste(path, HP, sep=""))
        HM <- load.bigWig(paste(path, HM, sep=""))

        H_meta_p <- metaprofile.bigWig(bed, HP, HM, step=stp)
        H_meta_m <- metaprofile.bigWig(bed, HM, HP, step=stp)

        N = length(H_meta_p$middle)
        x = 1:N*stp ## ((1:N) - N/2)* stp

        plot.metaprofile(H_meta_p, minus.profile=H_meta_m, X0= 0, ...)
        abline(v=halfWindow, lty="dotted")
}

if(0) { ## Examples below. if(0) designed to present this from runnig. 

 ## Load dREG-HD sites.
 hd  <- read.table("Bickmore.enhancers.hg18.bed.gz"); ## Actually ... it's hg19
 hd1 <- hd[grep("grp1",hd$V8),1:3]
 hd2 <- hd[grep("grp2",hd$V8),1:3]
 hd3 <- hd[grep("grp3",hd$V8),1:3]
 hd4 <- hd[grep("grp4",hd$V8),1:3]
 
 
 pdf("Marks.Enh.1.2.pdf")
 
 dist <- 10000
 
 pth= "/local/storage/data/hg19/k562/histones/"
 metaPlot(hd1, "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig", main="H3K27ac type1", path=pth, stp=10, halfWindow=300)
 metaPlot(hd2, "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig", main="H3K27ac type2", path=pth, stp=10, halfWindow=300)

 pth= "/local/storage/data/hg19/k562/groseq_tss/"
 metaPlot_str(hd1, "groseq_tss_wTAP_plus.bigWig", "groseq_tss_wTAP_minus.bigWig", main="GROcap type1", path=pth, stp=10, halfWindow=300)
 metaPlot_str(hd1, "groseq_tss_noTAP_plus.bigWig", "groseq_tss_noTAP_minus.bigWig", main="GROcap noTAP type1", path=pth, stp=10, halfWindow=300)
 metaPlot_str(hd2, "groseq_tss_wTAP_plus.bigWig", "groseq_tss_wTAP_minus.bigWig", main="GROcap type2", path=pth, stp=10, halfWindow=300)
 metaPlot_str(hd2, "groseq_tss_noTAP_plus.bigWig", "groseq_tss_noTAP_minus.bigWig", main="GROcap noTAP type2", path=pth, stp=10, halfWindow=300)

 pth= "/local/storage/data/hg19/k562/proseq/"
 metaPlot_str(hd1, "K562_unt.sort.bed.gz_plus.bw", "K562_unt.sort.bed.gz_minus.bw", main="PRO-seq type1", path=pth, stp=10, halfWindow=300)
 metaPlot_str(hd2, "K562_unt.sort.bed.gz_plus.bw", "K562_unt.sort.bed.gz_minus.bw", main="PRO-seq type2", path=pth, stp=10, halfWindow=300)
 
 dev.off()
}
