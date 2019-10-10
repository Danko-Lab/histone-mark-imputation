require(bigWig)
options(scipen=999) # Prevent the use of scientific notation in output (not compatible with BED format)

writeBin <- function(file_exp_bw, file_imp_bw, prefix, th= 0.9) {
	## Read in bigWigs.
	m_exp <- load.bigWig(file_exp_bw)
	m_pred <- load.bigWig(file_imp_bw)

	## Now get histograms in 200 bp windows on chr1
	exp.sig <- list()
	for(i in 1:NROW(m_pred$chroms)) { ## NOTE the use of m_pred here. Exp contains chrX. Here we're just looking for autosomes.
	exp.sig[[m_exp$chroms[i]]] <- step.bpQuery.bigWig(m_exp, m_exp$chroms[i], 0, m_exp$chromSizes[i], 200)
	}

	imp.sig <- list()
	for(i in 1:NROW(m_pred$chroms)) {
	imp.sig[[m_pred$chroms[i]]] <- step.bpQuery.bigWig(m_pred, m_pred$chroms[i], 0, m_pred$chromSizes[i], 200)
	}

	## Get top x%
	th.exp <- quantile(unlist(exp.sig), th)
	th.imp <- quantile(unlist(imp.sig), th)

	bed.exp <- data.frame()
	for(i in 1:NROW(exp.sig)){
	bed.exp <- rbind(bed.exp, data.frame(chrom= rep(names(exp.sig)[i]), chromStart= ((which(exp.sig[[i]] > th.exp)-1)* 200),  chromEnd= ((which(exp.sig[[i]] > th.exp))* 200)))
	}

	bed.imp <- data.frame()
	for(i in 1:NROW(imp.sig)){
	bed.imp <- rbind(bed.imp, data.frame(chrom= rep(names(imp.sig)[i]), chromStart= ((which(imp.sig[[i]] > th.imp)-1)* 200),  chromEnd= ((which(imp.sig[[i]] > th.imp))* 200)))
	}

	write.table(bed.exp, file = pipe(paste("sort-bed - | bedops --merge - | bgzip > ",prefix,".exp.bed.gz", sep="")), row.names=FALSE, sep="\t", quote=FALSE, col.names=FALSE)
	write.table(bed.imp, file = pipe(paste("sort-bed - | bedops --merge - | bgzip > ",prefix,".imp.bed.gz", sep="")), row.names=FALSE, sep="\t", quote=FALSE, col.names=FALSE)
}

file_K562_H3k27ac_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig"
file_K562_H3k27me3_bw = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig"
file_K562_H3k36me3_bw = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig"
file_K562_H3k9me3_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig"
file_K562_H3k4me1_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig"
file_K562_H3k4me3_bw  = "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig"
file_GM_H3k27ac_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig"
file_GM_H3k27me3_bw = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27me3StdSig.bigWig"
file_GM_H3k36me3_bw = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k36me3StdSig.bigWig"
file_GM_H3k9me3_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k9me3StdSig.bigWig"
file_GM_H3k4me1_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me1StdSig.bigWig"
file_GM_H3k4me3_bw  = "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me3StdSig.bigWig"

file_K562_H3k27ac_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-k562/raw.H3k27ac.S1.V3.G1.bw"
file_K562_H3k27me3_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-k562/raw.H3k27me3.S1.V3.G1.bw"
file_K562_H3k36me3_pred_raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-k562/raw.H3k36me3.S1.V3.G1.bw"
file_K562_H3k9me3_pred_raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-k562/raw.H3k9me3.S1.V2.G1.bw"
file_K562_H3k4me1_pred_raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-k562/raw.H3k4me1.S1.V2.G1.bw"
file_K562_H3k4me3_pred_raw  <- "/local/workdir/zw355/proj/prj15-histone/pred-k562/raw.H3k4me3.S1.V3.G1.bw"
file_GM_H3k4me1_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-gm/raw.H3k4me1.S1.V2.GM.bw"
file_GM_H3k4me3_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-gm/raw.H3k4me3.S1.V3.GM.bw"
file_GM_H3k9me3_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-gm/raw.H3k9me3.S1.V2.GM.bw"
file_GM_H3k36me3_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-gm/raw.H3k36me3.S1.V3.GM.bw"
file_GM_H3k27me3_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-gm/raw.H3k27me3.S1.V3.GM.bw"
file_GM_H3k27ac_pred_raw <- "/local/workdir/zw355/proj/prj15-histone/pred-gm/raw.H3k27ac.S1.V3.GM.bw"

thBroad = 0.9
thNarrow = 0.99

writeBin(file_K562_H3k27ac_bw, file_K562_H3k27ac_pred_raw, "k562.H3k27ac", thNarrow)
writeBin(file_K562_H3k4me3_bw, file_K562_H3k4me3_pred_raw, "k562.H3k4me3", thNarrow)
writeBin(file_K562_H3k27me3_bw, file_K562_H3k27me3_pred_raw, "k562.H3k27me3", thBroad)
writeBin(file_K562_H3k36me3_bw, file_K562_H3k36me3_pred_raw, "k562.H3k36me3", thBroad)
writeBin(file_K562_H3k4me1_bw, file_K562_H3k4me1_pred_raw, "k562.H3k4me1", thBroad)
writeBin(file_K562_H3k9me3_bw, file_K562_H3k9me3_pred_raw, "k562.H3k9me3", thNarrow)

writeBin(file_GM_H3k27ac_bw, file_GM_H3k27ac_pred_raw, "GM.H3k27ac", thNarrow)
writeBin(file_GM_H3k4me3_bw, file_GM_H3k4me3_pred_raw, "GM.H3k4me3", thNarrow)
writeBin(file_GM_H3k27me3_bw, file_GM_H3k27me3_pred_raw, "GM.H3k27me3", thBroad)
writeBin(file_GM_H3k36me3_bw, file_GM_H3k36me3_pred_raw, "GM.H3k36me3", thBroad)
writeBin(file_GM_H3k4me1_bw, file_GM_H3k4me1_pred_raw, "GM.H3k4me1", thBroad)
writeBin(file_GM_H3k9me3_bw, file_GM_H3k9me3_pred_raw, "GM.H3k9me3", thNarrow)
