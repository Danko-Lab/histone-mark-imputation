##
## Creates heatmaps on human dREG-HD sites.

require(bigWig)
require(pheatmap)
require(RColorBrewer)

writeHeatmap<- function(bed, hMarkFile, name, hm_order= NULL, subs= NULL, breaks= NULL,
							cols= NULL, dist= 5000, step=25,
							path="/local/storage/data/hg19/cd4/epiRoadmap_histone/") {
	## Load mark.
	hMark <- load.bigWig(paste(path, hMarkFile, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")

	## Get a matrix of counts.
	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed, dist, dist), step=step, abs.value=TRUE)
	hmat <- matrix(unlist(hCountMatrix), nrow= NROW(bed), byrow=TRUE)# log(matrix(unlist(hCountMatrix), nrow= NROW(bed), byrow=TRUE)+1)
	if(is.null(hm_order)) {
	  hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -10):(NCOL(hmat)/2 +10)]), decreasing=TRUE)
	}
	hmat <- hmat[hm_order,]

	## Average by rows of 10.
	navg <- 20 ## Average every navg rows
	avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
	hmat <- avgMat

	## Write out a heatmap.
	if(is.null(breaks)) {
		bk <- seq(min(hmat), max(hmat), 0.01)
	} else {
		bk <- breaks
		hmat[hmat > max(breaks)] <- max(breaks)
	}

	if(is.null(cols)) {
		hmcols <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
	} else {
		hmcols <- colorRampPalette(cols)(length(bk)-1) # red
	}

	png(paste(name,".png",sep=""), width=400, height = 800)     # width and height are in pixels
		pheatmap( hmat, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE )
	dev.off()
	
	return(list(order= hm_order, breaks= bk))
}

writeColorScale<- function(breaks, cols, name) {
        ## Write out a heatmap.
        if(is.null(breaks)) {
                bk <- seq(min(hmat), max(hmat), 0.01)
        } else {
                bk <- breaks
        }

        if(is.null(cols)) {
                hmcols <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
        } else {
                hmcols <- colorRampPalette(cols)(length(bk)-1) # red
        }

	png(paste(name,".scale.png",sep=""), width=500, height=300)
		pheatmap(rev(bk), cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=TRUE, legend_breaks= quantile(bk), legend_labels= signif(quantile(bk),3), show_rownames=FALSE, show_colnames=FALSE)
	dev.off()
}

enhancers <- read.table("dreg.TF.dHS-.bed")
promoters <- read.table("dreg.TF.dHS+.bed")

hd <- rbind(enhancers, promoters)

ord   <- order(hd$V5, decreasing=TRUE)
ord_e <- order(enhancers$V5, decreasing=TRUE)
ord_p <- order(promoters$V5, decreasing=TRUE)

#sup <- writeHeatmap(hd[,1:5], "7895_5598_44886_HTYVLBGXY_K562_1_AGGCAGAA.bw", "ATAC-rep-1", hm_order= ord, path="./")
sup <- writeHeatmap(promoters[,1:5], "7895_5598_44886_HTYVLBGXY_K562_1_AGGCAGAA.bw", "ATAC-rep-1-DHS+", hm_order= ord_p, breaks= c(10, 50, 100, 200, seq(300, 1000, 10)), cols= brewer.pal(50, "Purples"), path="./", dist= 10000)
sup <- writeHeatmap(enhancers[,1:5], "7895_5598_44886_HTYVLBGXY_K562_1_AGGCAGAA.bw", "ATAC-rep-1-DHS-", hm_order= ord_e, breaks= sup$breaks, cols= brewer.pal(50, "Purples"), path="./", dist= 10000)
writeColorScale(sup$breaks, cols=brewer.pal(50, "Purples"), "ATAC-rep-1")

#sup <- writeHeatmap(hd[,1:5], "7895_5598_44887_HTYVLBGXY_K562_2_TCCTGAGC.bw", "ATAC-rep-2", hm_order= ord, path="./")
sup <- writeHeatmap(promoters[,1:5], "7895_5598_44887_HTYVLBGXY_K562_2_TCCTGAGC.bw", "ATAC-rep-2-DHS+", hm_order= ord_p, breaks= c(10, 50, 100, 200, seq(300, 1000, 10)), cols= brewer.pal(50, "Purples"), path="./", dist=10000)
sup <- writeHeatmap(enhancers[,1:5], "7895_5598_44887_HTYVLBGXY_K562_2_TCCTGAGC.bw", "ATAC-rep-2-DHS-", hm_order= ord_e, breaks= sup$breaks, cols= brewer.pal(50, "Purples"), path="./", dist=10000)
writeColorScale(sup$breaks, cols=brewer.pal(50, "Purples"), "ATAC-rep-2")

q("no")

sup <- writeHeatmap(hd[,1:5], "H3K27ac.bw", "H3K27ac", hm_order= ord)
sup <- writeHeatmap(enhancers, "H3K27ac.bw", "H3K27ac-e", hm_order= ord_e, breaks= sup$breaks)
sup <- writeHeatmap(promoters, "H3K27ac.bw", "H3K27ac-p", hm_order= ord_p, breaks= sup$breaks)
writeColorScale(sup$breaks, cols=NULL, "H3K27ac")

sup <- writeHeatmap(hd[,1:5], "H3K4me1.bw", "H3K4me1", hm_order= ord)
sup <- writeHeatmap(enhancers, "H3K4me1.bw", "H3K4me1-e", hm_order= ord_e, breaks= sup$breaks)
sup <- writeHeatmap(promoters, "H3K4me1.bw", "H3K4me1-p", hm_order= ord_p, breaks= sup$breaks)
writeColorScale(sup$breaks, cols=NULL, "H3K4me1")

sup <- writeHeatmap(hd[,1:5], "H3K4me3.bw", "H3K4me3", hm_order= ord)
sup <- writeHeatmap(enhancers, "H3K4me3.bw", "H3K4me3-e", hm_order= ord_e, breaks= sup$breaks)
sup <- writeHeatmap(promoters, "H3K4me3.bw", "H3K4me3-p", hm_order= ord_p, breaks= sup$breaks)
writeColorScale(sup$breaks, cols=NULL, "H3K4me3")

pth= "/local/storage/projects/NHP/AllData/All_Merge/"
sup <- writeHeatmap(hd[,1:5], "H-U_plus.bw", "PROseq.plus", hm_order= ord, cols=c("white","#fe0000"), path=pth)
sup <- writeHeatmap(enhancers, "H-U_plus.bw", "PROseq.plus-e", hm_order= ord_e, cols=c("white","#fe0000"), path=pth, breaks= sup$breaks)
sup <- writeHeatmap(promoters, "H-U_plus.bw", "PROseq.plus-p", hm_order= ord_p, cols=c("white","#fe0000"), path=pth, breaks= sup$breaks)
writeColorScale(sup$breaks, cols=c("white","#fe0000"), "PROseq.plus")

sup <- writeHeatmap(hd[,1:5], "H-U_minus.bw", "PROseq.minus", hm_order= ord, cols=c("white","#0000fe"), path=pth)
sup <- writeHeatmap(enhancers, "H-U_minus.bw", "PROseq.minus-e", hm_order= ord_e, cols=c("white","#0000fe"), path=pth, breaks= sup$breaks)
sup <- writeHeatmap(promoters, "H-U_minus.bw", "PROseq.minus-p", hm_order= ord_p, cols=c("white","#0000fe"), path=pth, breaks= sup$breaks)
writeColorScale(sup$breaks, cols=c("white","#0000fe"), "PROseq.minus")

pth= "/local/storage/data/hg19/cd4/dnase1fp/"
sup <- writeHeatmap(hd[,1:5], "dnase1.5prime.bw", "DNase1", hm_order= ord, path=pth)
sup <- writeHeatmap(enhancers, "dnase1.5prime.bw", "DNase1-e", hm_order= ord_e, breaks= sup$breaks, path=pth)
sup <- writeHeatmap(promoters, "dnase1.5prime.bw", "DNase1-p", hm_order= ord_p, breaks= sup$breaks, path=pth)
writeColorScale(sup$breaks, cols=NULL, "DNase1")


