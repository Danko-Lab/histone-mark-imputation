# PCA to visualize the F1 atlas.

require(bigWig)
step=100

## Read all of the data into a matrix.
tissues <- c("BN", "GI", "HT", "KD", "LG", "LV", "SK", "SP", "ST")
marks   <- c("H3k27ac.S1.V3.bw", "H3k27me3.S1.V3.bw", "H3k36me3.S1.V3.bw", "H3k4me1.S1.V2.bw", "H3k4me2.S1.V2.bw", "H3k4me3.S1.V3.bw", "H3k9ac.S1.V2.bw", "H3k9me3.S1.V2.bw", "H4k20me1.S1.V3.bw")
active_c <- c(1, 5, 6, 7) #c(1*c(1:9),5*c(1:9),6*c(1:9),7*c(1:9)) #K4me3, K4me2, K27ac, K9ac

nc <- c()
for(t in tissues) {
  for(m in marks) {
    ## Read tissue-mark.
    bw <- load.bigWig(paste("../pred_F1/",t,".",m,sep=""))
    nc <- cbind(nc, step.bpQuery.bigWig(bw, chrom="chr1", start=0, end=197195432, step=step))
    colnames(nc)[NCOL(nc)] <- paste(t,".",m,sep="")
  }
}

## Now do PCA.
pca <- prcomp(nc, center=TRUE, scale=TRUE)
pca_a<-prcomp(nc[,active_c])

cols <- c(rep("#5dab4a",9), rep("#b460bd",9), rep("#b2b248",9), rep("#6980ce",9), rep("#d68d46",9), rep("#4aac8d",9), rep("#cb5242",9), rep("#7c7130",9), rep("#c75980",9))
pch <- c(rep(c(1:9),9))
data.frame(colnames(nc), cols, pch)

pdf("PC1.PC2.pdf")
 plot(y= pca$rotation[,1], x= pca$rotation[,2], col=cols, pch=pch, xlab="PC2", ylab="PC1")
 legend("topleft", pch=c(1:9), legend=marks)
 pairs(pca$rotation[,1:5], col=cols, pch=pch)
dev.off()

## Do UMAP on the PCA output.
#hm.umap <- umap(nc, init="pca")
hm.umap <- umap(pca$rotation[,1:40])

pdf("F1_mouse.UMAP.pdf")
 plot(y= hm.umap$layout[,2], x= hm.umap$layout[,1], col=cols, pch=pch, xlab="UMAP1", ylab="UMAP2")
dev.off()


