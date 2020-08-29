# PCA to visualize the F1 atlas.

require(bigWig)
step=100
sz=3

## Read all of the data into a matrix.
tissues <- c("BN", "GI", "HT", "KD", "LG", "LV", "SK", "SP", "ST")
marks   <- c("H3k27ac.S1.V3.bw", "H3k27me3.S1.V3.bw", "H3k36me3.S1.V3.bw", "H3k4me1.S1.V2.bw", "H3k4me2.S1.V2.bw", "H3k4me3.S1.V3.bw", "H3k9ac.S1.V2.bw", "H3k9me3.S1.V2.bw", "H4k20me1.S1.V3.bw")
offs<-c(9 *c(0:8)); active_c <- c(1+offs, 5+offs, 6+offs, 7+offs) #c(1+9*c(1:9),5+9*c(1:9),6+9*c(1:9),7+9*c(1:9)) #K4me3, K4me2, K27ac, K9ac # 

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
 plot(y= pca$rotation[,1], x= pca$rotation[,2], col=cols, pch=pch, xlab="PC2", ylab="PC1", cex=sz)
 legend("topleft", pch=c(1:9), legend=marks)
 pairs(pca$rotation[,1:5], col=cols, pch=pch)

 plot(y= pca_a$rotation[,1], x= pca_a$rotation[,2], col=cols[active_c], pch=pch[active_c], xlab="PC2", ylab="PC1", cex=sz)
 pairs(pca_a$rotation[,1:5], col=cols[active_c], pch=pch[active_c])
dev.off()

## Do UMAP on the PCA output.
require(umap)
#hm.umap <- umap(nc, init="pca")
hm.umap <- umap(pca$rotation[,1:40])
hm.umap_a <- umap(pca_a$rotation[,1:20])

pdf("F1_mouse.UMAP.pdf")
 plot(y= hm.umap$layout[,1], x= hm.umap$layout[,2], col=cols, pch=pch, xlab="UMAP2", ylab="UMAP1", cex=sz)
 plot(y= hm.umap_a$layout[,1], x= hm.umap_a$layout[,2], col=cols[active_c], pch=pch[active_c], xlab="UMAP2", ylab="UMAP1", cex=sz)
dev.off()


