##
## getCounts.R - Counts reads in each gene.
require(bigWig)

prefix <- "groseq_tss_wTAP"
tus <- read.table("/home/cgd24/cbsudanko/data/hg19/gm12878/groseq_tss/peaks/hg19.gm12878.new_hmm2b.post2.bed", header=FALSE)

countBigWig <- function(prefix, bed, rpkm=TRUE, path="/home/cgd24/cbsudanko/data/hg19/gm12878/groseq_tss/") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bigWig", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bigWig", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

## Gets RPKMs
rpkm <- countBigWig(prefix, tus, rpkm=TRUE)
tus[,5] <- rpkm


options(scipen=999)
write.table(tus, "hg19.gm12878.new_hmm2b.post2.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


