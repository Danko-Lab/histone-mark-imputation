library(dREG);
source("../script/hist.param.R");

## Read in refseq genes.
refGene <- read.table("/fs/cbsudanko/storage/projects/mcf7tamres/annotations/refGene.bed.gz")
refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene <- refGene[(refGene$V3-refGene$V2)>5000,]

bodies <- refGene
bodies <- bodies [bodies$V1!="chrY",]

bodies$V2[bodies$V6 == "+"] <-bodies$V2[bodies$V6 == "+"]+1000
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-1000

tss <- refGene
tss$V3[tss$V6 == "+"] <-tss$V2[tss$V6 == "+"]+1000
tss$V2[tss$V6 == "-"] <- tss$V3[tss$V6 == "-"]-1000

## Read in bigWigs.
require(bigWig)

H1_bw <- load.bigWig(paste(path.histone, file.H3k27ac.bw, sep="/"));
H2_bw <- load.bigWig(paste(path.histone, file.H3k9me3.bw, sep="/"));
H3_bw <- load.bigWig(paste(path.histone, file.H3k4me3.bw, sep="/"));
H4_bw <- load.bigWig(paste(path.histone, file.H3k4me1.bw, sep="/"));
H5_bw <- load.bigWig(paste(path.histone, file.H3k27me3.bw, sep="/"));
H6_bw <- load.bigWig(paste(path.histone, file.H4k20me1.bw, sep="/"));
H7_bw <- load.bigWig(paste(path.histone, file.H3k36me3.bw, sep="/"));
H8_bw <- load.bigWig(paste(path.histone, file.H3k4me2.bw, sep="/"));
H9_bw <- load.bigWig(paste(path.histone, file.H3k9ac.bw, sep="/"));
H10_bw <- load.bigWig(paste(path.histone, file.H3k122ac.bw, sep="/"));


## Count reads in each ...
H1  <- bed6.region.bpQuery.bigWig(H1_bw, H1_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H1_bw$basesCovered * H1_bw$mean ) * 1e6
H2  <- bed6.region.bpQuery.bigWig(H2_bw, H2_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H2_bw$basesCovered * H2_bw$mean ) * 1e6
H3  <- bed6.region.bpQuery.bigWig(H3_bw, H3_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H3_bw$basesCovered * H3_bw$mean ) * 1e6
H4  <- bed6.region.bpQuery.bigWig(H4_bw, H4_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H4_bw$basesCovered * H4_bw$mean ) * 1e6
H5  <- bed6.region.bpQuery.bigWig(H5_bw, H5_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H5_bw$basesCovered * H5_bw$mean ) * 1e6
H6  <- bed6.region.bpQuery.bigWig(H6_bw, H6_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H6_bw$basesCovered * H6_bw$mean ) * 1e6
H7  <- bed6.region.bpQuery.bigWig(H7_bw, H7_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H7_bw$basesCovered * H7_bw$mean ) * 1e6
H8  <- bed6.region.bpQuery.bigWig(H8_bw, H8_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H8_bw$basesCovered * H8_bw$mean ) * 1e6
H9  <- bed6.region.bpQuery.bigWig(H9_bw, H9_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H9_bw$basesCovered * H9_bw$mean ) * 1e6
H10 <- bed6.region.bpQuery.bigWig(H10_bw, H10_bw, bodies, abs.value = TRUE)/ 2 / (bodies$V3-bodies$V2) * 1000/ (H10_bw$basesCovered * H10_bw$mean ) * 1e6

HX <- cbind(H3k27ac=H1, H3k9me3=H2, H3k4me3=H3, H3k4me1=H4, H3k27me3=H5, H4k20me1=H6, H3k36me3=H7, H3k4me2=H8, H3k9ac=H9, H3k122ac=H10 )


if(1)
{
	source("/home/zw355/src/Rplot/levelplot2.R")
	source("/home/zw355/src/Rplot/drawCor.R")
	pdf("original-gene.pdf")
		drawCorWithPairs(HX, log=TRUE)
	dev.off()

	pdf("original-gene-strong.pdf")
		drawCorWithPairs( HX[,c("H3k4me3", "H3k9ac", "H3k4me1", "H3k27ac", "H3k4me2")] , log=TRUE)
	dev.off()

}