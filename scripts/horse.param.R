#path.histone    <- ""

file.horse.H3k27ac.bw   <- "H3K27ac_liver_merge.bw"
file.horse.H3k27me3.bw  <- "H3K27me3_liver_merge.bw"
file.horse.H3k4me1.bw   <- "H3K4me1_liver_merge.bw"
file.horse.H3k4me3.bw   <- "H3K4me3_liver_merge.bw"

file.horse.H3k27ac.peak <- "/local/workdir/zw355/proj/prj15-histone/horse_don/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz";
file.horse.H3k27me3.peak <-"/local/workdir/zw355/proj/prj15-histone/horse_don/wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak.gz";
file.horse.H3k4me1.peak <- "/local/workdir/zw355/proj/prj15-histone/horse_don/wgEncodeBroadHistoneK562H3k4me2StdPk.broadPeak.gz";
file.horse.H3k4me3.peak <- "/local/workdir/zw355/proj/prj15-histone/horse_don/wgEncodeBroadHistoneK562H3k4me3StdPk.broadPeak.gz";

file.horse.H3k27ac.bw2 <- ""
file.horse.H3k27me3.bw2 <- "wgEncodeUwHistoneK562H3k27me3StdRawRep1.bigWig"
file.horse.H3k4me3.bw2 <- "wgEncodeUwHistoneK562H3k4me3StdRawRep1.bigWig"
file.horse.H3k4me1.bw2 <- ""

file.proseq.horse <- c( "PROseq_liver_merge_R1.plus.bw", "PROseq_liver_merge_R1.minus.bw");

file.H3k27ac.positive  <- "H3K27ac_LiverAB_merge.narrowPeak.bed.gz"
file.H3k27ac.negative  <- "H3K27ac_LiverAB_negative.bed.gz"
file.H3k27me3.positive <- "H3K27me3_LiverAB_merge.narrowPeak.bed.gz"
file.H3k27me3.negative <- "H3K27me3_LiverAB_negative.bed.gz"
file.H3k4me3.positive  <- "H3K4me3_LiverAB_merge.narrowPeak.bed.gz"
file.H3k4me3.negative  <- "H3K4me3_LiverAB_negative.bed.gz"
file.H3k4me1.positive  <- "H3K4me1_LiverAB_merge.narrowPeak.bed.gz"
file.H3k4me1.negative  <- "H3K4me1_LiverAB_negative.bed.gz"

path.proseq <- "/local/workdir/zw355/proj/prj15-histone/horse_don/"
path.histone <- "/local/workdir/zw355/proj/prj15-histone/horse_don/"

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20) )
