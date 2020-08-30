## 
## R script to analyze simulated data.

mark <- c("H3k27ac", "H3k4me3")
I <- c("1.0", "0.1", "0.01", "0.001")
R <- c("1.0", "0.1", "0.01", "0.001")

require(bigWig)
require(ComplexHeatmap)

getData <- function(m, i, r) {
 filename <- paste("../pred-simulation/",m,".S1.V3.sim-",i,"-",r,"_chr6_ssto_hap7.bw", sep="")
 bw <- load.bigWig(filename)
 bed <- data.frame(chr= rep("chr6_ssto_hap7", 2), chromStart= c(998000, 1998000), chromEnd= c(1002000, 2002000))
 data <- bed.region.bpQuery.bigWig(bw, bed)
 return( data[1] - data[2] ) ## Subtract a comparable sized background region.
}


alldata <- data.frame(expand.grid(I= I, R= R, mark = mark), signal= 0)

for(i in 1:NROW(alldata)) {
  alldata$signal[i] <- getData(alldata$mark[i], alldata$I[i], alldata$R[i])
}

mat_k27ac <- matrix(0, nrow = NROW(I), ncol = NROW(R), dimnames = list(I, rev(R)))
m <- "H3k27ac"
for(i in 1:NROW(I)) {
 for(r in 1:NROW(R)) {
  mat_k27ac[i, NROW(R) - r + 1] <- alldata$signal[alldata$I == I[i] & alldata$R == R[r] & alldata$mark == m]
 }
}

mat_k4me3 <- matrix(0, nrow = NROW(I), ncol = NROW(R), dimnames = list(I, rev(R)))
m <- "H3k4me3"
for(i in 1:NROW(I)) {
 for(r in 1:NROW(R)) {
  mat_k4me3[i, NROW(R) - r + 1] <- alldata$signal[alldata$I == I[i] & alldata$R == R[r] & alldata$mark == m]
 }
}


myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(100)
myBreaks <- seq(0, 15000, length.out = 100)

pdf("~/transfer/heatmaps.pdf")
Heatmap(mat_k27ac, cluster_rows = FALSE, cluster_columns = FALSE,
    col = colorRamp2(myBreaks, myCol),
    row_title = 'Initiation Rate [1/(s, cell)]',
    row_names_side = 'left', 
    column_title = 'Pause Release Rate [1/(s, cell)]',
    column_title_side = 'bottom', 
    heatmap_legend_param = list(
      title = 'H3k27ac',
      color_bar = 'continuous',
      legend_direction = 'vertical',
      legend_width = unit(8, 'cm'),
      legend_height = unit(5.0, 'cm'),
      title_position = 'topcenter',
      title_gp=gpar(fontsize = 30, fontface = 'bold'),
      labels_gp=gpar(fontsize = 24, fontface = 'bold'))
 )

Heatmap(mat_k4me3, cluster_rows = FALSE, cluster_columns = FALSE,
    col = colorRamp2(myBreaks, myCol),
    row_title = 'Initiation Rate [1/(s, cell)]',
    row_names_side = 'left',
    column_title = 'Pause Release Rate [1/(s, cell)]',
    column_title_side = 'bottom',
    heatmap_legend_param = list(
      title = 'H3k4me3',
      color_bar = 'continuous',
      legend_direction = 'vertical',
      legend_width = unit(8, 'cm'),
      legend_height = unit(5.0, 'cm'),
      title_position = 'topcenter',
      title_gp=gpar(fontsize = 30, fontface = 'bold'),
      labels_gp=gpar(fontsize = 24, fontface = 'bold'))
 )

dev.off()

