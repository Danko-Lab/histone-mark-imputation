## 
## R script to analyze simulated data.

mark <- c("H3k27ac", "H3k9ac","H3k122ac", "H3k4me3", "H3k4me2", "H3k4me1", "H3k27me3")
I <- c("1.0", "0.1", "0.01", "0.001")
R <- c("1.0", "0.1", "0.01", "0.001")

require(bigWig)
require(ComplexHeatmap)
require(circlize) # For colorRamp2

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

## Is the slope higher for initiation or pausing?
ad <- alldata
ad$I <- as.double(as.character(alldata$I))
ad$R <- as.double(as.character(alldata$R))

#write.table(ad, "~/transfer/ad.rflat") ## For transfering to my PC
 
## Plot out scatterplots over I and R.
myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(NROW(I))

pdf("~/transfer/scatterplots.pdf", height=15, width=5)
 par(mfrow=c(7,2))
 #plot(ad$signal~ad$I, pch=19, col=ad$mark)
 #plot(ad$signal~ad$R, pch=19, col=ad$mark)

 for(m in mark) {

 xlim=c(min(ad$R), max(ad$R)); 
 ylim=c(min(ad$signal[ad$mark == m]), max(ad$signal[ad$mark == m]))
 plot(-1, -1, ylim=ylim, xlim=xlim, ylab=paste("Signal", m), xlab= "Pause Release Rate [1/(s, cell)]", log="x")
 count=1
 for(i in unique(ad$I)) {
  points(ad$signal[ad$I == i & ad$mark == m]~ad$R[ad$I == i & ad$mark == m], col= myCol[count], type="b")
  count=count+1
 }

 xlim=c(min(ad$I), max(ad$I));
 plot(-1, -1, ylim=ylim, xlim=xlim, ylab=paste("Signal", m), xlab= "Initiation Rate [1/(s, cell)]", log="x")
 count=1
 for(i in unique(ad$R)) {
  points(ad$signal[ad$R == i & ad$mark == m]~ad$I[ad$R == i & ad$mark == m], col= myCol[count], type="b")
  count=count+1
 }

 
 }

dev.off()

h3k4me3_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k4me3",])
h3k27ac_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k27ac",])
h3k9ac_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k9ac",])
h3k122ac_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k122ac",])
h3k4me2_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k4me2",])
h3k4me1_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k4me1",])
h3k27me3_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k27me3",])


summary(h3k4me3_r)
summary(h3k27ac_r)
summary(h3k9ac_r)
summary(h3k122ac_r)
summary(h3k4me2_r)
summary(h3k4me1_r)
summary(h3k27me3_r)


confint(h3k4me3_r)
confint(h3k27ac_r)
confint(h3k9ac_r)
confint(h3k122ac_r)
confint(h3k4me2_r)
confint(h3k4me1_r)
confint(h3k27me3_r)

#Print out a figure!
## Key note: This was NOT working on CBSUGRU. I pulled the data into my PC and ran this code there.
require(jtools)

pdf("Coefficients.Fit.pdf", width=4, height=2)
plot_summs(h3k4me1_r, h3k4me2_r, h3k4me3_r, h3k27me3_r, plot.distributions = TRUE, model.names=c("H3k4me1", "H3k4me2", "H3k4me3", "H3k27me3"), colors=c("#b38807", "#fb9a99", "#b31b1b", "#4b00b3"))
plot_summs(h3k27ac_r, h3k9ac_r, h3k122ac_r, plot.distributions = TRUE, model.names=c("H3k27ac", "H3k9ac", "H3k122ac"), colors=c("#007138", "#b30086", "#007171"))

plot_summs(h3k4me1_r, h3k4me2_r, h3k4me3_r, h3k27ac_r, h3k9ac_r, h3k122ac_r, plot.distributions = TRUE)
plot_summs(h3k4me3_r, h3k27ac_r, plot.distributions = TRUE, model.names=c("H3k4me3", "H3k27ac"))
dev.off()





