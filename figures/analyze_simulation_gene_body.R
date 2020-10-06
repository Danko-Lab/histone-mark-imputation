## 
## R script to analyze simulated data.

mark <- c("H3k4me1", "H3k27me3", "H3k36me3")
I <- c("1.0", "0.1", "0.01", "0.001")
R <- c("1.0", "0.1", "0.01", "0.001")

require(bigWig)
require(ComplexHeatmap)
require(circlize) # For colorRamp2

getData <- function(m, i, r) {
 filename <- paste("../pred-simulation/",m,".S1.V3.sim-",i,"-",r,"_chr6_ssto_hap7.bw", sep="")
 bw <- load.bigWig(filename)
 bed <- data.frame(chr= rep("chr6_ssto_hap7", 2), chromStart= c(1002000, 2002000), chromEnd= c(1025000, 2025000))
 data <- bed.region.bpQuery.bigWig(bw, bed)
 return( data[1] - data[2] ) ## Subtract a comparable sized background region.
}

alldata <- data.frame(expand.grid(I= I, R= R, mark = mark), signal= 0)

for(i in 1:NROW(alldata)) {
  alldata$signal[i] <- getData(alldata$mark[i], alldata$I[i], alldata$R[i])
}

myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(100)
myBreaks <- seq(0, 15000, length.out = 100)

## Is the slope higher for initiation or pausing?
ad <- alldata
ad$I <- as.double(as.character(alldata$I))
ad$R <- as.double(as.character(alldata$R))

#write.table(ad, "~/transfer/ad.rflat") ## For transfering to my PC
 
## Plot out scatterplots over I and R.
myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(NROW(I))

pdf("~/transfer/scatterplots-GB.pdf", height=15, width=5)
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


# write.table(ad, "~/transfer/ad.gb.rflat")
# ad <- read.table("ad.gb.rflat")


h3k36me3_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k36me3",])
h3k4me1_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k4me1",])
h3k27me3_r <- glm(scale(signal)~log(I)+log(R), data=ad[ad$mark == "H3k27me3",])


summary(h3k36me3_r)
summary(h3k4me1_r)
summary(h3k27me3_r)

confint(h3k36me3_r)
confint(h3k4me1_r)
confint(h3k27me3_r)

#Print out a figure!
## Key note: This was NOT working on CBSUGRU. I pulled the data into my PC and ran this code there.
require(jtools)

pdf("Coefficients.Fit-GB.pdf", width=4, height=2)
plot_summs(h3k4me1_r, h3k36me3_r, h3k27me3_r, plot.distributions = TRUE, model.names=c("H3k4me1", "H3k36me3", "H3k27me3"), colors=c("#b38807", "#47b304", "#4b00b3"))
dev.off()





