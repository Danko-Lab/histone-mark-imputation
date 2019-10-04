# /local/storage/projects/prophaseI/metaplots/ScaledMetaPlotFunctions.r

library(bigWig)
library(lattice)

rowMins <- function(x) {
  sapply(1:NROW(x), function(i) {min(x[i,])})
}
rowMaxs <- function(x) {
  sapply(1:NROW(x), function(i) {max(x[i,])})
}
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}
rowWMeans <- function(x, w) {
  sapply(1:NROW(x), function(i) {weighted.mean(x[i,], w)})
}

#############################################################################################################
## Functions for making scaled meta plots

# mround is used to round the length (in bases) to a number divisible by the 'base' argument
mround = function(x, base = 500){ 
  #base*ceiling(x/base)# base = 30 was chosen so that there will be a total of 60 windows (30 on each side of the midpoint of scaled window)
  base*floor(x/base)
} # rather than rounding up, floor can be used to round down and prevent enchroaching on the unscaled region 

# scale_params assesses the length of your gene, returning a buffer length (1/2 the size of the region to be scaled), a step (window size for scaled region), and a scale facor (windowsize/10)
scale_params = function(start, end, step=50, buffer=2000){
  length = end-start
  if (length <1000){
    buffer =  buffer
    step1  = step
    scaleFact = 1
  }
  else{
    length = ((end-start)-2*buffer)/2
    buffer = mround(length, base = 500)
    step1 = (buffer/500)
    scaleFact = step1/step # scale factor is divided by your scaled window count to adjust as if it were a 10bp window (i.e. window size = 5bp will have a scale factor of 2)
  }
  return(c(buffer,step1,scaleFact))
}

collect.many.scaled = function (bed, bigWig.plus, bigWig.minus, flank=10000, FLbuffer = 2000, step=50, 
                                do.sum = T) 
{
  midPoint = (bed[, 2] + bed[, 3])/2
  TSS = bed[,2]
  END = bed[,3]
  prstart = (TSS - flank)
  prend = (TSS + FLbuffer)
  finstart = (END - FLbuffer)
  finend = (END + flank)
  windowSize = 2*((flank + FLbuffer)/step)+1000
  #windowSize = 320
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      midscaleParams = scale_params(bed[i,2], bed[i,3], step, buffer = FLbuffer)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint[i] - midbuf)
      midend = (midPoint[i] + midbuf)
      prrow = step.bpQuery.bigWig(bigWig, chrom, prstart[i], prend[i], step, strand= "+") 
      midrow = step.bpQuery.bigWig(bigWig, chrom, midstart, midend, midstep, strand= "+")/ midScaleFact
      endrow = step.bpQuery.bigWig(bigWig, chrom, finstart[i], finend[i], step, strand="+")
      row = c(prrow, midrow, endrow)
      if (length(row) == windowSize){
        result[i, ] = abs(row)
      }
    }
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      midscaleParams = scale_params(bed[i,2], bed[i,3], step, buffer = FLbuffer)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint[i] - midbuf)
      midend = (midPoint[i] + midbuf)
      prrow = step.bpQuery.bigWig(bigWig, chrom, prstart[i], prend[i], step, strand="-", abs.value=TRUE)
      midrow =  step.bpQuery.bigWig(bigWig, chrom, midstart, midend, midstep, strand="-", abs.value=TRUE)/midScaleFact
      endrow = step.bpQuery.bigWig(bigWig, chrom, finstart[i], finend[i], step, strand="-", abs.value=TRUE)
      row = c(prrow, midrow, endrow)
      if (length(row) == windowSize){
        result[i, ] = abs(rev(row))
      }
    }
  }
  rownames(result) = bed[,4]
  result[complete.cases(result), ]
}

meta.subsample.scaled = function (bed, bigWig.plus, bigWig.minus, flank=10000, buffer = 2000, step=50, do.sum = T) 
{
  values = collect.many.scaled(bed, bigWig.plus, bigWig.minus, flank=flank, buffer, step, 
                               do.sum = do.sum)
  N = dim(values)[1]
  nPermut = 1000
  sampleFrac = 0.1
  #windowSize = ((prend-prstart)+(midend-midstart)+(finend-finstart))/step
  #windowSize = 320
  windowSize = 2*(flank + buffer)/step + 1000
  result = matrix(nrow = nPermut, ncol = windowSize)
  M = as.integer(round(N * sampleFrac, 0))
  #values = collect.many.scaled(bed, bigWig.plus, bigWig.minus, flank=flank, buffer, step, 
  #do.sum = do.sum)
  for (i in 1:nPermut) {
    idx <- sample(N, size = M, replace = T)
    result[i, ] = colSums(values[idx, ])/M
  }
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.875))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.125))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.5))
  return(list(result, ci9, ci1, ci5))
}
#####################################################################################################################
## functions for reading in a list of bigwig files and formatting for plotting. 
#example of bigwig list
#wigset = rbind(c("sample1_plus.bw", "sample1_minus.bw", "sample1", spikeCounts1/100000),
              # c("sample2_plus.bw", "sample2_minus.bw", "sample2", spikeCounts2/100000)
## function for loading a list of bigwig files (formatted as above with 4 columns per sample)
load.wigset <- function(wigset, wigset_row) {
  file = wigset[wigset_row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(file)
  file = wigset[wigset_row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(file)
  return(list(wig.p, wig.m, wigset[wigset_row, 3], wigset[wigset_row, 4]))
}

meta.normalize <- function(result, scaleFactor) {
  lapply(result, function(res) res * scaleFactor)
}

## formatting Scaled meta plots for lattice
FormatScaledMetaData_lattice <- function(wigset, bed){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 9, nrow = 0), stringsAsFactors = F)
  sampleNames = c()
  for (i in 1:N){
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample scaled meta-plot data ...\n")
    meta= meta.subsample.scaled(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], g_flank, g_buffer, g_step, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1)#/as.numeric(wigs[[4]]))
    sampleNames = c(sampleNames, wigs[[3]])
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    xAxis = seq(from = 1, to = dim(meta[[1]])[2])
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample")
  return(df)
}

#####################################################################################################################
## plotting profiles with lattice
### Note: lattice_Scaledmeta.proSeq can/should be adjusted based on desired grouping factors. 
require(lattice)
my.panel.bands <-
  function(x, y, upper, lower,
           fill, col,
           subscripts, ..., font, fontface)
  {
    upper <- upper[subscripts]
    lower <- lower[subscripts]
    panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                  col = fill, border = FALSE,
                  ...)
  }

analyze <- function(bwp, bwm, cat="data") {

  HP <- load.bigWig( bwp )
  HM <- load.bigWig( bwm )
#  size <- HP$mean*HP$basesCovered + abs(HM$mean*HM$basesCovered)

  # Gene body normalize.
  body <- bed; body$TXSTART[bed$TXSTRAND == "+"] = body$TXSTART[bed$TXSTRAND == "+"] + 500; body$TXEND[bed$TXSTRAND == "-"] = body$TXEND[bed$TXSTRAND == "-"] -500
  size <- sum(bed6.region.bpQuery.bigWig(HP, HM, body, abs.value=TRUE))

  wigset = rbind(c(bwp, bwm, cat, size/1e6))
  data <- FormatScaledMetaData_lattice(wigset, bed)

  return(data)
}


lattice_Scaledmeta.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                                     xlab = "Distance to gene boundaries (bp)",  ylab = "Median PRO-seq intensity"){
  pdf(paste(fig_dir, filename), width = 35, height = 10)
  result <- xyplot(mean ~ x, data = df,
                  group = factor(sample, labels = histone.family), 
                  scales = list(tck=c(1,0),alternating = c(1,1),
                                 x=list(relation='free',axs='i', labels = c('-5000', 'TSS', "+2000", "-2000", "CPS", "30000"), at = c(1, 601, 641, 1641, 1681, 2280)),
#                                y=list(log = 10, equispaced.log = FALSE)), ##Logscale
                                 y=list(relation='free',axs='i')), ##Nonlog

                   key=list( corner=c(0.98,0.85),
                             padding.text=3,
                            text=list(histone.family),col=rep("black", 9), cex=0.8, font=2),
                            rectangles=list(col=histone.color),#, rgb(0,0.5,0.5), rgb(0.5,0,0.5)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   xlim = c(500, 2280),
                   col = histone.color, # rgb(0,0.5,0.5,0.7), rgb(0.5,0,0.5,0.7)),  ## just add or take away colors to match number of samples in each panel
                   fill = histone.color, #, rgb(0,0.5,0.5,0), rgb(0.5,0,0.5,0)),
                   ylab = ylab,
                   xlab = xlab,
                   main = main,
                   aspect=0.5,
                   lwd=3, ##
                   par.settings = list(strip.background=list(col="lightgrey"),par.xlab.text=list(cex=1.1,font=2),par.ylab.text=list(cex=1.1,font=2),axis.text=list(cex=1)),
                   par.strip.text=list(cex=1.2,font=2),
                   upper = df$upper, ##Nonlog
                   lower = df$lower, ##Nonlog
#                   upper = log(df$upper,10), ##Logscale
#                   lower = log(df$lower,10), ##Logscale
                   panel = function(x, y, ...){
                     panel.grid(h=-1, v=-1)
                     panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     panel.xyplot(x, y, ...)
                   }) + latticeExtra::layer(panel.abline(v = c(601,641,1641,1681), lty = 2, lwd = 2, col = 'gray'))
  print(result)
  dev.off()
  print(result)
}


## Globals
fig_dir = "./"

g_flank= 30000 #10000
g_buffer=2000
g_step=50

## Create BED file.
bed <- read.table("/fs/cbsudanko/storage/projects/prophaseI/annotations/tuSelecter/final_tus.txt", header=TRUE);
bed <- bed[(bed$TXEND - bed$TXSTART)>10000,]
bed <- bed[bed$TXTYPE=="protein_coding",]
bed <- bed[!bed$TXNAME=="Untranscribed",]
exclude_chromosome = "_|chrM|chrY|chrX|rRNA"; 
bed <- bed [ grep( exclude_chromosome , bed$TXCHROM, invert=TRUE), ]

#bed <- bed[bed$TXSTRAND == "+",]

source("../scripts/hist.param.R");
## Create wigset for all nuclei samples.
H3k27ac    = analyze( file.gm.H3k27ac.bw,  file.gm.H3k27ac.bw );
H3k27me3   = analyze( file.gm.H3k27me3.bw, file.gm.H3k27me3.bw);
H3k36me3   = analyze( file.gm.H3k36me3.bw, file.gm.H3k36me3.bw);
H3k4me1    = analyze( file.gm.H3k4me1.bw,  file.gm.H3k4me1.bw );
H3k4me2    = analyze( file.gm.H3k4me2.bw,  file.gm.H3k4me2.bw );
H3k4me3    = analyze( file.gm.H3k4me3.bw,  file.gm.H3k4me3.bw );
H3k9ac     = analyze( file.gm.H3k9ac.bw,   file.gm.H3k9ac.bw  );
H3k9me3    = analyze( file.gm.H3k9me3.bw,  file.gm.H3k9me3.bw );
H4k20me1   = analyze( file.gm.H4k20me1.bw, file.gm.H4k20me1.bw);

PH3k27ac    = analyze( file.gm.H3k27ac.pred.bw,  file.gm.H3k27ac.pred.bw ); PH3k27ac[,c(2:4)]  <- PH3k27ac[,c(2:4)]*10
PH3k27me3   = analyze( file.gm.H3k27me3.pred.bw, file.gm.H3k27me3.pred.bw); PH3k27me3[,c(2:4)] <- PH3k27me3[,c(2:4)]*10
PH3k36me3   = analyze( file.gm.H3k36me3.pred.bw, file.gm.H3k36me3.pred.bw); PH3k36me3[,c(2:4)] <- PH3k36me3[,c(2:4)]*10
PH3k4me1    = analyze( file.gm.H3k4me1.pred.bw,  file.gm.H3k4me1.pred.bw ); PH3k4me1[,c(2:4)]  <- PH3k4me1[,c(2:4)]*10
PH3k4me2    = analyze( file.gm.H3k4me2.pred.bw,  file.gm.H3k4me2.pred.bw ); PH3k4me2[,c(2:4)]  <- PH3k4me2[,c(2:4)]*10
PH3k4me3    = analyze( file.gm.H3k4me3.pred.bw,  file.gm.H3k4me3.pred.bw ); PH3k4me3[,c(2:4)]  <- PH3k4me3[,c(2:4)]*10
PH3k9ac     = analyze( file.gm.H3k9ac.pred.bw,   file.gm.H3k9ac.pred.bw  ); PH3k9ac[,c(2:4)]   <- PH3k9ac[,c(2:4)]*10
PH3k9me3    = analyze( file.gm.H3k9me3.pred.bw,  file.gm.H3k9me3.pred.bw ); PH3k9me3[,c(2:4)]  <- PH3k9me3[,c(2:4)]*10
PH4k20me1   = analyze( file.gm.H4k20me1.pred.bw, file.gm.H4k20me1.pred.bw); PH4k20me1[,c(2:4)] <- PH4k20me1[,c(2:4)]*10


histone.family <- c("H3k27ac", "H3k27me3", "H3k36me3", "H3k4me1","H3k4me2","H3k4me3", "H3k9ac","H3k9me3","H4k20me1")
histone.color  <- c( rgb(0,0.25,0.25), rgb(0,0,0.25), rgb(0.25,0,0), rgb(0,0.5,0.5), rgb(0,0,0.5), rgb(0.5,0,0), rgb(0,0.75,0.75), rgb(0,0,0.75), rgb(0.75,0,0) )

save.image("ScaledMetaPlotFunctions.RData")

## Combine for plotting.
scaled_df_H3k27ac  <- rbind(H3k27ac,  PH3k27ac );  scaled_df_H3k27ac[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H3k27me3 <- rbind(H3k27me3, PH3k27me3 ); scaled_df_H3k27me3[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H3k36me3 <- rbind(H3k36me3, PH3k36me3);  scaled_df_H3k36me3[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H3k4me1  <- rbind(H3k4me1,  PH3k4me1 );  scaled_df_H3k4me1[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H3k4me2  <- rbind(H3k4me2,  PH3k4me2);   scaled_df_H3k4me2[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H3k4me3  <- rbind(H3k4me3,  PH3k4me3 );  scaled_df_H3k4me3[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H3k9ac   <- rbind(H3k9ac,   PH3k9ac  );  scaled_df_H3k9ac[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H3k9me3  <- rbind(H3k9me3,  PH3k9me3 );  scaled_df_H3k9me3[,5] <- c(rep(1,nelem), rep(2,nelem) );
scaled_df_H4k20me1 <- rbind(H4k20me1, PH4k20me1 ); scaled_df_H4k20me1[,5] <- c(rep(1,nelem), rep(2,nelem) );

#nelem <- NROW(scaled_df)/9 ## Since computed separately, have to change the categorical variable for groups.
#scaled_df[,5] <- c(rep(1,nelem), rep(2,nelem), rep(3,nelem), rep(4,nelem), rep(5,nelem), rep(6,nelem), rep(7,nelem), rep(8,nelem), rep(9,nelem))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-.25.pdf", scaled_df, ylim=c(0.0,.25))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-.5.pdf", scaled_df, ylim=c(0.0,.5))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-3.pdf", scaled_df, ylim=c(0.0,3))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-30.pdf", scaled_df, ylim=c(0.0,30))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-50.pdf", scaled_df, ylim=c(0.0,50))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-100.pdf", scaled_df, ylim=c(0.0,100))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-150.pdf", scaled_df, ylim=c(0.0,150))
#lattice_Scaledmeta.proSeq(filename="Jurkat.nc1.LenNorm-0-250.pdf", scaled_df, ylim=c(0.0,150))



lattice_Scaledmeta.2group = function(filename= "test.pdf", df, df.col=c(rgb(0,0.5,0.5), rgb(0,0,0.5)), ylim = c(0, 10), main = NULL, 
                                     xlab = "Distance to gene boundaries (bp)",  ylab = "Median PRO-seq intensity")
{
  pdf(paste(fig_dir, filename), width = 35, height = 10)
  result <- xyplot(mean ~ x, data = df,
                  group = factor(sample, labels = c("Experiment", "Prediction")), 
                  scales = list(tck=c(1,0),alternating = c(1,1),
                                 x=list(relation='free',axs='i', labels = c('-5000', 'TSS', "+2000", "-2000", "CPS", "30000"), at = c(1, 601, 641, 1641, 1681, 2280)),
#                                y=list(log = 10, equispaced.log = FALSE)), ##Logscale
                                 y=list(relation='free',axs='i')), ##Nonlog

                   key=list( corner=c(0.98,0.85),
                             padding.text=3,
                            text=list(c("Experiment", "Prediction")),col=rep("black", 9), cex=0.8, font=2),
                            rectangles=list(col=df.col),
                   type = 'l', 
                   ylim = ylim,
                   xlim = c(500, 2280),
                   col = df.col, 
                   fill = df.col, 
                   ylab = ylab,
                   xlab = xlab,
                   main = main,
                   aspect=0.5,
                   lwd=3,
                   par.settings = list(strip.background=list(col="lightgrey"),par.xlab.text=list(cex=1.1,font=2),par.ylab.text=list(cex=1.1,font=2),axis.text=list(cex=1)),
                   par.strip.text=list(cex=1.2,font=2),
                   upper = df$upper, 
                   lower = df$lower, 
                   panel = function(x, y, ...){
                     panel.grid(h=-1, v=-1)
                     panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     panel.xyplot(x, y, ...)
                   }) + latticeExtra::layer(panel.abline(v = c(601,641,1641,1681), lty = 2, lwd = 2, col = 'gray'))
  print(result)
  dev.off()
  print(result)
}

lattice_Scaledmeta.2group(filename="Metaplot.H3k27ac.150.pdf",  scaled_df_H3k27ac,  ylab ="H3k27ac",  ylim=c(0.0,150))
lattice_Scaledmeta.2group(filename="Metaplot.H3k27me3.150.pdf", scaled_df_H3k27me3, ylab ="H3k27me3", ylim=c(0.0,150))
lattice_Scaledmeta.2group(filename="Metaplot.H3k36me3.150.pdf", scaled_df_H3k36me3, ylab ="H3k36me3", ylim=c(0.0,150))
lattice_Scaledmeta.2group(filename="Metaplot.H3k4me1.150.pdf",  scaled_df_H3k4me1,  ylab ="H3k4me1",  ylim=c(0.0,150))
lattice_Scaledmeta.2group(filename="Metaplot.H3k4me2.150.pdf",  scaled_df_H3k4me2,  ylab ="H3k4me2",  ylim=c(0.0,150))
lattice_Scaledmeta.2group(filename="Metaplot.H3k4me3.150.pdf",  scaled_df_H3k4me3,  ylab ="H3k4me3",  ylim=c(0.0,100))
lattice_Scaledmeta.2group(filename="Metaplot.H3k9ac.150.pdf",   scaled_df_H3k9ac,   ylab ="H3k9ac",   ylim=c(0.0,150))
lattice_Scaledmeta.2group(filename="Metaplot.H3k9me3.150.pdf",  scaled_df_H3k9me3,  ylab ="H3k9me3",  ylim=c(0.0,150))
lattice_Scaledmeta.2group(filename="Metaplot.H4k20me1.150.pdf", scaled_df_H4k20me1, ylab ="H4k20me1", ylim=c(0.0,150))




