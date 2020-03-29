##
## Creates heatmaps on human dREG-HD sites.

require(bigWig)
require(pheatmap)
require(RColorBrewer)
library(gtable)
library(grid)
library(sqldf);

options(scipen =99)

rowMax <- function(x) { sapply(1:NROW(x), function(i) {return(max(x[i,], na.rm=TRUE))}) }

transform.bed <- function(anchor, bed, upstreamWindow, downstreamWindow) {
  start = anchor - upstreamWindow
  end = anchor + downstreamWindow + 1

  N = dim(bed)[2]
  if (N >= 6) { # use strands
    is_minus = bed[,6] == '-'
    start[is_minus] = anchor[is_minus] - downstreamWindow
    end[is_minus] = anchor[is_minus] + upstreamWindow + 1
  }
  start = as.integer(start)
  end = as.integer(end)

  res = NULL
  if (N > 3)
    res = data.frame(bed[,1], start, end, bed[,4:N])
  else
    res = data.frame(bed[,1], start, end)

  colnames(res) <- colnames(bed)
  rownames(res) <- rownames(bed)

  return(res)
}

center.bed <- function(bed, upstreamWindow, downstreamWindow) {
  anchor = bed[,2] + ((bed[,3] - bed[,2]) %/% 2)

  transform.bed(anchor, bed, upstreamWindow, downstreamWindow)
}


dist_X_center <- function( pred.bed, file.X )
{
    options("scipen"=100, "digits"=4);

    tmp.bed <- tempfile(fileext=".bed" )
    write.table( pred.bed[,c(1:3)], file=tmp.bed, row.names=F, col.names=F, quote=F, sep="\t");
    system( paste("sort-bed ",  tmp.bed, " > ", tmp.bed, ".sorted", sep="") );

    tb.X <- read.table( file.X, header=F);
    tb.X[,2] <- round((tb.X[,2] + tb.X[,3])/2)
    tb.X[,3] <- tb.X[,2] + 1;
    tmp2.bed <- tempfile(fileext=".bed")
    write.table(tb.X[, c(1:3)], file=tmp2.bed, row.names=F, col.names=F, quote=F, sep="\t");

    tb.close <-  read.table( file = pipe(paste("sort-bed ",  tmp2.bed, " | bedtools closest -d -a ", tmp.bed, ".sorted -b - -t first", sep="") ) );
    return(tb.close);
}

dist_X_close <- function( pred.bed, file.X )
{
    options("scipen"=100, "digits"=4);

    tmp.bed <- tempfile(fileext=".bed" )
    write.table( pred.bed[,c(1:3)], file=tmp.bed, row.names=F, col.names=F, quote=F, sep="\t");
    system( paste("sort-bed ",  tmp.bed, " > ", tmp.bed, ".sorted", sep="") );

    tb.X <- read.table( file.X, header=F);
    tb.X0 <- tb.X[tb.X$V5=="-",];
    tb.X1 <- tb.X[tb.X$V5=="+",];
    tb.X0$V2 <- tb.X0$V3 - 1;
    tb.X1$V3 <- tb.X1$V2 + 1;

    tmp2.bed <- tempfile(fileext=".bed")
    write.table( rbind(tb.X0[, c(1:3)], tb.X1[, c(1:3)]),  file=tmp2.bed, row.names=F, col.names=F, quote=F, sep="\t");
    tb.close <-  read.table( file = pipe(paste("sort-bed ",  tmp2.bed, " | bedtools closest -d -a ", tmp.bed, ".sorted -b - -t first", sep="") ) );

    return(tb.close);
}

draw_legend <- function(bk, hmcols)
{
    par(mar=c(10,0,0,0), plt=c(0.15, 1.00, 0.4, 0.8 ));
    plot(NA,NA, type="n", xlim=c(1, NROW(bk)*0.85/0.7), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n",bty="n" );
    for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
    axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),0), cex.axis=3, cex=3, tick=FALSE );
}

pheatmap_2layer <- function( hmat.plus, hmat.minus, breaks.plus=NULL, breaks.minus=NULL, cluster_rows = FALSE, cluster_cols = FALSE, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
{
    get_color_map <- function(pMat, mMat)
    {
        cMap <- matrix("", nrow=NROW(pMat), ncol=NCOL(pMat));

        if(is.null(breaks.minus)) breaks.minus <- hmat.minus;
        mMat <-  ( mMat -  min(breaks.minus)) / (max(breaks.minus) -  min(breaks.minus));
        if ( sum(mMat<0)>0 ) mMat[mMat<0] <- 0;
        if ( sum(mMat>1)>0 ) mMat[mMat>1] <- 1;

        if(is.null(breaks.plus)) breaks.plus <- hmat.plus;
        pMat <-  ( pMat -  min(breaks.plus)) / (max(breaks.plus) -  min(breaks.plus));
        if ( sum(pMat<0)>0 ) pMat[pMat<0] <- 0;
        if ( sum(pMat>1)>0 ) pMat[pMat>1] <- 1;

        for(i in 1:NROW(pMat))
        for(j in 1:NCOL(pMat))
        {
            B <- max( mMat[i,j], 1-pMat[i,j] );
            G <- 1 - max( mMat[i,j], pMat[i,j] );
            R <- max( 1- mMat[i,j], pMat[i,j] );

            cMap[i,j] <- rgb(R, G, B);
        }
        return(cMap);
    }

    gt <- pheatmap( hmat.plus, cluster_rows = cluster_rows, cluster_cols = cluster_cols, col= rgb( 0, 0, 1 - breaks.plus/max(breaks.plus)), breaks = breaks.plus, legend=legend, show_rownames=show_rownames, show_colnames=show_colnames, silent=silent );
    gt$gtable$grobs[[1]]$children[[1]]$gp$fill <- get_color_map( hmat.plus, hmat.minus);

    return(gt);
}



writeProseqHeatmap<- function(bed, file.plus.bw, file.minus.bw, png.name, subs= NULL, breaks= NULL,
                            cols= NULL, dist= 1000, step=25)
{
    navg <- 20 ## Average every navg rows


    ## Load GROseq data .
    hPlus <- load.bigWig(file.plus.bw)
    hMinus <- load.bigWig(file.minus.bw)

    ## Get a matrix of counts.
    hPlusMatrix <- bed.step.bpQuery.bigWig(hPlus, center.bed(bed[,c(1,2,2)], dist, dist), step=step, abs.value=TRUE)
    hPlusMatrix.rev <- lapply(1:NROW(bed), function(i){ rev(hPlusMatrix[[i]]) })
    hPlus.log <- log(do.call("rbind", hPlusMatrix )+1);
    hPlus.rev.log <- log(do.call("rbind", hPlusMatrix.rev )+1);

    hMinusMatrix <- bed.step.bpQuery.bigWig(hMinus, center.bed(bed[,c(1,2,2)], dist, dist), step=step, abs.value=TRUE)
    hMinusMatrix.rev <- lapply(1:NROW(bed), function(i){ rev(hMinusMatrix [[i]]) })
    hMinus.log <- log(do.call("rbind", hMinusMatrix)+1);
    hMinus.rev.log <- log(do.call("rbind", hMinusMatrix.rev )+1);

    #            + TSS           - TSS
    # plus       plus bigWig      reverse(minus bigWig)       
    # minus      minus bigwig     reverse(plus bigWig)       
    
    hmat.plus  <- do.call("rbind", lapply(1:NROW(bed), function(i){if(bed[i,6]=="-") hMinus.rev.log[i,] else hPlus.log[i,] }))
    hmat.minus <- do.call("rbind", lapply(1:NROW(bed), function(i){if(bed[i,6]=="-") hPlus.rev.log[i,] else hMinus.log[i,]})) 

    ## Average by rows of 10.
    navg <- 20 ## Average every navg rows
    avgMat <- t(sapply(1:floor(NROW(hmat.plus)/navg), function(x) {colMeans(hmat.plus[((x-1)*navg+1):min(NROW(hmat.plus),(x*navg)),])}))
    hmat.plus <- avgMat

    navg <- 20 ## Average every navg rows
    avgMat <- t(sapply(1:floor(NROW(hmat.minus)/navg), function(x) {colMeans(hmat.minus[((x-1)*navg+1):min(NROW(hmat.minus),(x*navg)),])}))
    hmat.minus <- avgMat

    ## Write out a heatmap.
    if(is.null(breaks)) {
        bk <- seq(min(hmat.plus), max(hmat.plus)+0.01, 0.01)
    } else {
        bk <- breaks
    }

    if(is.null(cols)) {
        hmcols.plus <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
    } else {
        hmcols.plus <- colorRampPalette(cols)(length(bk)-1) # red
    }
    breaks.plus <- bk;

    if(is.null(breaks)) {
        bk <- seq(min(hmat.minus), max(hmat.minus)+0.01, 0.01)
    } else {
        bk <- breaks
    }

    if(is.null(cols)) {
        hmcols.minus <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
    } else {
        hmcols.minus <- colorRampPalette(cols)(length(bk)-1) # red
    }
    breaks.minus <- bk;

    png(paste(png.name,".trans.png",sep=""), width=450, height = 800*1.55 )

    n.all <- NROW(hmat.plus) ;
    lay.heights <- c(0.85, 0.15);
    lay.widths  <- c(0.5, 0.5)
    layout(matrix(c(1, 1, 2, 3 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)

    ##part1:
    gt <- pheatmap_2layer( hmat.plus, hmat.minus, breaks.plus, breaks.minus );

    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = c(1,2))
    grid.draw(gt$gtable)
    popViewport()

    ##part 2:
    ## draw colorScale for Heatmap
    bk <- breaks.minus;
    hmcols <- rgb(1, 1-breaks.minus/max(breaks.minus), 1-breaks.minus/max(breaks.minus) );
    {
        par(mar=c(5,0,0,0), plt=c(0.15, 0.85, 0.4, 0.8 ));
        plot(NA,NA, type="n", xlim=c(1, NROW(bk)), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n" );
        for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
        axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),1), cex.axis=3, cex=3, tick=FALSE );
    }

    ##part 4:
    ## draw colorScale for peak track
    bk <- breaks.plus;
    hmcols <- rgb(1-breaks.plus/max(breaks.plus), 1-breaks.plus/max(breaks.plus), 1 );
    {
        par(mar=c(5,0,0,0), plt=c(0.15, 0.85, 0.4, 0.8 ));
        plot(NA,NA, type="n", xlim=c(1, NROW(bk)), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n" );
        for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
        axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),1), cex.axis=3, cex=3, tick=FALSE );
    }
    dev.off();


    ## start new PNG file
    png(paste(png.name,".2blk.png",sep=""), width=450*1.5, height = 800*1.2 )

    lay.heights <- c(0.85, 0.15);
    lay.widths  <- c(0.48, 0.04, 0.48 )
    layout(matrix(c(1, 5, 2, 3, 5, 4 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)

    ##part 1
    gt <- pheatmap( hmat.plus, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.plus, breaks = breaks.plus, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1 );
    grid.draw(gt$gtable)
    popViewport()

    ##part 2:
    gt <- pheatmap( hmat.minus, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.minus, breaks = breaks.minus, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3)
    grid.draw(gt$gtable)
    popViewport()

    ##part 3:
    ## draw colorScale for peak track
    draw_legend(breaks.plus, hmcols.plus );
    
    ##part 4:
    ## draw colorScale for Heatmap
    draw_legend(breaks.minus, hmcols.minus );

    dev.off(); 

    return()
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
    pheatmap(rev(bk), cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=TRUE, legend_breaks= quantile(bk), legend_labels= signif(exp(quantile(bk)),3), show_rownames=FALSE, show_colnames=FALSE)
    dev.off()
}

if(1)
{
   file.TSS.bed <- "/workdir/zw355/proj/prj15-histone/pred-gm/hg19.gm12878.TIDs.bed"
   file.plus.bw  <- "/local/workdir/zw355/proj/prj10-dreg/GM12878/groseq_plus.bigWig"
   file.minus.bw <- "/local/workdir/zw355/proj/prj10-dreg/GM12878/groseq_minus.bigWig"

   tb <- read.table(file.TSS.bed)
   bed.Tss <- bed.region.temp<-c();
   win.size <- 0

   for(i in 1:((NROW(tb)/2-1)))
   {
       bed.region.temp <- c(bed.region.temp, i);

       r0 <- tb[i*2+c(-1:0),];
       r0.max <- r0[ which.max(r0[,5]) , ,drop=F];
       r0.min <- r0[ which.min(r0[,5]) , ,drop=F];
       r1 <- tb[i*2+c(1:2),];
       r1.max <- r1[ which.max(r1[,5]) , ,drop=F];
       r1.min <- r1[ which.min(r1[,5]) , ,drop=F];

       gap <- max(r0.max[,2:3], r0.min[,2:3] ) - min(r1.max[,2:3], r1.min[,2:3] );
       if(abs(gap)<=600)
       {   
           cat(i, gap, "\n");
           next;
       }
       else
       {
           #find the adjcent Tss region which distantce is grant than 600 
           r0 <- tb[bed.region.temp*2+c(-1:0),];
           
           #find the highest score as the maximum center
           r0.max <- r0[ which.max(r0[,5]) , ,drop=F];
           r0.max.center <- (r0.max [1, 3] + r0.max [1, 2])/2
           
           #find the lowest score  which is not same strand
           r0 <- r0[ r0[,6] != r0.max[1,6],, drop=F];
           r0.min <- r0[ which.min(r0[,5]), ,drop=F];
           r0.min.center <- (r0.min[1, 2]+r0.min[1, 3])/2
           
           # distance bw two points
           r0.dist <- r0.max.center - r0.min.center;

           # select the plus and minus strand
           r0 <- rbind(r0.max, r0.min);
           r0.plus <- r0[ r0[,6]=="+", ,drop=F];
           r0.minus <- r0[ r0[,6]=="-", ,drop=F];
        
           # if plus strand is great than minus strand ==> divergent=1 
           #if(r0.plus[1,3] > r0.minus[1,3])
           #   r.type = 1 #divergent
           #else
           #   #r.type = -1  #divergent 
        
           r.type = 0
           # if plus strand is great than minus strand ==> divergent=1 
           if(r0.plus[1,2] >= r0.minus[1,3]) r.type = 1 #divergent 
           if(r0.plus[1,3] <= r0.minus[1,2]) r.type = -1

           gene <- if(!is.na(r0.max[1,4])) r0.max[1,4] else r0.min[1,4];
    
           bed.Tss<-rbind(bed.Tss, data.frame(chr=r0.max[1,1], start=r0.max.center-0, 
               end=r0.max.center+0, 
               dist=r0.dist, 
               score=r0.max[1,5], 
               strand=r0.max[1,6],
               type=r.type,
               gene=gene ));
              
         bed.region.temp <-c();
      }
   }

   bed.Tss$dist <- abs(bed.Tss$dist)
   bed.Tss <- bed.Tss[order(bed.Tss$dist),]

   # only write heatmap for divergent regions
   sup <- writeProseqHeatmap( bed.Tss[bed.Tss$type != -1,], file.plus.bw, file.minus.bw, "GM-TID")

}

if(0)
{
    #MNase.path="/fs/cbsudanko/storage/data/hg19/k562/sydh_mnase/"
    #file.MNase.bw <- "wgEncodeSydhNsomeK562Sig.bigWig"

    file.dnase.bw = "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeOpenChromDnaseK562SigV2.bigWig";
    file.dnase.peak = "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz";
    file.protein.coding <- "./gencode.v19.protein.coding.bed"

    sup <- writeProseqHeatmap( bed.Tss, file.plus.bw, file.minus.bw, "GM-TID")
}
