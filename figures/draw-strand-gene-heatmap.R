source("/local/workdir/zw355/proj/prj10-dreg/paper-fig-201808/heatmaps.R");
source("../scripts/hist.param.R");
source("../scripts/hist.svm.com.R");


read_peak_overlap <- function(dregX, file.peak, navg = 20)
{
    file.tmp <- tempfile(fileext=".bed");
    write.table( dregX[,c(1:3)], file=file.tmp, quote=F, row.names=F, col.names=F, sep="\t");
    tb <- unique(read.table(pipe(paste("zcat ", file.peak, "| awk -v OFS='\t' '{print $1,$2,$3 }' - | bedtools intersect -a" , file.tmp, "-b - -loj")))[,c(1:6)]);
    tbovp <- unlist(lapply(1:NROW(tb), function(x){
        if(as.character(tb[x,4]) == ".") return(0)
        else
            return(1);
    }));

    tborg <- dregX[,c(1:3)];
    colnames(tborg) <- c("chr", "start", "end");
    tbovp <- data.frame(tb[,c(1:3)], ratio=tbovp);
    colnames(tbovp) <- c("chr", "start", "end", "ratio");
    tbovp <- sqldf("select chr, start, end, max(ratio) as ratio from tbovp group by chr, start, end")
    tbovp <- sqldf("select tborg.chr, tborg.start, tborg.end, tbovp.ratio from tborg left join tbovp on tborg.chr=tbovp.chr and tborg.start=tbovp.start");

    show(all(tbovp[,c(1:3)]==data.frame(dregX[,c(1:3)], stringsAsFactors=F)))
    ovp <- sapply(1:floor(NROW(tbovp)/navg), function(x) {mean(tbovp[((x-1)*navg+1):min(NROW(tbovp),(x*navg)),4])})

    ovp.min <- quantile( ovp, probs=0.05);
    ovp.max <- quantile( ovp, probs=0.95);
    
    ovp.norm <- (ovp - ovp.min)/(ovp.max- ovp.min);
    if (NROW(which(ovp.norm<0))>0) ovp.norm [ which(ovp.norm<0) ] <- 0;
    if (NROW(which(ovp.norm>1))>0) ovp.norm [ which(ovp.norm>1) ] <- 1;
    
    return(ovp.norm);
}

draw_legend <- function(bk, hmcols)
{
    par(mar=c(10,0,0,0), plt=c(0.15, 1.00, 0.4, 0.8 ));
    plot(NA,NA, type="n", xlim=c(1, NROW(bk)*0.85/0.7), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n",bty="n" );
    for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
    axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),0), cex.axis=3, cex=3, tick=FALSE );
}


read_read_mat0 <-function (file.bw.org, dregX, step, navg = 20, times=1)
{
    hmark <- load.bigWig( file.bw.org )

    hCountMatrix <- bed.step.bpQuery.bigWig( hmark, dregX[,c(1,2,3)], step=step, abs.value=TRUE)
    hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(dregX[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
    hmat <- log( times * matrix(unlist(hCountMatrix), nrow= NROW(dregX), byrow=TRUE)+1) ;

    avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
    
    unload.bigWig(hmark);
    
    return(avgMat);    
}

heatmap.gene<-function( df.bed.strand, file.plus.bw, file.minus.bw, file.bw.org, file.peak.org, file.bw.pred, file.peak.pred, file.png, 
                   subs = NULL, breaks = NULL, cols = NULL, step = 25) 
{
   if(is.na(file.bw.org) || is.na(file.bw.pred))
      return;

    dregX <- df.bed.strand 
    dregX <- dregX[grep("_|chrY|chrX", dregX[,1], invert=TRUE),]

    hmat.pred <- read_read_mat0 (file.bw.pred, dregX[,c(1:6)], step, times=10)
    hmat.org  <- read_read_mat0 (file.bw.org, dregX[,c(1:6)],  step, times=1)

    ## Write out a heatmap.
    if(is.null(breaks)) {
        bk.org <- seq(min(hmat.org), max(hmat.org), 0.01)
    } else {
        bk.org <- breaks
    }

    if(is.null(cols)) {
        hmcols.org <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk.org)-1))
    } else {
        hmcols.org <- colorRampPalette(cols)(length(bk.org)-1) # red
    }

    if(is.null(breaks)) {
        bk.pred <- seq(min(hmat.pred), max(hmat.pred), 0.01)
    } else {
        bk.pred <- breaks
    }

    if(is.null(cols)) {
        hmcols.pred <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk.pred)-1))
    } else {
        hmcols.pred <- colorRampPalette(cols)(length(bk.pred)-1) # red
    }

    ## start new PNG file
    png(file.png, width=450*1.5, height = 800*1.2 )

    lay.heights <- c(0.85, 0.15);
    lay.widths  <- c(0.48, 0.04, 0.48 )
    layout(matrix(c(1, 5, 2, 3, 5, 4 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)


    ##part 1:
    gt <- pheatmap( hmat.org, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.org, breaks = bk.org, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
    #gt$gtable$grobs[[1]]$children[[1]]$gp$fill[,20:21] <- "#000000";
    gt$gtable$grobs[[1]]$children[[1]]$gp$fill[,20] <- "#000000";
    grid.draw(gt$gtable)
    popViewport()

    ##part 2
    gt <- pheatmap( hmat.pred, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.pred, breaks = bk.pred, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3 );
    #gt$gtable$grobs[[1]]$children[[1]]$gp$fill[,20:21] <- "#000000";
    gt$gtable$grobs[[1]]$children[[1]]$gp$fill[,20] <- "#000000";
    grid.draw(gt$gtable)
    popViewport()


    ##part 3:
    ## draw colorScale for peak track
    draw_legend(bk.org, hmcols.org );
    
    ##part 4:
    ## draw colorScale for Heatmap
    draw_legend(bk.pred, hmcols.pred );

    dev.off(); 

    invisible()
  
} 

file.TSS.bed <- "/workdir/zw355/proj/prj15-histone/pred-k562/K562_DivergentPairs.bed"
file.plus.bw  <- "/local/workdir/zw355/proj/prj10-dreg/k562/K562_unt.sort.bed.gz_plus.bw"
file.minus.bw <- "/local/workdir/zw355/proj/prj10-dreg/k562/K562_unt.sort.bed.gz_minus.bw"
file.gene <- "./hg19.gene.table.txt"

bed.Tss <- read.table(file.TSS.bed, skip=1)
bed.Tss <- bed.Tss[!is.na(bed.Tss[,4]),]

bed.Tss1 <- bed.Tss[bed.Tss$V6=="+",]
bed.Tss1[,2] <- bed.Tss1[,3] - 2000
bed.Tss1[,3] <- bed.Tss1[,3] + 50000

bed.Tss0 <- bed.Tss[bed.Tss$V6=="-",]
bed.Tss0[,3] <- bed.Tss0[,2] + 2000
bed.Tss0[,2] <- bed.Tss0[,2] - 50000

tb0 <- rbind(bed.Tss1, bed.Tss0);
tbg <- read.table(file.gene, sep=" ");
tbg$length <- tbg[,4]-tbg[,3] 
library(sqldf)
tb0 <- sqldf("select tb0.V1, tb0.V2, tb0.V3, tb0.V4, tb0.V5, tb0.V6, tbg.length from tb0 left join tbg on tb0.V4=tbg.V1 order by tbg.length ");
tb0 <- unique(tb0[!is.na(tb0$length),])

if(0)
{
write.bed(tb0, file="Alex.Mnase.Tss.sup.bed", compress=FALSE)

file.png <- "Alex-H3K36me3-strand-heatmap.png"
heatmap.gene(tb0, file.plus.bw, file.minus.bw, 
        file.Alex.Mnase.H3k36me3.bw, file.Alex.Mnase.H3k36me3.peak, 
        "../pred-alex/raw.Alex-Mnase-H3k36me3.pred.TSS.sup.bw", file.Alex.Mnase.H3k36me3.peak, 
        file.png, step = 100);

file.png <- "Alex-H3K79me3-strand-heatmap.png"
heatmap.gene(tb0, file.plus.bw, file.minus.bw, 
        file.Alex.Mnase.H3k79me3.bw, file.Alex.Mnase.H3k79me3.peak, 
        "../pred-alex/raw.Alex-Mnase-H3k79me3.pred.TSS.sup.bw", file.Alex.Mnase.H3k79me3.peak, 
        file.png, step = 100);
}

if(1)
{
file.png <- "GM-H3K36me3-strand-heatmap.png"
heatmap.gene(tb0, file.plus.bw, file.minus.bw, 
        file.gm.H3k36me3.bw, file.gm.H3k36me3.peak, 
        file.gm.H3k36me3.pred.raw, file.gm.H3k36me3.peak, 
        file.png, step = 100);

#file.png <- "GM-H3K79me3-strand-heatmap.png"
#heatmap.gene(tb0, file.plus.bw, file.minus.bw,  
#        file.gm.H3k79me3.bw, file.gm.H3k79me3.peak, 
#        file.gm.H3k79me3.pred.raw, file.gm.H3k79me3.peak, 
#        file.png, step = 100);
}


