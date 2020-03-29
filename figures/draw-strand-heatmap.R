source("/local/workdir/zw355/proj/prj10-dreg/paper-fig-201808/heatmaps.R");
source("/local/workdir/zw355/proj/prj15-histone/scripts/hist.param.R");

read_read_mat <- function(file.bw.pred, dregX, dist, step, navg = 20, times=1 )
{
    hmark.pred <- load.bigWig( file.bw.pred )
    ## Get a matrix of counts.
    hCountMatrix <- bed.step.bpQuery.bigWig(hmark.pred, center.bed(dregX[,c(1,2,3)], dist, dist), step=step, abs.value=TRUE)
    hmat <- log( times * matrix(unlist(hCountMatrix), nrow= NROW(dregX), byrow=TRUE)+1) ;

    #hm_order <- order( dregX$score, decreasing=TRUE);
    #if(is.null(hm_order)) {
    #  hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -10):(NCOL(hmat)/2 +10)]), decreasing=TRUE)
      #  hmat <- hmat[hm_order,]
    #}

    ## Average by rows of 20.
    ## navg : Average every navg rows
    avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
    hmat.pred <- avgMat
}

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


read_read_mat <-function (file.bw.org, dregX, dist, step, navg = 20, times=1)
{
    hmark <- load.bigWig( file.bw.org )
    hCountMatrix <- bed.step.bpQuery.bigWig( hmark, center.bed(dregX[,c(1,2,3)], dist, dist), step=step, abs.value=TRUE)
    
    hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(dregX[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
    hmat <- log( times * matrix(unlist(hCountMatrix), nrow= NROW(dregX), byrow=TRUE)+1) ;
    avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
    
    unload.bigWig(hmark);
    
    return(avgMat);    
}

heatmap<-function( df.bed.strand, file.plus.bw, file.minus.bw, file.bw.org, file.peak.org, file.bw.pred, file.peak.pred, file.png, 
                   subs = NULL, breaks = NULL, cols = NULL, dist = 1000, step = 25) 
{
   if(is.na(file.bw.org) || is.na(file.bw.pred))
      return;

    dregX <- df.bed.strand 
    dregX <- dregX[grep("_|chrY|chrX", dregX[,1], invert=TRUE),]
    dregX <- dregX[ dregX[,2]> dist,]

    hmat.org <- read_read_mat (file.bw.org, dregX, dist, step, times=1)
    hmat.pred <- read_read_mat (file.bw.pred, dregX, dist, step, times=10)

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

    ##part 2
    gt <- pheatmap( hmat.pred, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.pred, breaks = bk.pred, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3 );
    grid.draw(gt$gtable)
    popViewport()

    ##part 1:
    gt <- pheatmap( hmat.org, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.org, breaks = bk.org, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
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


#K562 

file.TSS.bed <- "/workdir/zw355/proj/prj15-histone/pred-k562/K562_DivergentPairs.bed"
file.plus.bw  <- "/local/workdir/zw355/proj/prj10-dreg/k562/K562_unt.sort.bed.gz_plus.bw"
file.minus.bw <- "/local/workdir/zw355/proj/prj10-dreg/k562/K562_unt.sort.bed.gz_minus.bw"

tb <- read.table(file.TSS.bed, skip=1)
bed.Tss <- bed.region.temp<-c();

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
        r0 <- tb[bed.region.temp*2+c(-1:0),];
        r0.max <- r0[ which.max(r0[,5]) , ,drop=F];
        r0.max.center <- (r0.max [1, 3] + r0.max [1, 2])/2

        r0 <- r0[ r0[,6] != r0.max[1,6],, drop=F];
        r0.min <- r0[ which.min(r0[,5]), ,drop=F];
        r0.min.center <- (r0.min[1, 2]+r0.min[1, 3])/2

        r0.dist <- r0.max.center - r0.min.center;
        
        r0 <- rbind(r0.max, r0.min);
        r0.plus <- r0[ r0[,6]=="+", ,drop=F];
        r0.minus <- r0[ r0[,6]=="-", ,drop=F];
        
        if(r0.plus[1,3] > r0.minus[1,3])
           r.type = 1 #divergent
        else
           r.type = -1  #divergent 
        
        gene <- if(!is.na(r0.max[1,4])) r0.max[1,4] else r0.min[1,4];
    
        bed.Tss<-rbind(bed.Tss, data.frame(chr=r0.max[1,1], start=r0.max.center-1000, 
              end=r0.max.center+1000, 
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

if(0)
{
file.png <- "Alex-H3K27ac-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.Alex.Mnase.H3k27ac.bw, file.Alex.Mnase.H3k27ac.peak, 
        "../pred-alex/raw.Alex-Mnase-H3k27ac.pred.TSS.bw", file.Alex.Mnase.H3k27ac.peak, 
        file.png);

file.png <- "Alex-H3K4me1-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.Alex.Mnase.H3k4me1.bw, file.Alex.Mnase.H3k4me1.peak, 
        "../pred-alex/raw.Alex-Mnase-H3k4me1.pred.TSS.bw", file.Alex.Mnase.H3k4me1.peak, 
        file.png);

file.png <- "Alex-H3K4me2-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.Alex.Mnase.H3k4me2.bw, file.Alex.Mnase.H3k4me2.peak, 
        "../pred-alex/raw.Alex-Mnase-H3k4me2.pred.TSS.bw", file.Alex.Mnase.H3k4me2.peak, 
        file.png);

file.png <- "Alex-H3K4me3-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.Alex.Mnase.H3k4me3.bw, file.Alex.Mnase.H3k4me3.peak, 
        "../pred-alex/raw.Alex-Mnase-H3k4me3.pred.TSS.bw", file.Alex.Mnase.H3k4me3.peak, 
        file.png);
}

if(0)
{
file.png <- "GM-H3K27ac-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.gm.H3k27ac.bw, file.gm.H3k27ac.peak, 
        file.gm.H3k27ac.pred.raw, file.gm.H3k27ac.peak, 
        file.png);

file.png <- "GM-H3K4me1-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.gm.H3k4me1.bw, file.gm.H3k4me1.peak, 
        file.gm.H3k4me1.pred.raw, file.gm.H3k4me1.peak, 
        file.png);

file.png <- "GM-H3K4me2-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.gm.H3k4me2.bw, file.GM.H3k4me2.peak, 
        file.gm.H3k4me2.pred.raw, file.GM.H3k4me2.peak, 
        file.png);

file.png <- "GM-H3K4me3-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.gm.H3k4me3.bw, file.gm.H3k4me3.peak, 
        file.gm.H3k4me3.pred.raw, file.gm.H3k4me3.peak, 
        file.png);
}

if(1)
{
    file.png <- "K562-H3K4me1-strand-heatmap.png"
    heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k4me1.bw, file.k562.H3k4me1.peak, 
        file.k562.H3k4me1.pred.raw, file.k562.H3k4me1.peak, 
        file.png);
        
    file.png <- "K562-H3K4me2-strand-heatmap.png"
    heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k4me2.bw, file.k562.H3k4me2.peak, 
        file.k562.H3k4me2.pred.raw, file.k562.H3k4me2.peak, 
        file.png);

    file.png <- "K562-H3K4me3-strand-heatmap.png"
    heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k4me3.bw, file.k562.H3k4me3.peak, 
        file.k562.H3k4me3.pred.raw, file.k562.H3k4me3.peak, 
        file.png);

    file.png <- "K562-H3K27ac-strand-heatmap.png"
    heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k27ac.bw, file.k562.H3k27ac.peak, 
        file.k562.H3k27ac.pred.raw, file.k562.H3k27ac.peak, 
        file.png);
}

        
if(0)
{
file.png <- "K562-H3K122ac-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k122ac.bw, file.k562.H3k122ac.peak, 
        file.k562.H3k122ac.pred.bw, file.k562.H3k122ac.peak, 
        file.png);


file.png <- "K562-H3K27me3-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k27me3.bw, file.k562.H3k27me3.peak, 
        file.k562.H3k27me3.pred.bw, file.k562.H3k27me3.peak, 
        file.png);


file.png <- "K562-H3K36me3-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k36me3.bw, file.k562.H3k36me3.peak, 
        file.k562.H3k36me3.pred.bw, file.k562.H3k36me3.peak, 
        file.png);


file.png <- "K562-H3K9ac-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k9ac.bw, file.k562.H3k9ac.peak, 
        file.k562.H3k9ac.pred.bw, file.k562.H3k9ac.peak, 
        file.png);
 
file.png <- "K562-H3K9me3-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H3k9me3.bw, file.k562.H3k9me3.peak, 
        file.k562.H3k9me3.pred.bw, file.k562.H3k9me3.peak, 
        file.png);


file.png <- "K562-H4K20me1-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.k562.H4k20me1.bw, file.k562.H4k20me1.peak, 
        file.k562.H4k20me1.pred.bw, file.k562.H4k20me1.peak, 
        file.png);

file.png <- "Alex-H3K27ac-strand-heatmap.png"
heatmap(bed.Tss, file.plus.bw, file.minus.bw, 
        file.Alex.Mnase.H3k27ac.bw, file.Alex.Mnase.H3k27ac.peak, 
        file.Alex.Mnase.H3k27ac.bw, file.Alex.Mnase.H3k27ac.peak, 
        file.png);

}
