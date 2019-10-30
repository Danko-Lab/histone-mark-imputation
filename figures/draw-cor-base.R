
library(sqldf);
library(bigWig)
library(dREG)
library(data.table);
library(snowfall);
library(minerva);
library(parallel);
library(philentropy)

source("../scripts/hist.svm.com.R");
source("../scripts/hist.param.R");

calculate_cor <- function (y.exp, y.pred , win.size=2000 )
{
   if(win.size!=10)
   {
       n.per.range <- round( win.size/10 );

       y.range.exp <-  unlist(lapply(1:floor(NROW(y.exp)/n.per.range),  function(i) {sum(y.exp [(i-1)*n.per.range+c(1:n.per.range)])}));
       y.range.pred <- unlist(lapply(1:floor(NROW(y.pred)/n.per.range), function(i) {sum(y.pred[(i-1)*n.per.range+c(1:n.per.range)])}));
   }
   else
   {
        y.range.exp <- y.exp;
        y.range.pred <- y.pred;
   }

   return(list(r.pearson = cor( y.range.exp, y.range.pred, method = "pearson" ),
      r.spearman = cor( y.range.exp, y.range.pred, method = "spearman" ),
      r.mad = mad( y.range.exp -  y.range.pred ) )  );
}

get_histone_peak_max<-function( file.peak, file.histone)
{
   bw.hist <- load.bigWig(file.histone);
   bedA <- read.table(file.peak)[,c(1:3)];

   y <- try( unlist(bed.region.bpQuery.bigWig( bw=bw.hist, bed= bedA, op="max") ) );
   if(class(y)=="try-error")
   {
       stop("error");
     return(NULL);
   }
   unload.bigWig(bw.hist);
   return(data.frame(bedA, y));
}

get_histone_readX<-function(bedA, file.histone, ref.bed=NULL, method="mean", normalized=TRUE, ncores=1)
{
   bedA <- bedA[,c(1:3)];
   na.idx <- which( is.na(bedA[,1]) |  is.na(bedA[,2]) |  is.na(bedA[,3])  )
   if(length(na.idx)>0)
         stop("NA in bedA\n");

   ## This part occupies too many memory by the mclapply threads
   ## ABANDON
   if(0)
   {
      y <- unlist( mclapply(1:ceiling(NROW(bedA)/1000000), function(i) {
         bw.hist <- load.bigWig(file.histone);
         start <- 1+(i-1)*1000000;
         end <- if(i*1000000>NROW(bedA)) NROW(bedA) else i*1000000;
         y <- try( unlist(bed.region.bpQuery.bigWig( bw=bw.hist, bed= bedA[start:end,]) ) );
         if(class(y)=="try-error")
         {
            stop("error");
            return(NULL);
         }
         unload.bigWig(bw.hist);
         return(y);
         }, mc.cores=5) );
   }
   
   snow.task.size <- 1000000
   fun.snow<-function(i)
   {
      library(bigWig)
      bw.hist <- load.bigWig(file.histone);
       start <- 1+(i-1)*snow.task.size;
       end <- if(i*snow.task.size>NROW(bedA)) NROW(bedA) else i*snow.task.size;
      y <- try( unlist(bed.region.bpQuery.bigWig( bw=bw.hist, bed= bedA[start:end,]) ) );
      if(class(y)=="try-error")
      {
           stop("error");
         return(NULL);
        }
      unload.bigWig(bw.hist);
      return(y);
   } 

   if(ncores>1)
   {
     library(snowfall);
     sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
     sfExport("bedA", "file.histone","snow.task.size" );

     fun <- as.function(fun.snow);
     environment(fun)<-globalenv();

     y <- unlist(sfLapply( 1:ceiling(NROW(bedA)/snow.task.size), fun));
     sfStop();
   }
   else
     y <- unlist(lapply( 1:ceiling(NROW(bedA)/snow.task.size), fun.snow));
     
   if(!is.null(ref.bed))
   {
      file.map <- tempfile(fileext=".bed")
      write.table( cbind(bedA, y), file=file.map, quote=F, row.names=F, col.names=F, sep="\t");
      file.ref <- tempfile(fileext=".bed")
      write.table( ref.bed, file=file.map, quote=F, row.names=F, col.names=F, sep="\t");
      tb <- read.table(pipe( paste("bedmap --echo --mean ", file.ref,  file.map )));
      y <- tb[,4];

      if(normalized)
         y <- y / (bw.hist$basesCovered * bw.hist$mean ) * 1e9
   }

   gc();
   
#cat("Y size=", NROW(y), "\n");    
   return(y);
}

draw_trend <- function( out.prefix, file.region.bed, file.exp.bw, file.pred.bw, ref.gene.bed=NULL, scatter=TRUE, scatter.log=NULL, win.size=2000, est.range=c(10, 20, 50, 100, seq(200, 10000, 200)) )
{
    cat("source 1 file=", file.bed, "\n");
    cat("source 2 file=", file.bw1, "\n");
    cat("source 3 file=", file.bw2, "\n");

    tb.s1 <- read.table(file.bed);
    if(is.null(file.bw2))
    {
        y.exp <- tb.s1[,4];
        y.bw1 <- get_histone_read( tb.s1[, c(1:3)], file.bw1, ref.gene.bed);
    }
    else
    {
        y.bw1 <- get_histone_read( tb.s1[, c(1:3)], file.bw1, ref.gene.bed);
        y.exp <- get_histone_read( tb.s1[, c(1:3)], file.bw2, ref.gene.bed);
    }

    if(scatter)
    {
        r <- rbindlist( mclapply(est.range, function(x) { calculate_cor( y.exp, y.bw1, x) }, mc.cores=5  ) )

        pdf(paste(out.prefix, "trend", "pdf", sep="."), width=6, height=6);
        plot(est.range, unlist(r[,1]), type="l", col="red", cex=0.6, xlim=range(est.range), ylim=c(-0.2,1) )
        lines(est.range, unlist(r[,2]), cex=0.6, col="blue")
        lines(est.range, unlist(r[,3]/max(r[,3])), cex=0.6, col="black", lty="22");
        
        legend("topleft", c("pearson", "spearman", "mad"), col=c("red", "blue", "black"), lty=c("solid", "solid", "22") );
        dev.off();
    }

}


dens_scatter_plot <- function(x1, x2, uselog=FALSE, n=256, cex=2, ...) {

  x1 <- x1+1;
  x2 <- x2+1;
  df <- data.frame(x1, x2)

  if(uselog) {
    x <- densCols(log(x1,10),log(x2,10), colramp=colorRampPalette(c("black", "white")))
  } else {
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  }

  # log or non-log style use same color scheme
  #x <- densCols(log(x1,10),log(x2,10), colramp=colorRampPalette(c("black", "white")))
  
  df$dens <- col2rgb(x)[1,] + 1L
  ## Map densities to colors
  cols <- colorRampPalette(c("light gray", "#000099", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(48)
  cols <- c(cols, rep( "#FF1F00", n-48) );
  df$col <- cols[df$dens]
   
  par(mar=par("mar")*1.3); 
  par.mgp <- par("mgp")
  par.mgp[1] <- par.mgp[1]*1.2;
  par(mgp=par.mgp); 
  
  ## Plot it, reordering rows so that densest points are plotted on top
  if(uselog) {
    plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=cex, log="xy", ...)
    args <- list(...)
    if(is.null(args[['xlim']]))
    {
      usr <- par("usr");
      par(new=TRUE)                                                                                                                                                                                                                                                              
      plot(1,1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(usr[1], usr[2]), ylim=c(usr[3], usr[4]));
    }
    #else
    #{ 
    #  par(new=TRUE)                                                                                                                                                                                                                                                              
    #  plot(1,1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=args[['xlim']], ylim=args[['ylim']]);
    #}
    segments(-2,-2, 100000, 100000, col="black");
  } else {
    plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=cex, ...)
    segments(0,0, 1000000, 1000000, col="black");
  }
  
}

draw_cor <- function( out.prefix, file.pred.bw, file.exp.bw, 
            file.ctrl.bw = NULL,
            file.exp.peak = NULL, 
            file.unmap.bed = NULL,  
            file.black.bed = NULL, 
            file.cell.remove = NULL,
            chr = NULL, 
            gplot = TRUE, 
            overwrite=TRUE,
            save.rdata=FALSE,
            ref.gene.bed = NULL, 
            scatter.log = NULL, 
            win.size = 20000, 
            xlab = "Histone", ylab = "Predict" , xlim = NULL, ylim = NULL, 
            file.draw.peak = NULL, peak.ext = 100)
{
    cat("param 0 chr = ", chr, "\n");
    cat("param 1 win.size = ", win.size, "\n");
    cat("param 2 prefix = ", out.prefix, "\n");
    cat("param 3 file (file.exp.bw)=", file.exp.bw, "\n");
    cat("param 4 file (file.pred.bw)=", file.pred.bw, "\n");
    cat("param 5 file (file.ctrl.bw)=", file.ctrl.bw, "\n");
    cat("param 6 file (file.exp.peaks)=", file.exp.peak, "\n");
    cat("param 7 file (file.draw.peak)=", file.draw.peak, "\n");

    if (is.na(file.exp.bw) || is.na(file.pred.bw))
       return(NA);

    fig.prefix <- if (is.null(chr)) paste(out.prefix, "all", sep=".") else  paste(out.prefix, chr, sep="."); 
    if (!is.null(file.draw.peak)) fig.prefix <- paste(fig.prefix, "peak", sep=".")
    
    file.png <- paste(fig.prefix, paste0("w", round(win.size/1000,3), "k"), sep=".") 
    if (!is.null(file.draw.peak)) file.png <- paste(file.png, paste0("ex", peak.ext), sep=".") 

    file.log.png <- paste(file.png, "log.png", sep=".")
    file.nonlog.png <- paste(file.png, "png", sep=".")
    if (gplot && !overwrite && file.exists(file.log.png))
    {
       cat(file.log.png, "is existing.\n");
       return(NA);
    }   

    ## get the huge-mapped regions
    file.tmp.hugemap <- NULL; 
    if(is.null(file.exp.peak) || is.na(file.exp.peak) )
       file.tmp.hugemap <- detect_hugemap_reg(file.exp.peak, file.exp.bw)
    
    ## merge all remvoved regions
    file.removed.reg <- NULL;
    if(!( is.null(file.tmp.hugemap) && is.null(file.unmap.bed) && is.null(file.black.bed)))
       file.removed.reg <-  merge_removed_regs( file.unmap.bed, file.black.bed, file.cell.remove, file.tmp.hugemap )

    ## get all window regions based on window size and chromosome info.
    tb.chr <- get_chrom_info( file.exp.bw )
    tb.s1 <- as.data.frame(rbindlist(apply(tb.chr, 1, function(x){ 
           seq.frag=1:floor(as.numeric(x[2])/win.size); 
           bed <- data.frame( chr=as.character(x[1]), start=(seq.frag-1)*win.size+1, stop=seq.frag*win.size);
           bed <- bed [ bed$start>=100000 & bed$stop <= (as.numeric(x[2])-100000),];
           # only for test Alex data...
           # bed <- bed [ bed$start>=23,830,526,]; #size of chr. 22 51,304,566
           bed; } )));

    ## only keep chr1-22 
    exclude_chromosome = "_|chrM|chrY|chrX|rRNA|chr13_random"; 
    tb.s1 <- tb.s1 [ grep( exclude_chromosome , tb.s1[,1], invert=TRUE), ]

    ## only keep specified chromosome. 
    if(!is.null(chr))
       tb.s1 <- tb.s1 [ tb.s1[,1]==chr,] 
cat("window count1 =", NROW(tb.s1), "\n");
    
    if(NROW(tb.s1)==0)
       stop("No windows!");

    ## substract the removed regions.
    tbs0 <-  data.frame(tb.s1, idx=1:NROW(tb.s1))
    if(!is.null(file.removed.reg))
    {
       file.tb.tmp <- write.temp.bed( tbs0 );
       tb.s1 <- read.table(pipe(paste("bedtools subtract  -a ", file.tb.tmp, "-b", file.removed.reg )))
       unlink(file.tb.tmp);
    }
    else
       tb.s1 <- data.frame(tb.s1, idx=1:NROW(tb.s1))
    unlink(file.removed.reg);
    unlink(file.tmp.hugemap);
cat("window count2 =", NROW(tb.s1), "\n")    

    ## get the intersect with the peak regions if only peak regions are draws.
    if(!is.null(file.draw.peak))
    {
       tb.peaks <- read.table( file.draw.peak ); 
       tb.peaks[,2] <- tb.peaks[,2] - peak.ext;
       if (NROW(which(tb.peaks[,2]<0)) >0 )
          tb.peaks[which(tb.peaks[,2]<0),2] <- 0
       tb.peaks[,3] <- tb.peaks[,3] + peak.ext;
          
       file.tb.tmp<- write.temp.bed( tb.s1 );
       file.tb.peaks <- write.temp.bed( tb.peaks );
       tb.s1 <- read.table(pipe(paste("bedtools intersect -a ", file.tb.tmp, "-b", file.tb.peaks, "-wa"))); 
       unlink(file.tb.tmp);
       unlink(file.tb.peaks);
       
       cat("window count3 =", NROW(tb.s1), "\n")    
    }


    y.win.exp0 <- y.win.exp  <- get_histone_readX( tb.s1[, c(1:3)], file.exp.bw, ref.gene.bed)/10;
    y.win.pred <- get_histone_readX( tb.s1[, c(1:3)], file.pred.bw, ref.gene.bed);

    ## get histone reads from experiment data and imputed data. 
    tby <- data.frame(ID=tb.s1[,4], exp=y.win.exp, pred=y.win.pred);
    tby0 <- sqldf("select ID, sum(exp) as exp, sum(pred) as pred  from tby group by ID order by ID");
    bigtb1 <- sqldf("select chr, start, stop, tbs0.idx, exp, pred from tbs0, tby0 where tbs0.idx==tby0.ID order by tbs0.idx")
   
    # ex. bigtb1 
    #     chr  start   stop idx exp pred
    # 1 chr22 100001 110000   1   0    0
    # 2 chr22 110001 120000   2   0    0
    # 3 chr22 120001 130000   3   0    0
    # 4 chr22 130001 140000   4   0    0
    # 5 chr22 140001 150000   5   0    0

    ## the experiment data for each region
    y.win.exp <- bigtb1$exp;
    ## the imputed data for each region
    y.win.pred <- bigtb1$pred;
    
    bigtb1$ratio <- (bigtb1$exp+1)/(bigtb1$pred+1)
    ## for test
    ## show( bigtb1[y.win.exp>20000 & bigtb1$ratio>4,]);

    ##2019/09/19 Remove this filter   
    idx.rem <- c();
    if(0)
    {
    if(!is.null(file.ctrl.bw) && !is.na(file.ctrl.bw) )
    {
        y.win.ctrl  <- get_histone_readX( tb.s1[, c(1:3)], file.ctrl.bw, ref.gene.bed)/10;
        tby <- data.frame(ID=tb.s1[,4], exp=y.win.exp0, ctrl=y.win.ctrl);
        tby0 <- sqldf("select ID, sum(exp) as exp, sum(ctrl) as ctrl from tby group by ID order by ID");
        bigtb0 <- sqldf("select chr, start, stop, tbs0.idx, exp, ctrl, (exp+1)/(ctrl+1) as ratio from tbs0, tby0 where tbs0.idx==tby0.ID order by tbs0.idx")

        # ex. bigtb0 
        #     chr  start   stop idx org ctrl ratio
        # 1 chr22 100001 110000   1   0    0     1
        # 2 chr22 110001 120000   2   0    0     1
        # 3 chr22 120001 130000   3   0    0     1
        # 4 chr22 130001 140000   4   0    0     1
        # 5 chr22 140001 150000   5   0    0     1
        # 6 chr22 150001 160000   6   0    0     1
        
        q1=0.99 # 0.05 by default, other 0.99
        q2=0.99;
        tbdel <- bigtb0[ bigtb0$exp>quantile(bigtb0$exp, q1) & bigtb0$ratio>quantile(bigtb0$ratio, q2),]
        #print( tbdel );

        idx.rem <- which( (bigtb0$exp>quantile(bigtb0$exp, q1) & bigtb0$ratio>quantile(bigtb0$ratio, q2)) );

        if(NROW(idx.rem)>0)
        {
           #print(idx.rem);    
           #show( bigtb1[idx.rem,]);
           y.win.exp <- bigtb1$exp[-idx.rem];
           y.win.pred <- bigtb1$pred[-idx.rem];
        }
        else
        {
           y.win.exp <- bigtb1$exp
           y.win.pred <- bigtb1$pred
        }   
    }
    }

cat("window count 4=", NROW(y.win.pred), "\n")    
    r.cor <- comp_cor_coefficient (y.win.exp, y.win.pred)
    str.cor <- paste("COR=", round(r.cor[1],3),  round(r.cor[2], 3), round(r.cor[3], 3), "JSD=", round(r.cor[4] ,3) );
    cat(str.cor, "\n\n");
 
    if(save.rdata)
       save( y.win.exp, y.win.pred, file=paste0(file.nonlog.png, ".rdata"));
    if(gplot)
    {
        png(file.nonlog.png, width=900, height=900)
        if(is.null(xlim)) xlim <- c(min(y.win.exp), quantile(y.win.exp, 0.999) );
        if(is.null(ylim)) ylim <- c(min(y.win.pred), quantile(y.win.pred, 0.999) );
        dens_scatter_plot ( y.win.exp, y.win.pred, cex=1.5, cex.axis=2, cex.lab=2, main=paste(fig.prefix, "(", str.cor, ")"), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
        dev.off();
cat("PNG 1=", file.nonlog.png, "\n");    

        png(file.log.png, width=900, height=900)
        dens_scatter_plot ( y.win.exp, y.win.pred, uselog=TRUE, cex=1.5, cex.axis=2, cex.lab=2, main=paste(fig.prefix, "(", str.cor, ")"), xlab=paste0("log(", xlab, ")"), ylab=paste0("log(", ylab, ")"));
        dev.off();
cat("PNG 2=", file.log.png, "\n");    
    }

    return(data.frame(plot=out.prefix, rem.count=NROW(idx.rem), pearson=r.cor[1], spearman=r.cor[2], mad=r.cor[3], JSD=r.cor[4], win.size=win.size, chr=if(is.null(chr))"all" else chr, peak.ext = peak.ext))
}

comp_cor_coefficient <- function(y.exp, y.pred, exclZeros=FALSE)
{
    if(exclZeros) {
      indx <- y.exp > 0 & y.pred > 0
    }
    else{ 
      indx <- rep(TRUE, NROW(y.exp))
    }
    r.cor1 = cor( y.exp[indx], y.pred[indx], method = "pearson" );
    r.cor2 = cor( y.exp[indx], y.pred[indx], method = "spearman" );
    r.mad = mad( y.exp -  y.pred );

    r.gmic = 0.0;
    #if(NROW(y.exp)<35000)
    #  r.gmic = mine( y.exp,  y.pred, n.cores=12, )$MIC;

    x <- rbind(y.exp/sum(y.exp), y.pred/sum(y.pred))
    r.JSD <- distance(x, method = "jensen-shannon");

    return(c(pearson=r.cor1, spearman=r.cor2, mad=r.mad, JSD=r.JSD))
}

detect_hugemap_reg<-function(file.exp.peak, file.exp.bw)
{
   peak_max <- get_histone_peak_max( file.exp.peak, file.exp.bw );
   y_max <- quantile(peak_max[,4],0.98);
   peak_spike <- peak_max[peak_max[,4] > y_max,] ;

   bw.hist <- load.bigWig(file.exp.bw);
   max_count <- unlist(lapply(1:NROW(peak_spike), function(i) {
         bedv <- step.bpQuery.bigWig( bw.hist, as.character(peak_spike[i,1]),peak_spike[i,2], peak_spike[i,3],1)
         bedv_max <- sort(abs(bedv[-NROW(bedv)] - bedv[-1]), decreasing=T);
         return( (bedv_max[1]+bedv_max[2])/max(bedv)  );
   }));
   
   unload.bigWig( bw.hist)

   bed.missmap <- data.frame(peak_spike, max_count);
   file.temp  <- write.temp.bed( bed.missmap[bed.missmap[,5]>1, 1:3 ], compress=FALSE )
   
   return(file.temp);
}

merge_removed_regs <- function( file.unmap.bed, file.black.bed, file.cell.remove, file.exp.hugemap )
{
   file.temp.black  <- write.temp.bed(read.table(file.black.bed)[,c(1:3)], compress=FALSE )
   
   file.temp.unmap <-"";
   if(!is.null(file.unmap.bed) && !is.na(file.unmap.bed) )
   {
       tb.unmap.bed <- read.table(file.unmap.bed)[,c(1:3)];
       tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
       file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
   }

   bed.regs <- read.table(pipe(paste("cat ", file.temp.black, file.temp.unmap, file.cell.remove, file.exp.hugemap, " | awk '{ if ($1!=\"rRNA\") print $0 }' - | sort-bed - | bedtools merge -i - ")))
   file.removed.reg <- write.temp.bed(bed.regs);
   unlink(file.temp.black);
   unlink(file.temp.unmap);

   return(file.removed.reg);
}


