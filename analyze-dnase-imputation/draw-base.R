
library(sqldf);
library(bigWig)
library(dREG)
library(data.table);
library(snowfall);
library(minerva);
library(parallel);
library(philentropy)

#source("/workdir/zw355/proj/prj15-histone/scripts/hist.svm.main.R");
source("/workdir/zw355/proj/prj15-histone/scripts/hist.svm.com.R");
source("/workdir/zw355/proj/prj15-histone/scripts/hist.param.R");

calculate_cor <- function (y.org, y.pred , range=2000 )
{
   if(range!=10)
   {
	   n.per.range <- round( range/10 );

	   y.range.org <-  unlist(lapply(1:floor(NROW(y.org)/n.per.range),  function(i) {sum(y.org [(i-1)*n.per.range+c(1:n.per.range)])}));
	   y.range.pred <- unlist(lapply(1:floor(NROW(y.pred)/n.per.range), function(i) {sum(y.pred[(i-1)*n.per.range+c(1:n.per.range)])}));
   }
   else
   {
		y.range.org <- y.org;
		y.range.pred <- y.pred;
   }

   return(list(r.pearson = cor( y.range.org, y.range.pred, method = "pearson" ),
      r.spearman = cor( y.range.org, y.range.pred, method = "spearman" ),
      r.mad = mad( y.range.org -  y.range.pred ) )  );
}

get_histone_peak<-function( file.peak, file.histone)
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

get_histone_readX<-function(bedA, file.histone, ref.bed=NULL, method="mean", normalized=TRUE)
{
   bedA <- bedA[,c(1:3)];
   na.idx <- which( is.na(bedA[,1]) |  is.na(bedA[,2]) |  is.na(bedA[,3])  )
   if(length(na.idx)>0)
   	  stop("NA in bedA\n");

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

cat("Y size=", NROW(y), "\n");    
   return(y);
}

draw_trend <- function( out.prefix, file.region.bed, file.org.bw, file.pred.bw, ref.gene.bed=NULL, scatter=TRUE, scatter.log=NULL, range=2000, est.range=c(10, 20, 50, 100, seq(200, 10000, 200)) )
{
	cat("source 1 file=", file.bed, "\n");
	cat("source 2 file=", file.bw1, "\n");
	cat("source 3 file=", file.bw2, "\n");

	tb.s1 <- read.table(file.bed);
	if(is.null(file.bw2))
	{
		y.org <- tb.s1[,4];
		y.bw1 <- get_histone_read( tb.s1[, c(1:3)], file.bw1, ref.gene.bed);
	}
	else
	{
		y.bw1 <- get_histone_read( tb.s1[, c(1:3)], file.bw1, ref.gene.bed);
		y.org <- get_histone_read( tb.s1[, c(1:3)], file.bw2, ref.gene.bed);
	}

	if(scatter)
	{
	    r <- rbindlist( mclapply(est.range, function(x) { calculate_cor( y.org, y.bw1, x) }, mc.cores=5  ) )

		pdf(paste(out.prefix, "trend", "pdf", sep="."), width=6, height=6);
		plot(est.range, unlist(r[,1]), type="l", col="red", cex=0.6, xlim=range(est.range), ylim=c(-0.2,1) )
		lines(est.range, unlist(r[,2]), cex=0.6, col="blue")
		lines(est.range, unlist(r[,3]/max(r[,3])), cex=0.6, col="black", lty="22");
		
		legend("topleft", c("pearson", "spearman", "mad"), col=c("red", "blue", "black"), lty=c("solid", "solid", "22") );
		dev.off();
	}

}


get.chrom.info <- function( file.bw.plus )
{
    bw.plus <- load.bigWig(file.bw.plus);

    chrom <-  cbind( bw.plus$chroms, bw.plus$chromSizes)
    chr.size <- unlist( lapply( unique(chrom[,1]), function(chr){max( as.numeric( chrom[which(chrom[,1]==chr),2])) } ) );

    df.bed <- data.frame( V1=unique(chrom[,1]), V2=chr.size );
    return(df.bed);
}

densScatterplot <- function(x1, x2, uselog=FALSE, n=256, cex=2, ...) {

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


draw_cor_signal <- function( out.prefix, file.org.bw, file.pred.bw, file.ctrl.bw=NULL, chr=NULL, ref.gene.bed=NULL, scatter.log=NULL, range=20000, xlab="Histone", ylab="Predict" , xlim=c(0,50000), ylim=c(0, 50000), file.peaks=NULL, peak.mode=0, peak.ext=100, out.bed.file=NULL)
{
	cat("source 1 file=", out.prefix, "\n");
	cat("source 2 file (file.org.bw)=", file.org.bw, "\n");
	cat("source 3 file (file.pred.bw)=", file.pred.bw, "\n");
	cat("source 4 file (file.ctrl.bw)=", file.ctrl.bw, "\n");
	cat("source 5 file (file.peaks)=", file.peaks, "\n");
	cat("writing data to=", out.bed.file, "\n");

    if (is.na(file.org.bw) || is.na(file.pred.bw))
       return(NA);
#browser();
	file.png <- if (is.null(chr))  
	       paste(out.prefix, "all", "png", sep=".")
	   else  
	       paste(out.prefix, chr, "png", sep=".");

    if (!g_noplot && file.exists(file.png))
    {
       cat(file.png , "is existing.\n");
       return(NA);
    }   

    tb.chr <- get.chrom.info( file.org.bw )
    tb.s1 <- do.call("rbind", apply(tb.chr, 1, function(x){ 
           seq.frag=1:floor(as.numeric(x[2])/range); 
           bed <- data.frame( chr=as.character(x[1]), start=(seq.frag-1)*range+1, stop=seq.frag*range);
           bed <- bed [ bed$start>=100000 & bed$stop <= (as.numeric(x[2])-100000),];
           # for HAlex data...
           # bed <- bed [ bed$start>=23,830,526,]; #size of chr. 22 51,304,566
           bed; } ));

    exclude_chromosome = "_|chrM|chrY|chrX|rRNA"; 
    tb.s1 <- tb.s1 [ grep( exclude_chromosome , tb.s1[,1], invert=TRUE), ]


    tbs0 <-  data.frame(tb.s1, idx=1:NROW(tb.s1))
    file.tb.s1 <- write.temp.bed( tbs0 );
    tb.s1 <- read.table(pipe(paste("bedtools subtract  -a ", file.tb.s1, "-b", file.removed.reg )));
    unlink(file.tb.s1);
cat("regions 1=", NROW(tb.s1), "\n")    

    if(!is.null(file.peaks))
    {
       if (peak.mode==0)
       {  
          file.tb.s1 <- write.temp.bed( tb.s1 );
          tb.s1 <- read.table(pipe(paste("bedtools intersect -a ", file.tb.s1, "-b", file.peaks, "-wa"))); 
       }
       else
       {  
          tb.peaks <- read.table( file.peaks ); 
          tb.peaks[,2] <- tb.peaks[,2] - peak.ext;
          if (NROW(which(tb.peaks[,2]<0)) >0 )
             tb.peaks[which(tb.peaks[,2]<0),2] <- 0
          tb.peaks[,3] <- tb.peaks[,3] + peak.ext;
          
          file.tb.s1 <- write.temp.bed( tb.s1 );
          file.tb.peaks <- write.temp.bed( tb.peaks );
          tb.s1 <- read.table(pipe(paste("bedtools intersect -a ", file.tb.s1, "-b", file.tb.peaks, "-wa"))); 
       }
    }

    if(!is.null(chr))
       tb.s1 <- tb.s1 [ tb.s1[,1]==chr,] 

    exclude_chromosome = "_|chrM|chrY|chrX|rRNA"; 
    tb.s1 <- tb.s1 [ grep( exclude_chromosome , tb.s1[,1], invert=TRUE), ]

cat("regions 2=", NROW(tb.s1), "\n")    

    y.range.org0 <- y.range.org  <- get_histone_readX( tb.s1[, c(1:3)], file.org.bw, ref.gene.bed)/10;
    y.range.pred <- get_histone_readX( tb.s1[, c(1:3)], file.pred.bw, ref.gene.bed);

    tby <- data.frame(ID=tb.s1[,4], org=y.range.org, pred=y.range.pred);
    tby0 <- sqldf("select ID, sum(org) as org, sum(pred) as pred  from tby group by ID order by ID");
    bigtb1 <- sqldf("select chr, start, stop, tbs0.idx, org, pred from tbs0, tby0 where tbs0.idx==tby0.ID order by tbs0.idx")
    y.range.org <- bigtb1$org;
    y.range.pred <- bigtb1$pred;
    
    bigtb1$ratio <- (bigtb1$org+1)/(bigtb1$pred+1)
#show( bigtb1[y.range.org>20000 & bigtb1$ratio>4,]);

    idx.rem <- c();
    if(!is.null(file.ctrl.bw) && !is.na(file.ctrl.bw) )
    {
		y.range.ctrl  <- get_histone_readX( tb.s1[, c(1:3)], file.ctrl.bw, ref.gene.bed)/10;
		tby <- data.frame(ID=tb.s1[,4], org=y.range.org0, ctrl=y.range.ctrl);
        tby0 <- sqldf("select ID, sum(org) as org, sum(ctrl) as ctrl from tby group by ID order by ID");
        bigtb0 <- sqldf("select chr, start, stop, tbs0.idx, org, ctrl, (org+1)/(ctrl+1) as ratio from tbs0, tby0 where tbs0.idx==tby0.ID order by tbs0.idx")
        
        q1=0.99 # 0.05 by default, other 0.99
        q2=0.99;
        tbdel <- bigtb0[ bigtb0$org>quantile(bigtb0$org, q1) & bigtb0$ratio>quantile(bigtb0$ratio, q2),]
	    #print( tbdel );

	    idx.rem <- which( (bigtb0$org>quantile(bigtb0$org, q1) & bigtb0$ratio>quantile(bigtb0$ratio, q2)) );

        if(NROW(idx.rem)>0)
        {
           #print(idx.rem);	
           #show( bigtb1[idx.rem,]);
           y.range.org <- bigtb1$org[-idx.rem];
           y.range.pred <- bigtb1$pred[-idx.rem];
        }
        else
        {
           y.range.org <- bigtb1$org
           y.range.pred <- bigtb1$pred
        }   

	}

cat("regions 3=", NROW(y.range.pred), "\n")    
cat("q(pred, 0.25)", quantile(y.range.pred, 0.25), "\n")    

	r.cor1 = cor( y.range.org, y.range.pred, method = "pearson" );
	r.cor2 = cor( y.range.org, y.range.pred, method = "spearman" );
    r.mad = mad( y.range.org -  y.range.pred );

    r.gmic = 0.0;
    #if(NROW(y.range.org)<35000)
    #  r.gmic = mine( y.range.org,  y.range.pred, n.cores=12, )$MIC;

    x <- rbind(y.range.org/sum(y.range.org), y.range.pred/sum(y.range.pred))
    r.JSD <- distance(x, method = "jensen-shannon");

    str.cor <- paste("COR=", round(r.cor1,3),  round(r.cor2,3), round(r.mad,3), "JSD=", round(r.JSD ,3) );

    if(!g_noplot)
    {
		source("/home/zw355/src/Rplot/denScatter.R");
		densScatterplot( y.range.org, y.range.pred, cex=1.5, cex.axis=2, cex.lab=2, main=paste(out.prefix, "(", str.cor, ")"), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);

		png(file.png, width=900, height=900)
		densScatterplot( y.range.org, y.range.pred, cex=1.5, cex.axis=2, cex.lab=2, main=paste(out.prefix, "(", str.cor, ")"), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
		dev.off();

#		if (is.null(chr))
#			png(paste(out.prefix, "all", "log", "png", sep="."), width=900, height=900)
#		else
#			png(paste(out.prefix, chr, "log", "png", sep="."), width=900, height=900);

#		densScatterplot( y.range.org, y.range.pred, uselog=FALSE, cex=1.5, cex.axis=2, cex.lab=2, main=paste(out.prefix, "(", str.cor, ")"), xlab=paste0("log(", xlab, ")"), ylab=paste0("log(", ylab, ")"))#, xlim=c(0.01, log(max(xlim), 10)), ylim=c(0.01, log(max(ylim),10)));
#		dev.off();
	}

  if(!is.null(out.bed.file)) {
    ## Generate a data.frame with the coordinates, experimental, and predicted values used for the analysis.
    print(paste("Writing out BED file.", NROW(tb.s1), NROW(y.range.org), NROW(y.range.pred), NROW(bigtb1)))
    print(head(bigtb1))
    #dataBedToWrite <- cbind(tb.s1, experiment= y.range.org, imputed= y.range.pred);
    write.table(bigtb1, file=out.bed.file, quote=F, row.names=F, col.names=F, sep="\t");

  }
		
	return(data.frame(plot=out.prefix, rem.count=NROW(idx.rem), pearson=r.cor1, spearman=r.cor2, mad=r.mad, JSD=r.JSD))
	
}

detec_spike_reg<-function(file.org.peak, file.org.bw)
{
   peak_max <- get_histone_peak( file.org.peak, file.org.bw );
   y_max <- quantile(peak_max[,4],0.98);
   peak_spike <- peak_max[peak_max[,4] > y_max,] ;

   bw.hist <- load.bigWig(file.org.bw);
   max_count <- unlist(lapply(1:NROW(peak_spike), function(i) {
   	  bedv <- step.bpQuery.bigWig( bw.hist, as.character(peak_spike[i,1]),peak_spike[i,2], peak_spike[i,3],1)
   	  bedv_max <- sort(abs(bedv[-NROW(bedv)] - bedv[-1]), decreasing=T);
   	  return( (bedv_max[1]+bedv_max[2])/max(bedv)  );
   }));
   
   unload.bigWig( bw.hist)
   
   return(data.frame(peak_spike, max_count));
}

add_spikes_to_removed_regs <- function(file.peak, file.bw)
{
	bed_spike <- detec_spike_reg(file.peak, file.bw)
	file.temp.spike  <- write.temp.bed( bed_spike[bed_spike[,5]>1,1:3], compress=FALSE )
	file.removed.reg <<- write.temp.bed(read.table(pipe(paste("cat ", file.temp.black, file.temp.unmap , file.temp.spike, " | sort-bed - | bedtools merge -i - "))));
}


