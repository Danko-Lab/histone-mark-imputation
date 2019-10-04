source("draw-cor-base.R");
source("draw-cor-param.R");

if(0)
{
    NCORES <- 5;

    df.GM.chr22 <- c();
    for(zoom_rate in c(10, 1, 0.1))
    {
       L=mclapply( 1:NROW(df.GM.info), function(i) {
	  rx <- draw_cor( paste( df.GM.info[i,"Cell"], df.GM.info[i,"Mark"], df.GM.info[i,"Ver"], sep="." ),   
			     df.GM.info[i,"file.pred.bw"], 
			     df.GM.info[i,"file.exp.bw"], 
			     file.ctrl.bw = df.GM.info[i,"file.ctrl.bw"],  
			     file.exp.peak = df.GM.info[i,"file.exp.peak"],  
			     file.unmap.bed = df.GM.info[i,"file.unmap.bed"],  
			     file.black.bed = df.GM.info[i,"file.black.bed"],  
			     file.cell.remove = NULL,
			     win.size = 1000*zoom_rate,
			     chr = "chr22", 
			     gplot = TRUE, 
			     overwrite=TRUE,
			     save.rdata=TRUE);}, mc.cores=NCORES  )
	df.GM.chr22 <- rbind(df.GM.chr22, data.frame( do.call("rbind", L), win=1000*zoom_rate ) )
    }
}


dens_scatter_plot2 <- function(x1, x2, uselog=FALSE, n=256, cex=2, ...) {

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

draw_scatter_subplot<-function(fig.rdata, log=TRUE)
{
    #y.win.exp, y.win.pred
    load(file=fig.rdata);
   
    if(!log)
    {
        xlim <- c(min(y.win.exp), quantile(y.win.exp, 0.999) );
        ylim <- c(min(y.win.pred), quantile(y.win.pred, 0.999) );
        dens_scatter_plot2 ( y.win.exp, y.win.pred, 
                            cex=1.5, cex.lab=2, cex.axis=2, 
                            main="", xlab="", ylab="", xlim=xlim, ylim=ylim,xaxt="n", yaxt="n");
    }
    else
    {
        dens_scatter_plot2 ( y.win.exp, y.win.pred, uselog=TRUE, 
                            cex=1.5, cex.axis=2, cex.lab=2, 
                            main="", xlab="", ylab="", xaxt="n", yaxt="n" );
    }
}

if(0)
{
    bLog <- TRUE;

    if(bLog)
       png("draw-gm-scatter-log.png", width=3.5*900, height=9.5*900)
    else    
       png("draw-gm-scatter-nonlog.png", width=3.5*900, height=9.5*900)

    lay.heights <- c( rep(1,9), 0.5);
    lay.widths  <- c( 0.5, rep(1,3) )
    layout(matrix(c(1:40), ncol=4, byrow=F), widths=lay.widths, heights=lay.heights)

    for(i in 1:NROW(df.GM.info) ) 
    {
	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
	plot(1, type="n", xaxt = 'n', yaxt = 'n', bty = 'n', xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
	text(8, 5, df.GM.info[i,"Mark"], cex=10, srt=90)
    }

    plot(1, type="n", xaxt = 'n', yaxt = 'n', bty = 'n', xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))

    for(zoom_rate in c(0.1, 1, 10))
    {
	for(i in 1:NROW(df.GM.info))
       {
	   fig.rdata <- paste( df.GM.info[i,"Cell"], df.GM.info[i,"Mark"], df.GM.info[i,"Ver"], "chr22", paste0("w", zoom_rate , "k"), "png.rdata",sep=".")
	   draw_scatter_subplot(fig.rdata, bLog);
       }
       plot(1, type="n", xaxt = 'n', yaxt = 'n', bty = 'n', xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
       text(5, 8, zoom_rate* 1000, cex=10)
    }   
    dev.off();
}

if(1)
{
    bLog <- TRUE;
    if(bLog)
       png("draw-gm-scatter-log-3x3.png", width=(0.5*4+3*2)*900, height=(0.5+3)*900)

    lay.heights <- c( 1, 1, 1, 0.5);
    lay.widths  <- c( 0.25, 1, 1, 0.25, 1, 1, 0.25, 1, 1, 0.25 )
    layout(matrix(c(1, 10, 11,  4, 16, 17,  7, 22, 23, 35,
	     2, 12, 13,  5, 18, 19,  8, 24, 25, 35,
	     3, 14, 15,  6, 20, 21,  9, 26, 27, 35,
	    34, 28, 29, 34, 30, 31, 34, 32, 33, 35), nrow=4, byrow=T), widths=lay.widths, heights=lay.heights)

    for(i in 1:NROW(df.GM.info) ) 
    {
	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
	plot(1, type="n", xaxt = 'n', yaxt = 'n', bty = 'n', xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
	text(8, 5, df.GM.info[i,"Mark"], cex=10, srt=90)
    }

    for(i in 1:NROW(df.GM.info))
    {
      for(zoom_rate in c(1, 10))
      {
	   fig.rdata <- paste( df.GM.info[i,"Cell"], df.GM.info[i,"Mark"], df.GM.info[i,"Ver"], "chr22", paste0("w", zoom_rate , "k"), "png.rdata",sep=".")
	   draw_scatter_subplot(fig.rdata, bLog);
       }
    }   

    for(i in 1:3)
    for(zoom_rate in c(1, 10))
    {
       plot(1, type="n", xaxt = 'n', yaxt = 'n', bty = 'n', xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
       text(5, 8, zoom_rate* 1000, cex=10)
    }

    plot(1, type="n", xaxt = 'n', yaxt = 'n', bty = 'n', xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
    dev.off();
}
